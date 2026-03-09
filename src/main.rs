// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

extern crate clap;
extern crate csv;
extern crate failure;
extern crate rayon;
extern crate rust_htslib;
extern crate simplelog;
extern crate tempfile;
extern crate terminal_size;
#[macro_use]
extern crate log;
extern crate human_panic;

use clap::{App, Arg};
use failure::Error;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use simplelog::*;
use std::cmp;
use std::collections::HashSet;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::{Path, PathBuf};
use std::process;
use tempfile::tempdir;
use terminal_size::{terminal_size, Width};

fn get_args() -> clap::App<'static, 'static> {
    let args = App::new("subset-bam")
        .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
        .version("1.1.0")
        .author("Ian Fiddes <ian.fiddes@10xgenomics.com>, Wyatt McDonnell <wyatt.mcdonnell@10xgenomics.com>")
        .about("Subsetting 10x Genomics BAM files")
        .arg(Arg::with_name("bam")
             .short("b")
             .long("bam")
             .value_name("FILE")
             .multiple(true)
             .help("Cellranger BAM/CRAM file. Repeat this argument to process multiple files.")
             .required(true))
        .arg(Arg::with_name("cell_barcodes")
             .short("c")
             .long("cell-barcodes")
             .value_name("FILE")
             .help("File with cell barcodes to be extracted.")
             .required(true))
        .arg(Arg::with_name("out_bam")
             .short("o")
             .long("out-bam")
             .value_name("OUTPUT_FILE")
             .help("Output BAM.")
             .required(true))
        .arg(Arg::with_name("log_level")
             .long("log-level")
             .possible_values(&["info", "debug", "error"])
             .default_value("error")
             .help("Logging level."))
        .arg(Arg::with_name("cores")
             .long("cores")
             .default_value("1")
             .value_name("INTEGER")
             .help("Number of cores to use. If larger than 1, will write BAM subsets to temporary files before merging."))
        .arg(Arg::with_name("bam_tag")
             .long("bam-tag")
             .multiple(true)
             .takes_value(true)
             .help("Subset alignments based on one or more tags. If omitted, defaults to CB."));
    args
}

pub struct Locus {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

pub struct Metrics {
    pub total_reads: usize,
    pub barcoded_reads: usize,
    pub kept_reads: usize,
}

pub struct ChunkArgs<'a> {
    barcode_tuples: &'a HashSet<Vec<Vec<u8>>>,
    bam_tags: &'a [String],
    bam_index: usize,
    i: usize,
    bam_file: &'a str,
    tmp_dir: &'a Path,
    virtual_start: Option<i64>,
    virtual_stop: Option<i64>,
}

pub struct ChunkOuts {
    metrics: Metrics,
    out_bam_file: PathBuf,
}

pub struct ProcessedBamOuts {
    metrics: Metrics,
    out_bam_file: PathBuf,
}

fn main() {
    //setup_panic!();  // pretty panics for users
    let mut cli_args = Vec::new();
    for arg in std::env::args_os() {
        cli_args.push(arg.into_string().unwrap());
    }
    _main(cli_args);
}

fn _main(cli_args: Vec<String>) {
    let args = get_args().get_matches_from(cli_args);
    let bam_files: Vec<String> = args
        .values_of("bam")
        .expect("You must provide at least one BAM/CRAM file")
        .map(|s| s.to_string())
        .collect();
    let cell_barcodes = args
        .value_of("cell_barcodes")
        .expect("You must provide a cell barcodes file");
    let out_bam_file = args
        .value_of("out_bam")
        .expect("You must provide a path to write the new BAM file");
    let ll = args.value_of("log_level").unwrap();
    let cores = args
        .value_of("cores")
        .unwrap_or_default()
        .parse::<u64>()
        .expect("Failed to convert cores to integer");
    if cores == 0 {
        error!("--cores must be >= 1");
        process::exit(1);
    }
    let bam_tags: Vec<String> = match args.values_of("bam_tag") {
        Some(vals) => vals.map(|s| s.to_string()).collect(),
        None => vec!["CB".to_string()],
    };

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => {
            println!("Log level not valid");
            process::exit(1);
        }
    };
    let _ = SimpleLogger::init(ll, Config::default());

    check_inputs_exist(&bam_files, cell_barcodes, out_bam_file);
    check_bam_headers_compatible(&bam_files);
    let barcode_tuples = load_barcode_tuples(&cell_barcodes, bam_tags.len()).unwrap();
    let tmp_dir = tempdir().unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cores as usize)
        .build()
        .unwrap();

    // combine metrics
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };

    let mut filtered_bams = Vec::new();
    for (bam_index, bam_file) in bam_files.iter().enumerate() {
        let r = process_single_bam(
            bam_file,
            bam_index,
            &barcode_tuples,
            &bam_tags,
            &cores,
            tmp_dir.path(),
            &pool,
        );
        add_metrics(&mut metrics, &r.metrics);
        filtered_bams.push(r.out_bam_file);
    }

    if metrics.kept_reads == 0 {
        error!("Zero alignments were kept. Does your BAM contain the cell barcodes and/or tag you chose?");
        process::exit(2);
    }

    if filtered_bams.len() == 1 {
        fs::copy(&filtered_bams[0], out_bam_file).unwrap();
    } else {
        info!(
            "Merging {} filtered BAM files into final output",
            filtered_bams.len()
        );
        let final_merge_inputs: Vec<&PathBuf> = filtered_bams.iter().collect();
        merge_bams(final_merge_inputs, Path::new(out_bam_file));
    }

    info!("Done!");
    info!(
        "Visited {} alignments, found {} with barcodes and kept {}",
        metrics.total_reads, metrics.barcoded_reads, metrics.kept_reads
    );
}

fn add_metrics(metrics: &mut Metrics, m: &Metrics) {
    metrics.total_reads += m.total_reads;
    metrics.barcoded_reads += m.barcoded_reads;
    metrics.kept_reads += m.kept_reads;
}

pub fn process_single_bam<'a>(
    bam_file: &'a str,
    bam_index: usize,
    barcode_tuples: &'a HashSet<Vec<Vec<u8>>>,
    bam_tags: &'a [String],
    cores: &u64,
    tmp_dir: &'a Path,
    pool: &rayon::ThreadPool,
) -> ProcessedBamOuts {
    let virtual_offsets = bgzf_noffsets(bam_file, cores).unwrap();
    let mut chunks = Vec::new();
    for (i, (virtual_start, virtual_stop)) in virtual_offsets.iter().enumerate() {
        let c = ChunkArgs {
            barcode_tuples: barcode_tuples,
            bam_tags: bam_tags,
            bam_index: bam_index,
            i: i,
            bam_file: bam_file,
            tmp_dir: tmp_dir,
            virtual_start: *virtual_start,
            virtual_stop: *virtual_stop,
        };
        chunks.push(c);
    }
    let results: Vec<_> = pool.install(|| {
        chunks
            .par_iter()
            .map(|chunk| slice_bam_chunk(chunk))
            .collect()
    });

    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };
    let mut tmp_bams = Vec::new();
    for c in results.iter() {
        add_metrics(&mut metrics, &c.metrics);
        tmp_bams.push(&c.out_bam_file);
    }

    let out_bam_file = tmp_dir.join(format!("filtered_{}.bam", bam_index));
    if *cores == 1 {
        fs::copy(tmp_bams[0], &out_bam_file).unwrap();
    } else {
        info!("Merging {} BAM chunks for {}", tmp_bams.len(), bam_file);
        merge_bams(tmp_bams, &out_bam_file);
    }

    ProcessedBamOuts {
        metrics: metrics,
        out_bam_file: out_bam_file,
    }
}

pub fn check_inputs_exist(bam_files: &[String], cell_barcodes: &str, out_bam_path: &str) {
    for path in bam_files {
        if !Path::new(path).exists() {
            error!("File {} does not exist", path);
            process::exit(1);
        }
    }
    if !Path::new(cell_barcodes).exists() {
        error!("File {} does not exist", cell_barcodes);
        process::exit(1);
    }
    for bam_file in bam_files {
        check_index_exists(bam_file);
    }
    let path = Path::new(out_bam_path);
    if path.exists() {
        error!("Output path already exists");
        process::exit(1);
    }
    if path.is_dir() {
        error!("Output path is a directory");
        process::exit(1);
    }
    let _parent_dir = path.parent();
    if _parent_dir.is_none() {
        error!("Unable to parse directory from {}", out_bam_path);
        process::exit(1);
    }
    let parent_dir = _parent_dir.unwrap();
    if (parent_dir.to_str().unwrap().len() > 0) & !parent_dir.exists() {
        error!("Output directory {:?} does not exist", parent_dir);
        process::exit(1);
    }
}

pub fn check_index_exists(bam_file: &str) {
    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = bam_file.to_owned() + ".bai";
            if !Path::new(&bai).exists() {
                error!("BAM index {} does not exist", bai);
                process::exit(1);
            }
        }
        "cram" => {
            let crai = bam_file.to_owned() + ".crai";
            if !Path::new(&crai).exists() {
                error!("CRAM index {} does not exist", crai);
                process::exit(1);
            }
        }
        &_ => {
            error!(
                "BAM file {} did not end in .bam or .cram. Unable to validate",
                bam_file
            );
            process::exit(1);
        }
    }
}

pub fn check_bam_headers_compatible(bam_files: &[String]) {
    use rust_htslib::bam::Read;
    if bam_files.len() <= 1 {
        return;
    }
    let bam = bam::Reader::from_path(&bam_files[0]).unwrap();
    let template_header = bam.header().as_bytes().to_vec();
    for bam_file in bam_files.iter().skip(1) {
        let bam = bam::Reader::from_path(bam_file).unwrap();
        if bam.header().as_bytes() != template_header.as_slice() {
            error!(
                "BAM headers are incompatible between {} and {}. Unable to merge into one output BAM.",
                bam_files[0], bam_file
            );
            process::exit(1);
        }
    }
}

pub fn load_barcode_tuples(
    filename: impl AsRef<Path>,
    expected_columns: usize,
) -> Result<HashSet<Vec<Vec<u8>>>, Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(filename.as_ref())?;
    let mut bc_set = HashSet::new();
    for (row_index, row) in rdr.records().enumerate() {
        let record = row?;
        if record.len() != expected_columns {
            error!(
                "Barcode row {} has {} columns but {} bam-tag values were provided.",
                row_index + 1,
                record.len(),
                expected_columns
            );
            process::exit(1);
        }
        let mut tuple = Vec::with_capacity(expected_columns);
        for col in 0..expected_columns {
            let field = record.get(col).unwrap().trim();
            if field.is_empty() {
                error!(
                    "Barcode row {} column {} is empty; each bam-tag column requires a barcode value.",
                    row_index + 1,
                    col + 1
                );
                process::exit(1);
            }
            tuple.push(field.as_bytes().to_vec());
        }
        bc_set.insert(tuple);
    }
    let num_bcs = bc_set.len();
    if num_bcs == 0 {
        error!("Loaded 0 barcode tuples. Is your barcode file gzipped or empty?");
        process::exit(1);
    }
    debug!("Loaded {} barcode tuples", num_bcs);
    Ok(bc_set)
}

pub fn get_barcode_tuple(rec: &Record, bam_tags: &[String]) -> Option<Vec<Vec<u8>>> {
    let mut tuple = Vec::with_capacity(bam_tags.len());
    for bam_tag in bam_tags.iter() {
        match rec.aux(bam_tag.as_bytes()) {
            Ok(Aux::String(value)) => tuple.push(value.as_bytes().to_vec()),
            _ => return None,
        }
    }
    Some(tuple)
}

pub fn load_writer(bam: &bam::Reader, out_bam_path: &Path) -> Result<bam::Writer, Error> {
    use rust_htslib::bam::Read; // collides with fs::Read
    let hdr = rust_htslib::bam::Header::from_template(bam.header());
    let out_handle = bam::Writer::from_path(out_bam_path, &hdr, bam::Format::Bam)?;
    Ok(out_handle)
}

pub fn bgzf_noffsets(
    bam_path: &str,
    num_chunks: &u64,
) -> Result<Vec<(Option<i64>, Option<i64>)>, Error> {
    fn vec_diff(input: &Vec<u64>) -> Vec<u64> {
        let vals = input.iter();
        let next_vals = input.iter().skip(1);

        vals.zip(next_vals).map(|(cur, next)| next - cur).collect()
    }

    // if we only have one thread, this is easy
    if *num_chunks == 1 as u64 {
        let final_offsets = vec![(None, None)];
        return Ok(final_offsets);
    }

    let bam_bytes = fs::metadata(bam_path)?.len();
    let mut initial_offsets = Vec::new();
    let step_size = bam_bytes / num_chunks;
    for n in 1..*num_chunks {
        initial_offsets.push((step_size * n) as u64);
    }

    let num_bytes = if initial_offsets.len() > 1 {
        let diff = vec_diff(&initial_offsets);
        let m = diff.iter().max().unwrap();
        cmp::min(1 << 16, *m)
    } else {
        1 << 16
    };

    // linear search to the right of each possible offset until
    // a valid virtual offset is found
    let mut adjusted_offsets = Vec::new();
    let mut fp = fs::File::open(bam_path)?;
    for offset in initial_offsets {
        fp.seek(io::SeekFrom::Start(offset))?;
        let mut buffer = [0; 2 << 16];
        fp.read(&mut buffer)?;
        for i in 0..num_bytes {
            if is_valid_bgzf_block(&buffer[i as usize..]) {
                adjusted_offsets.push(offset + i);
                break;
            }
        }
    }
    // bit-shift and produce start/stop intervals
    let mut final_offsets = Vec::new();

    // handle special case where we only found one offset
    if adjusted_offsets.len() == 1 {
        final_offsets.push((None, None));
        return Ok(final_offsets);
    }

    final_offsets.push((None, Some(((adjusted_offsets[1]) as i64) << 16)));
    for n in 2..num_chunks - 1 {
        let n = n as usize;
        final_offsets.push((
            Some((adjusted_offsets[n - 1] as i64) << 16),
            Some((adjusted_offsets[n] as i64) << 16),
        ));
    }
    final_offsets.push((
        Some(((adjusted_offsets[adjusted_offsets.len() - 1]) as i64) << 16),
        None,
    ));
    Ok(final_offsets)
}

pub fn is_valid_bgzf_block(block: &[u8]) -> bool {
    // look for the bgzip magic characters \x1f\x8b\x08\x04
    // TODO: is this sufficient?
    if block.len() < 18 {
        return false;
    }
    if (block[0] != 31) | (block[1] != 139) | (block[2] != 8) | (block[3] != 4) {
        return false;
    }
    true
}

pub fn slice_bam_chunk(args: &ChunkArgs) -> ChunkOuts {
    let mut bam = bam::Reader::from_path(args.bam_file).unwrap();
    let out_bam_file = args
        .tmp_dir
        .join(format!("{}_{}.bam", args.bam_index, args.i));
    let mut out_bam = load_writer(&bam, &out_bam_file).unwrap();
    let mut metrics = Metrics {
        total_reads: 0,
        barcoded_reads: 0,
        kept_reads: 0,
    };
    for r in bam.iter_chunk(args.virtual_start, args.virtual_stop) {
        let rec = r.unwrap();
        metrics.total_reads += 1;
        let barcode = get_barcode_tuple(&rec, args.bam_tags);
        if let Some(tuple) = barcode {
            metrics.barcoded_reads += 1;
            if args.barcode_tuples.contains(&tuple) {
                metrics.kept_reads += 1;
                out_bam.write(&rec).unwrap();
            }
        }
    }
    let r = ChunkOuts {
        metrics: metrics,
        out_bam_file: out_bam_file,
    };
    info!("Chunk {} is done", args.i);
    r
}

pub fn merge_bams(tmp_bams: Vec<&PathBuf>, out_bam_file: &Path) {
    use rust_htslib::bam::Read; // collides with fs::Read
    let bam = bam::Reader::from_path(tmp_bams[0]).unwrap();
    let mut out_bam = load_writer(&bam, out_bam_file).unwrap();
    for b in tmp_bams.iter() {
        let mut rdr = bam::Reader::from_path(b).unwrap();
        for _rec in rdr.records() {
            let rec = _rec.unwrap();
            out_bam.write(&rec).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as _;
    use std::io::{BufRead, BufReader, Write};
    use tempfile::tempdir;

    fn count_records(path: &str) -> usize {
        let mut bam = bam::Reader::from_path(path).unwrap();
        bam.records().count()
    }

    #[test]
    fn test_bam_single_core() {
        let tmp_dir = tempdir().unwrap();
        let out_file_1 = tmp_dir.path().join("result_1.bam");
        let out_file_2 = tmp_dir.path().join("result_2.bam");
        let out_file_1 = out_file_1.to_str().unwrap();
        let out_file_2 = out_file_2.to_str().unwrap();

        for out_file in [out_file_1, out_file_2].iter() {
            let mut cmds = Vec::new();
            for l in &[
                "subset-bam",
                "-b",
                "test/test.bam",
                "-c",
                "test/barcodes.csv",
                "-o",
                out_file,
                "--cores",
                "1",
            ] {
                cmds.push(l.to_string());
            }
            _main(cmds);
        }

        let bytes_1 = fs::read(out_file_1).unwrap();
        let bytes_2 = fs::read(out_file_2).unwrap();
        assert_eq!(bytes_1, bytes_2);
        assert!(count_records(out_file_1) > 0);
    }

    #[test]
    fn test_multi_bam_single_output() {
        let tmp_dir = tempdir().unwrap();
        let baseline_out = tmp_dir.path().join("baseline.bam");
        let merged_out = tmp_dir.path().join("merged.bam");
        let baseline_out = baseline_out.to_str().unwrap();
        let merged_out = merged_out.to_str().unwrap();

        let mut baseline_cmds = Vec::new();
        for l in &[
            "subset-bam",
            "-b",
            "test/test.bam",
            "-c",
            "test/barcodes.csv",
            "-o",
            baseline_out,
            "--cores",
            "1",
        ] {
            baseline_cmds.push(l.to_string());
        }
        _main(baseline_cmds);

        let mut merged_cmds = Vec::new();
        for l in &[
            "subset-bam",
            "-b",
            "test/test.bam",
            "-b",
            "test/test.bam",
            "-c",
            "test/barcodes.csv",
            "-o",
            merged_out,
            "--cores",
            "1",
        ] {
            merged_cmds.push(l.to_string());
        }
        _main(merged_cmds);

        let baseline_records = count_records(baseline_out);
        let merged_records = count_records(merged_out);
        assert_eq!(merged_records, baseline_records * 2);
    }

    #[test]
    fn test_multi_tag_tuple_match() {
        let tmp_dir = tempdir().unwrap();
        let baseline_out = tmp_dir.path().join("single_tag.bam");
        let multi_tag_out = tmp_dir.path().join("multi_tag.bam");
        let multi_tag_barcodes = tmp_dir.path().join("multi_tag_barcodes.csv");
        let baseline_out = baseline_out.to_str().unwrap();
        let multi_tag_out = multi_tag_out.to_str().unwrap();

        let barcode_fh = fs::File::open("test/barcodes.csv").unwrap();
        let mut writer = fs::File::create(&multi_tag_barcodes).unwrap();
        let reader = BufReader::new(barcode_fh);
        for line in reader.lines() {
            let bc = line.unwrap();
            writeln!(writer, "{0},{0}", bc).unwrap();
        }

        let mut baseline_cmds = Vec::new();
        for l in &[
            "subset-bam",
            "-b",
            "test/test.bam",
            "-c",
            "test/barcodes.csv",
            "-o",
            baseline_out,
            "--cores",
            "1",
        ] {
            baseline_cmds.push(l.to_string());
        }
        _main(baseline_cmds);

        let mut multi_tag_cmds = Vec::new();
        for l in &[
            "subset-bam",
            "-b",
            "test/test.bam",
            "--bam-tag",
            "CB",
            "--bam-tag",
            "CB",
            "-c",
            multi_tag_barcodes.to_str().unwrap(),
            "-o",
            multi_tag_out,
            "--cores",
            "1",
        ] {
            multi_tag_cmds.push(l.to_string());
        }
        _main(multi_tag_cmds);

        let baseline_records = count_records(baseline_out);
        let multi_tag_records = count_records(multi_tag_out);
        assert_eq!(baseline_records, multi_tag_records);
    }
}
