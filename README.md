# subset-bam

`subset-bam` is a tool to subset a 10x Genomics BAM file based on a tag, most commonly the cell barcode tag. The latest version is `v1.1.0` and can be found [on the releases page](https://github.com/10XGenomics/subset-bam/releases).

## Overview of how it works
`subset-bam` is a simple tool implemented in Rust that takes a 10x Genomics BAM file, a CSV file defining the subset of cells you want to isolate, and produces a new BAM file with only alignments associated with those cells.

In the subsetting process, temporary BAM files will be written to your temporary file (`$TMPDIR`) location before a final concatenation step. Please make sure this location is writeable has enough space to support this operation. If you cannot write to this location, you can set your `TMPDIR` variable to a path you can write by the command `export TMPDIR=/my/writeable/path`.

## Support
This tool is not officially supported. If you have any comments, please submit a GitHub issue.

## Installation

`subset-bam` has automatically generated downloadable binaries for generic linux and Mac OSX under the [releases page](https://github.com/10XGenomics/subset-bam/releases). The linux binaries are expected to work on [our supported Operating Systems](https://support.10xgenomics.com/os-support). 

## Compiling from source
`subset-bam` is standard Rust executable project, that works with stable Rust >=1.13. Install Rust through the standard channels, then type `cargo build --release`. The executable will appear at `target/release/subset-bam`. As usual it's important to use a release build to get good performance.

## Usage

`--bam (-b)`: Input 10x Genomics BAM/CRAM. This argument can be provided multiple times to process multiple files and merge the filtered output into a single BAM. Each input must have an index (`.bai`/`.crai`). REQUIRED.

`--cell-barcodes (-c)`: A barcode whitelist file. If one `--bam-tag` is provided (or default `CB` is used), this can be one barcode per line (backward compatible). If multiple `--bam-tag` values are provided, this file must have the same number of CSV columns as tags, where each row defines one tag-tuple that must match together. REQUIRED.

`--out-bam (-o)`: A path to write the subsetted BAM file to. REQUIRED.

`--cores`: Number of parallel cores to use. DEFAULT: 1.

`--log-level`: One of `info`, `error` or `debug`. Increasing levels of logging. DEFAULT: error.

`--bam-tag`: Can be provided multiple times. Reads are kept only if all requested tags are present and the full tuple of tag values appears in the barcode whitelist file. If omitted, defaults to `CB`.

## License
`subset-bam` is licensed under the [MIT license](http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.
