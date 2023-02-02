# nkmer

A tool for building kmer histograms from cumulative subsets of data from fastq.gz files

## Installation

Download and install [Rust and Cargo](https://www.rust-lang.org/tools/install).

Then, in this repo build the binary:

    cargo build --release

The executable is then at `target/debug/nkmer`

## Usage

To get full usage information, run

    nkmer --help

An example analysis would look like this:

    nkmer -k 21 -n 10 --histo-max 10000 -o Agalma-elegans agalma_1_R1.fastq.gz agalma_1_R2.fastq.gz agalma_2_R1.fastq.gz agalma_2_R2.fastq.gz

This ingests all the fastq.gz files and splits the reads into n chunks. It then generates 
the n cumulative histograms of kmer frequency across the chunks.

The output histogram files in this case would be:

    Agalma-elegans_k21_part0.histo
	Agalma-elegans_k21_part1.histo
	...
	Agalma-elegans_k21_part9.histo

Each histo includes the number of kmers with frequencies up to `--histo-max`.

## Implementation

`nkmer` uses an unsigned 64 bit integer representation of kmers, where each nucleotide is encoded as two bits. This limits the length of kmers to 31,
since k must be odd (to avoid pailindromes) and the two bits per nucleotide must fit within the integer length.

## Background

For an excellent summary of kmer tool optimization, see [Manekar and Sathe 2018](https://academic.oup.com/gigascience/article/7/12/giy125/5140149).

## Development

Some common tasks in development:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release

`zlib` dependency requires that cmake is installed and in the `PATH`. For example, on MacOS you need:

    export PATH="/Applications/CMake.app/Contents/bin:$PATH"

### Test data

For explroeing, testing, and debugging this tool, download the this set of  [E. coli genome reads](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR23291413&display=download) as a fastq and place it in a `data/` folder in the root of this repo.