# nkmer

A tool for building kmer histograms from cumulative subsets of data from fastq.gz files

## Installation

Download and install [Rust and Cargo](https://www.rust-lang.org/tools/install).

Then, in this repo build the binary:

    cargo build

The executable is then at `target/debug/nkmer`

## Usage

To get full usage information, run

    nkmer --help

An example analysis would look like this:

    nkmer -k 21 -n 10 --histo-max 10000 -o Agalma-elegans agalma_1_R1.fastq.gz agalma_1_R2.fastq.gz agalma_2_R1.fastq.gz agalma_2_R2.fastq.gz

This ingests all the fastq.gz files and extracts their kmers. It then splits the kmers into n chunks, and generates 
the n cumulative histograms of kmer frequency as each chunk is considered with all those already considered.

The output histogram files in this case would be:

    Agalma-elegans_k21_part0.histo
	Agalma-elegans_k21_part1.histo
	...
	Agalma-elegans_k21_part9.histo

Each histo includes the number of kmers with frequencies up to `--histo-max`.

## Implementation

`nkmer` uses an unsigned 64 bit integer representation of kmers, where each nucleotide is encoded as two bits. This limits the length of kmers to 31,
since k must be odd (to avoid pailindromes) and the two bits per nucleotide must fit within the integer length.