use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use flate2::read::GzDecoder;
use clap::Parser;
use std::path::Path;
use std::collections::HashSet;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use dashmap::DashMap;
use fxhash::FxHasher;
use rustc_hash::FxHashMap;
use nohash_hasher::NoHashHasher;

// Function that takes a DNA sequence as a string and returns reverse complement
fn revcomp (seq : &String) -> String {
    let mut revcomp = String::new();
    for c in seq.chars().rev() {
        match c {
            'A' => revcomp.push('T'),
            'C' => revcomp.push('G'),
            'G' => revcomp.push('C'),
            'T' => revcomp.push('A'),
            'N' => revcomp.push('N'),
            _ => panic!("Invalid base: {}", c),
        }
    }
    revcomp
}

fn kmer_revcomp (kmer : u64, k:u32) -> u64 {
    let mut revcomp = 0;
    for i in 0..k {
        let base = (kmer >> (2*i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    revcomp
}

#[inline(always)]
fn string2kmers (seq : &String, k:u32) -> Vec<u64> {

    // kmers are encoded as u64, with two bits per base
    // A = 00, C = 01, G = 10, T = 11

    // Create a mask that has 1s in the last 2*k bits
    let mask : u64 = (1 << (2*k)) - 1;

    let mut kmers : Vec<u64> = Vec::new();
    let mut frame : u64 = 0;  // read the bits for each base into the least significant end of this integer
    let mut revframe : u64 = 0;  // read the bits for complement into the least significant end of this integer
    let mut n_valid = 0; // number of valid bases in the frame
    for (_i, c) in seq.chars().enumerate() {
        let base = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            'N' => 4,
             _  => 5,
        };
        if base < 4 {
            frame = (frame << 2) | base;

            revframe = (revframe >> 2) | ((3 - base) << (2*(k-1)));

            n_valid += 1;
            if n_valid >= k {

                let forward = frame & mask;
                let reverse = revframe & mask;

                if forward < reverse {
                    kmers.push(forward);
                } else {
                    kmers.push(reverse);
                }
            }
        } else if base==4 {
            frame = 0;
            revframe = 0;
            n_valid = 0;
        } else {
            panic!("Invalid base: {}", c);
        }
    }
    kmers
}

fn count_histogram (counts : &HashMap<u64, u64>, histo_max : u64) -> Vec<u32> {
    let mut histo : Vec<u32> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for (_kmer, count) in counts {
        if *count <= histo_max {
            histo[*count as usize] += 1;
        }
        else {
            histo[histo_max as usize + 1] += 1;
        }
    }
    histo
}

fn histogram_string (histo : &Vec<u32>) -> String {
    let mut s = String::new();
    for (i, count) in histo.iter().enumerate() {
        if i > 0 {
            s.push_str(&format!("{} {}\n", i, count));
        }
    }
    s
}

/// Count k-mers in a set of fastq.gz files, with an option to assess cumulative subsets
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// k-mer length
    #[arg(short, default_value_t = 21)]
    k: u32,

    /// Maximum value for histogram
    #[arg(long, default_value_t = 10000)]
    histo_max: u64,

    /// Number of chunks to divide the data into
    #[arg(long, default_value_t = 10)]
    n: usize,

    /// Maximum number of reads to process
    #[arg(short, long, default_value_t = 0)]
    max_reads: u64,

    /// Name used for analysis output
    #[arg(short, long, default_value_t = String::from("sample") )]
    output: String,

    /// Input files, gzipped fastq
    #[arg(required = true)]
    input: Vec<String>,
}

fn get_fastq_stats(input_files : &Vec<String>) -> (u64, u64){
    let mut n_records_all:u64 = 0;
    let mut n_bases_all:u64 = 0;

    for file_name in input_files {
        let mut n_records:u64 = 0;
        let mut n_bases:u64 = 0;
        
        let start = std::time::Instant::now();
        let path = Path::new(&file_name);
        let file = File::open(path).expect("Ooops.");

        let mut line_n:u64 = 0;

        // From https://github.com/rust-lang/flate2-rs/issues/41#issuecomment-219058833
        // Handle gzip files with multiple blocks

        let mut reader = BufReader::new(file);
        loop {
            //loop over all possible gzip members
            match reader.fill_buf() {
                Ok(b) => if b.is_empty() { break },
                Err(e) => panic!("{}", e)
            }
        
            //decode the next member
            let gz = flate2::bufread::GzDecoder::new(&mut reader);
            let block_reader = BufReader::new(gz);
            for line in block_reader.lines() {
                if line_n % 4 == 1 {
                    let line = line.unwrap();
                    n_records += 1;
                    n_bases += line.len() as u64;
                }
                line_n += 1;
            }
        
        }

        n_records_all += n_records;
        n_bases_all += n_bases;

        println!("File {}: records {}, bases {}, mean read length {}", file_name, n_records, n_bases, n_bases as f64 / n_records as f64);
    }

    (n_records_all, n_bases_all)
}


fn main() {

    // Ingest all command line arguments as input files
    //let args: Vec<String> = std::env::args().collect();
    let args = Args::parse();

    assert!(args.k < 32, "k must be less than 32 due to use of 64 bit integers to encode kmers");

    println!("Reading input...");


    let mut n_records_all:u64 = 0;
    let mut n_bases_all:u64 = 0;

    if args.max_reads == 0 {
        println!("Getting input file statistics...");
        let start = std::time::Instant::now();
        let (n_records, n_bases) = get_fastq_stats(&args.input);
        let stop = std::time::Instant::now();
        println!("Time: {} ms", (stop - start).as_millis());
        
        n_records_all = n_records;
        n_bases_all = n_bases;
        println!("Total records {}, total bases {}, mean read length {}", n_records_all, n_bases_all, n_bases_all as f64 / n_records_all as f64)
    } else {
        n_records_all = args.max_reads;
        n_bases_all = 0;
        println!("Input file statistics not calculated, using max_reads = {} as the number of reads. Number of bases unknown.", args.max_reads);
    }

    let mut n_bases_processed:u64 = 0;
    let mut n_records_processed:u64 = 0;

    let chunk_size = n_records_all / args.n as u64;
    // https://docs.rs/nohash-hasher/latest/nohash_hasher/struct.NoHashHasher.html
    let mut kmer_counts : HashMap<u64, u64, BuildHasherDefault<NoHashHasher<u64>>> = HashMap::with_capacity_and_hasher(1000000, BuildHasherDefault::default());
    //let mut kmer_counts: FxHashMap<u64, u64> = FxHashMap::default();
    // let mut kmer_counts_dash : DashMap<u64, u64, BuildHasherDefault<FxHasher>> = DashMap::new();

    // Ingest the data and count kmers
    println!("Counting kmers...");
    let start = std::time::Instant::now();
    let mut chunk_i = 0;
    'processing_files: for file_name in args.input {
        
        let path = Path::new(&file_name);
        let file = File::open(path).expect("Ooops.");

        let mut line_n:u64 = 0;

        // From https://github.com/rust-lang/flate2-rs/issues/41#issuecomment-219058833
        // Handle gzip files with multiple blocks

        let mut reader = BufReader::new(file);
        loop {
            //loop over all possible gzip members
            match reader.fill_buf() {
                Ok(b) => if b.is_empty() { break },
                Err(e) => panic!("{}", e)
            }
        
            //decode the next member
            let gz = flate2::bufread::GzDecoder::new(&mut reader);
            let block_reader = BufReader::new(gz);
            for line in block_reader.lines() {
                if line_n % 4 == 1 {
                    let line = line.unwrap();
                    n_records_processed += 1;
                    n_bases_processed += line.len() as u64;

                    // Count the kmers
                    let kmers = string2kmers(&line, args.k);
                    for kmer in kmers {
                        let count = kmer_counts.entry(kmer).or_insert(0);
                        *count += 1;
                    }

                    // If we've processed enough records, write the output
                    if (n_records_processed % chunk_size == 0) & (n_records_processed > 0) {
                        println!("Cumulative chunk {} stats: {} reads, {} bases, {} unique kmers, {} kmers.", chunk_i, n_records_processed, n_bases_processed, kmer_counts.len(), kmer_counts.values().sum::<u64>());
                        let start = std::time::Instant::now();
                        
                        //let histo = count_histogram( &kmer_counts, args.histo_max);

                        let mut histo : Vec<u32> = vec![0; args.histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
                        for (_kmer, count) in kmer_counts {
                            if count <= args.histo_max {
                                histo[count as usize] += 1;
                            }
                            else {
                                histo[args.histo_max as usize + 1] += 1;
                            }
                        }

                        let out = histogram_string(&histo);
                
                        // Write the histogram to a file
                        let file_name =  &format!("{output}_k{k}_part{chunk_i}.histo", output=args.output, k=args.k);
                        let out_path = Path::new(&file_name);
                        let mut file = File::create(file_name).unwrap();
                        file.write_all(out.as_bytes());

                        chunk_i += 1;

                        if chunk_i >= args.n {
                            break 'processing_files;
                        }
                    }
                }
                line_n += 1;
            }
        }
    }

    println!("Total processed: records {}, bases {}", n_records_processed, n_bases_processed);
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tests(){
        assert_eq!(2+2,4);
    }

    #[test]
    fn test_revcomp(){
        let seq = String::from("ACGNTT");
        let rc = revcomp(&seq);
        assert_eq!(rc, "AANCGT");
    }

    #[test]
    fn test_strings2kmer(){
        let seq = String::from("ACGNTT");
        let kmer = string2kmers(&seq, 3);

        // The only valid 3-mer in this sequence is ACG
        // A = 00, C = 01, G = 10, T = 11
        // ACG = 000110
        // CGT = 011011 , need to check reverse complement
        // It is the smallest that will be retained
        assert_eq!(kmer[0], 0b0000000000000000000000000000000000000000000000000000000000000110);

        // Check reverse complement of the sequence, which should be the same
        let seq_rc = revcomp(&seq);
        let kmer_rc = string2kmers(&seq_rc, 3);
        assert_eq!(kmer_rc[0], 0b0000000000000000000000000000000000000000000000000000000000000110);
    }

    #[test]
    fn test_strings2kmer_revcomp(){
        let seq = String::from("ACGCTNCGTT");
        let kmers_f = string2kmers(&seq, 3);
        let kmers_r = string2kmers(&revcomp(&seq), 3);
        // let kmers_r = string2kmers(&String::from("ACTGCTNCGTT"), 3); // This should fail


        // Check that the kmers are the same
        // Convert the kmers to hash sets
        let kmers_f_set: HashSet<u64> = kmers_f.iter().cloned().collect();
        let kmers_r_set: HashSet<u64> = kmers_r.iter().cloned().collect();

        // Check that the sets are equal
        assert_eq!(kmers_f_set, kmers_r_set);


    }
}