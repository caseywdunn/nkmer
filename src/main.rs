use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use flate2::read::GzDecoder;
use clap::Parser;
use std::path::Path;

use std::collections::HashMap;

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
    return revcomp;
}

fn kmer_revcomp (kmer : u64, k:u32) -> u64 {
    let mut revcomp = 0;
    for i in 0..k {
        let base = (kmer >> (2*i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    return revcomp;
}

fn string2kmers (seq : &String, k:u32) -> Vec<u64> {

    // kmers are encoded as u64, with two bits per base
    // A = 00, C = 01, G = 10, T = 11

    // Create a mask that has 1s in the last 2*k bits
    let mask : u64 = (1 << (2*k)) - 1;

    let mut kmers : Vec<u64> = Vec::new();
    let mut frame : u64 = 0;  // read the bits for each based into the least significatnt end of this integer
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
            n_valid += 1;
            if n_valid >= k {
                kmers.push(frame & mask);
                //kmers.push(kmer_revcomp(frame & mask));
            }
        } else if base==4 {
            frame = 0;
            n_valid = 0;
        } else {
            panic!("Invalid base: {}", c);
        }
    }
    return kmers;
}

fn count_kmers (kmers : &Vec<u64>) -> HashMap<u64, u32> {
    let mut counts : HashMap<u64, u32> = HashMap::new();
    for kmer in kmers {
        let count = counts.entry(*kmer).or_insert(0);
        *count += 1;
    }
    return counts;
}

fn count_histogram (counts : &HashMap<u64, u32>, histo_max : u32) -> Vec<u32> {
    let mut histo : Vec<u32> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for (_kmer, count) in counts {
        if *count <= histo_max {
            histo[*count as usize] += 1;
        }
        else {
            histo[histo_max as usize + 1] += 1;
        }
    }
    return histo;
}

fn histogram_string (histo : &Vec<u32>) -> String {
    let mut s = String::new();
    for (i, count) in histo.iter().enumerate() {
        if i > 0 {
            s.push_str(&format!("{} {}\n", i, count));
        }
    }
    return s;
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
    histo_max: u32,

    /// Number of chunks to divide the data into
    #[arg(long, default_value_t = 10)]
    n: usize,

    /// Name used for analysis output
    #[arg(short, long, default_value_t = String::from("sample") )]
    output: String,

    /// Input files, gzipped fastq
    #[arg(required = true)]
    input: Vec<String>,
}


fn main() {

    // Ingest all command line arguments as input files
    //let args: Vec<String> = std::env::args().collect();
    let args = Args::parse();

    assert!(args.k < 32, "k must be less than 32 due to use of 64 bit integers to encode kmers");

    println!("Reading input...");
    // Iterate over all input files and parse all the records into a Vec
    let mut seqs: Vec<String> = Vec::new();
    let mut n_records_all = 0;
    let mut n_bases_all = 0;

    for file_name in args.input {
    //for file_name in &args[1..] {
        let mut n_records = 0;
        let mut n_bases = 0;
        
        let start = std::time::Instant::now();
        let path = Path::new(&file_name);
        let file = File::open(path).expect("Ooops.");
        let reader = BufReader::new(GzDecoder::new(file));

        let mut line_n = 0;

        for line in reader.lines() {
            if line_n % 4 == 1 {
                let line = line.unwrap();
                n_records += 1;
                n_bases += line.len();
                seqs.push(line);
            }
            line_n += 1;
        }

        let stop = std::time::Instant::now();

    
        println!("File {}: records {}, bases {}", file_name, n_records, n_bases);
        println!("Time: {} ms", (stop - start).as_millis());

        n_records_all += n_records;
        n_bases_all += n_bases;
    }

    println!("Total: records {}, bases {}", n_records_all, n_bases_all);

    let start = std::time::Instant::now();
    println!("Extracting kmers...");

    let mut kmers : Vec<u64> = Vec::new();
    for seq in seqs {
        // Call string2kmers on each sequence and append the result to kmers
        let seq_rc = revcomp(&seq);
        kmers.append(&mut string2kmers(&seq, args.k));
        kmers.append(&mut string2kmers(&seq_rc, args.k));
    }

    println!("Total kmers: {}", kmers.len());

    let stop = std::time::Instant::now();

    println!("Time: {} ms", (stop - start).as_millis());

    let start = std::time::Instant::now();
    println!("Counting all kmers...");

    //let counts = count_kmers(&kmers);
    //println!("Unique kmers: {}", counts.len());

    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());


    let start = std::time::Instant::now();
    println!("Counting kmers in chunks...");

    let chunk_size = kmers.len() / args.n;
    let mut cum_counts : HashMap<u64, u32> = HashMap::new();

    for i in 0..args.n {
        let start: usize = i * chunk_size;
        let end: usize = (i+1) * chunk_size;
        let counts = count_kmers(&kmers[start..end].to_vec());
        for (kmer, count) in counts {
            let cum_count = cum_counts.entry(kmer).or_insert(0);
            *cum_count += count;
        }

        let histo = count_histogram(&cum_counts, args.histo_max);
        println!("Part {}: unique kmers {}", i, cum_counts.len());
        let out = histogram_string(&histo);
        // println!("{out}");

        // Write the histogram to a file
        let file_name =  &format!("{output}_k{k}_part{i}.histo", output=args.output, k=args.k);
        let out_path = Path::new(&file_name);
        let mut file = File::create(file_name).unwrap();
        file.write_all(out.as_bytes());
    }

    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());


}

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
}