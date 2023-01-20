use std::array;
use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use flate2::read::GzDecoder;

use std::collections::HashMap;

const PARTS : u32 = 100;
const K : u32 = 21;

// Function that takes a DNA sequence as a string and returns reverse complement
fn revcomp (seq : String) -> String {
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

fn kmer_revcomp (kmer : u64) -> u64 {
    let mut revcomp = 0;
    for i in 0..K {
        let base = (kmer >> (2*i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    return revcomp;
}

fn string2kmers (seq : String) -> Vec<u64> {

    // kmers are encoded as u64, with two bits per base
    // A = 00, C = 01, G = 10, T = 11

    // Create a mask that has 1s in the last 2*K bits
    let mask : u64 = (1 << (2*K)) - 1;

    let mut kmers : Vec<u64> = Vec::new();
    let mut frame : u64 = 0;
    let mut n_valid = 0; // number of valid bases in the frame
    for (i, c) in seq.chars().enumerate() {
        let base = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            'N' => 4,
            _ => 5,
        };
        if (base < 4) {
            frame = (frame << 2) | base;
            n_valid += 1;
            if (n_valid >= K) {
                kmers.push(frame & mask);
                //kmers.push(kmer_revcomp(frame & mask));
            }
        } else if (base==4) {
            frame = 0;
            n_valid = 0;
        } else {
            panic!("Invalid base: {}", c);
        }
    }
    return kmers;
}
fn main() {

    // Ingest all command line arguments as input files
    let args: Vec<String> = std::env::args().collect();

    // Iterate over all input files and parse all the records into a Vec
    let mut seqs: Vec<String> = Vec::new();
    let mut n_records_all = 0;
    let mut n_bases_all = 0;

    for file_name in &args[1..] {
        let mut n_records = 0;
        let mut n_bases = 0;
        
        let start = std::time::Instant::now();
        let file = File::open(file_name).expect("Ooops.");
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
    println!("Extracting kmers");

    let mut kmers : Vec<u64> = Vec::new();
    for seq in seqs {
        // Call string2kmers on each sequence and append the result to kmers
        let seq_rc = revcomp(seq.clone());
        kmers.append(&mut string2kmers(seq));
        kmers.append(&mut string2kmers(seq_rc));
    }

    println!("Total kmers: {}", kmers.len());

    let stop = std::time::Instant::now();

    println!("Time: {} ms", (stop - start).as_millis());

    let start = std::time::Instant::now();
    println!("Counting kmers");

    let mut counts : HashMap<u64, u32> = HashMap::new();
    for kmer in kmers {
        let count = counts.entry(kmer).or_insert(0);
        *count += 1;
    }

    let stop = std::time::Instant::now();

    println!("Time: {} ms", (stop - start).as_millis());

    // print the number of unique kmers
    println!("Unique kmers: {}", counts.len());


}
