use std::array;
use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use flate2::read::GzDecoder;

use std::collections::HashMap;

const PARTS : u32 = 100;
const K : u32 = 21;

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

    let mut kmers : Vec<u64> = Vec::new();
    for seq in seqs {
        // Call string2kmers on each sequence and append the result to kmers
        kmers.append(&mut string2kmers(seq));
    }

    println!("Total kmers: {}", kmers.len());

    let stop = std::time::Instant::now();

    println!("Time: {} ms", (stop - start).as_millis());


}
