use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use flate2::read::GzDecoder;
use clap::Parser;
use std::path::Path;
use std::collections::HashSet;
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

        let mut line_n = 0;

        /* 
        let reader = BufReader::new(GzDecoder::new(file));
        for line in reader.lines() {
            if line_n % 4 == 1 {
                let line = line.unwrap();
                n_records += 1;
                n_bases += line.len();
                seqs.push(line);
            }
            line_n += 1;
        } */

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
                    n_bases += line.len();
                    seqs.push(line);
                }
                line_n += 1;
            }
        
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
        kmers.append(&mut string2kmers(&seq, args.k));
    }

    println!("Total kmers: {}", kmers.len());
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());

    /* let start = std::time::Instant::now();
    println!("Counting all kmers...");
    let counts_all = count_kmers(&kmers[..].to_vec());
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis()); */

    let start = std::time::Instant::now();
    println!("Counting kmers in chunks...");

    let chunk_size = kmers.len() / args.n;
    let mut cumulative_counts : HashMap<u64, u32> = HashMap::new();

    for i in 0..args.n {
        let start: usize = i * chunk_size;
        let end: usize = (i+1) * chunk_size;
        let counts = count_kmers(&kmers[start..end].to_vec());
        for (kmer, count) in counts {
            let cum_count = cumulative_counts.entry(kmer).or_insert(0);
            *cum_count += count;
        }

        let histo = count_histogram( &cumulative_counts, args.histo_max);
        println!("Part {}: unique kmers {}", i, cumulative_counts.len());
        let out = histogram_string(&histo);

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