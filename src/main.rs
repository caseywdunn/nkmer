use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use clap::Parser;
use std::path::Path;
use std::collections::HashMap;

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

fn main() {

    // Ingest all command line arguments as input files
    let args = Args::parse();

    assert!(args.k < 32, "k must be less than 32 due to use of 64 bit integers to encode kmers");
    assert!(args.k > 0, "k must be greater than 0");
    assert!(args.k % 2 == 1, "k must be odd");
    assert!(args.histo_max > 0, "histo_max must be greater than 0");
    assert!(args.n > 0, "n must be greater than 0");

    // Get the sequence vector
    println!("Reading input files...");
    let start = std::time::Instant::now();
    
    let mut n_records_all:u64 = 0;
    let mut n_bases_all:u64 = 0;
    let mut seqs: Vec<String>  = Vec::new();

    for file_name in args.input {
        let mut n_records:u64 = 0;
        let mut n_bases:u64 = 0;
        
        let _start = std::time::Instant::now();
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
                    seqs.push(line);
                }
                line_n += 1;
            }
        
        }

        n_records_all += n_records;
        n_bases_all += n_bases;

        println!("File {}: records {}, bases {}, mean read length {}", file_name, n_records, n_bases, n_bases as f64 / n_records as f64);
    }

    println!("Total processed: records {}, bases {}, mean read length {}", n_records_all, n_bases_all, n_bases_all as f64 / n_records_all as f64);
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());

    let chunk_size = seqs.len() / args.n;

    // Count the kmers in chunks
    println!("Counting kmers in each chunk...");
    let start = std::time::Instant::now();

    // Create a vector of hashmaps to store the counts
    let mut kmer_counts_chunks : Vec<HashMap<u64, u64>> = Vec::new();

    for part in 0..args.n{
        let mut kmer_counts_part : HashMap<u64, u64> = HashMap::new();

        let start: usize = part * chunk_size;
        let end: usize = (part + 1) * chunk_size;
        for seq in &seqs[start..end] {
            let kmers_line = string2kmers(seq, args.k);
            for kmer in kmers_line {
                let count = kmer_counts_part.entry(kmer).or_insert(0);
                *count += 1;
            }
        }
        kmer_counts_chunks.push(kmer_counts_part);
    }

    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());

    println!("Integrating chunks...");
    let start = std::time::Instant::now();
    let mut kmer_counts : HashMap<u64, u64> = HashMap::new();

    let mut chunk_i =0;
    
    let histo_max = args.histo_max;

    for kmer_count_chunk in kmer_counts_chunks {
        for (kmer, count) in kmer_count_chunk {
            let count_total = kmer_counts.entry(kmer).or_insert(0);
            *count_total += count;
        }

        let mut histo : Vec<u32> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
        for count in kmer_counts.values() {
            if count <= &histo_max {
                histo[*count as usize] += 1;
            }
            else {
                histo[args.histo_max as usize + 1] += 1;
            }
        }

        let out = histogram_string(&histo);

        // Write the histogram to a file
        let file_name =  &format!("{output}_k{k}_part{chunk_i}.histo", output=args.output, k=args.k);
        let _out_path = Path::new(&file_name);
        let mut file = File::create(file_name).unwrap();
        file.write_all(out.as_bytes()).map_err(|err| println!("{:?}", err)).ok();

        chunk_i += 1;
    }

    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

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