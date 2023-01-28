use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use clap::Parser;
use std::path::Path;
use std::collections::HashMap;

#[inline(always)]
fn string2kmers (seq : &str, k:u32) -> Vec<u64> {

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

        match base {
            0 | 1 | 2 | 3 => {
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
            },
            4 => {
                frame = 0;
                revframe = 0;
                n_valid = 0;
            },
            _ => panic!("Invalid base: {}", c),
        }
    }
    kmers
}

fn count_histogram (counts : &HashMap<u64, u64>, histo_max : u64) -> Vec<u32> {
    let mut histo : Vec<u32> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for count in counts.values() {
        if *count <= histo_max {
            histo[*count as usize] += 1;
        }
        else {
            histo[histo_max as usize + 1] += 1;
        }
    }
    histo
}

fn histogram_string (histo : &[u32]) -> String {
    let mut s = String::new();
    for (i, count) in histo.iter().enumerate() {
        if i > 0 {
            s.push_str(&format!("{} {}\n", i, count));
        }
    }
    s
}

fn get_fastq_stats(input_files : &Vec<String>, max_reads : u64) -> (u64, u64){
    let mut n_records_all:u64 = 0;
    let mut n_bases_all:u64 = 0;

    'processing_files: for file_name in input_files {
        let mut n_records:u64 = 0;
        let mut n_bases:u64 = 0;
        
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
                    n_records_all += 1;
                    if max_reads > 0 && n_records_all > max_reads {
                        println!("Reached maximum number of reads: {}", max_reads);
                        println!("File {}: records {}, bases {}, mean read length {}", file_name, n_records, n_bases, n_bases as f64 / n_records as f64);
                        break 'processing_files;
                    }
                    let line_len = line.len() as u64;
                    n_bases += line_len;
                    n_bases_all += line_len;
                }
                line_n += 1;
            }
        
        }

        println!("File {}: records {}, bases {}, mean read length {}", file_name, n_records, n_bases, n_bases as f64 / n_records as f64);
    }

    (n_records_all, n_bases_all)
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
    #[arg(short, default_value_t = 10)]
    n: usize,

    /// Maximum number of reads to process
    #[arg(short, long, default_value_t = 0)]
    max_reads: u64,

    /// Directory and filename prefix for analysis output, for example out_dir/Nanomia-bijuga
    #[arg(short, long, default_value_t = String::from("sample") )]
    output: String,

    /// Input files, gzipped fastq
    #[arg(required = true)]
    input: Vec<String>,
}

fn main() {

    // Ingest command line arguments
    let args = Args::parse();

    // Print the program name and version
    println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    // Print the arguments
    println!("{:?}", args);

    // Parse the output path and create directories if necessary
    let path = Path::new(&args.output);
    let out_name = path.file_name().unwrap().to_str().unwrap(); // This is the prefix of the output files
    let directory = path.parent().unwrap().to_str().unwrap();
    let _ = std::fs::create_dir_all(directory);

    // Check that the arguments are valid
    assert!(args.k < 32, "k must be less than 32 due to use of 64 bit integers to encode kmers");
    assert!(args.k > 0, "k must be greater than 0");
    assert!(args.k % 2 == 1, "k must be odd");
    assert!(args.histo_max > 0, "histo_max must be greater than 0");
    assert!(args.n > 0, "n must be greater than 0");

    // Get the sequence statistics
    println!("Getting input file statistics...");
    let start = std::time::Instant::now();
    let (n_records_all, n_bases_all) = get_fastq_stats(&args.input, args.max_reads);
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());
    println!("Total records {}, total bases {}, mean read length {}", n_records_all, n_bases_all, n_bases_all as f64 / n_records_all as f64);

    // Ingest the data and count kmers
    let mut n_bases_processed:u64 = 0;
    let mut n_records_processed:u64 = 0;

    let chunk_size = n_records_all / args.n as u64;
    let mut kmer_counts : HashMap<u64, u64> = HashMap::new();
    let start = std::time::Instant::now();
    let mut chunk_i = 0;
    'processing_files: for file_name in args.input {
        
        let path = Path::new(&file_name);
        let file = File::open(path).expect("Ooops.");

        let mut line_n = 0;

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
                        
                        let histo = count_histogram( &kmer_counts, args.histo_max);
                        let histo_text = histogram_string(&histo);
                
                        // Write the histogram to a file
                        let file_histo_name =  &format!("{out_name}_k{k}_part{chunk_i}.histo", k=args.k);
                        let file_histo_path = Path::new(directory).join(file_histo_name);
                        let mut file = File::create(file_histo_path).unwrap();
                        file.write_all(histo_text.as_bytes()).expect("Couldn't write output file");

                        // Write the stats to a file
                        let file_stats_name =  &format!("{out_name}_k{k}_part{chunk_i}.stats.tsv", k=args.k);
                        let file_stats_path = Path::new(directory).join(file_stats_name);
                        let mut file_stats = File::create(file_stats_path).unwrap();
                        writeln!(&mut file_stats, "sample	fastq_index	mean_readlen	num_reads	gigabases").unwrap();
                        writeln!(&mut file_stats, "{}   {}  {}  {}  {}", out_name, chunk_i, n_bases_processed/n_records_processed, n_records_processed, n_bases_processed as f64 / 1_000_000_000.).unwrap();

                        chunk_i += 1;
                    }

                    if n_records_processed >= (args.n as u64 * chunk_size) {
                        println!("  Skipping last {} reads after the last full chunk...", n_records_all - n_records_processed);
                        break 'processing_files;
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
    use std::collections::HashSet;

    // Function that takes a DNA sequence as a string and returns reverse complement
    fn revcomp (seq : &str) -> String {
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