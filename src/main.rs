use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use flate2::read::GzDecoder;

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
}
