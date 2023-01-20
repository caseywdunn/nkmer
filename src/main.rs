extern crate needletail;
use std::borrow::Cow;

use needletail::{parse_fastx_file, Sequence, FastxReader};

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn main() {

    // Ingest all command line arguments as input files
    let args: Vec<String> = std::env::args().collect();

    // Iterate over all input files and parse all the records into a Vec
    let mut seqs: Vec<FastxReader::FastxRecord<'static>> = Vec::new();
    let mut n_records_all = 0;
    let mut n_bases_all = 0;

    for filename in &args[1..] {
        //let records = parse_fastx_file(&filename).unwrap();
        let reader = parse_fastx_file(&filename).expect("valid path/file");
        let mut n_records = 0;
        let mut n_bases = 0;
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");

            // Increment the number of records
            n_records += 1;
            n_bases += seqrec.num_bases();
            seqs.push(seqrec);
            //let seq = seqrec.seq();
            //seqs.push(seqrec.seq().to_string());
            //seqs.push(seq);
            //print_type_of(&seqrec);
        }
        println!("File {}: records {}, bases {}", filename, n_records, n_bases);

        n_records_all += n_records;
        n_bases_all += n_bases;
    }


    println!("Total: records {}, bases {}", n_records_all, n_bases_all);
}
