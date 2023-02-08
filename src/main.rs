use clap::Parser;
use csv::WriterBuilder;
use find_peaks::PeakFinder;
use ndarray::{Array1, Array2};
use ndarray_csv::Array2Writer;
use plotters::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

// Takes a histo vector slice and counts the total number of kmers
fn count_kmers_in_histo(histo: &[i32]) -> i32 {
    let mut total_kmers = 0;
    for (i, count) in histo.iter().enumerate() {
        total_kmers += count * i as i32;
    }
    total_kmers
}

fn infer_stats(
    histo_chunks: &Array2<u64>,
    _n_bases_chunks: &Array1<u64>,
    _n_records_chunks: &Array1<u64>,
    k: u32,
    chunk_i: usize,
    directory: &str,
    out_name: &str,
) {
    // Create the histogram and save it to a file file
    // Based on histogram example from https://plotters-rs.github.io/book/basic/basic_data_plotting.html

    let histo_64 = histo_chunks.row(chunk_i).to_vec();
    let histo = histo_64.iter().map(|v| *v as i32).collect::<Vec<i32>>();

    // Find the peaks in the histogram
    let mut peak_finder = PeakFinder::new(&histo);
    peak_finder.with_min_prominence(1000);
    let peaks = peak_finder.find_peaks();

    // Determine the x_max of the plot based on peak location
    // Set the xmax
    let x_max = 50;

    // Run analyses depending on the number of peaks found
    if peaks.is_empty() {
        println!("Peaks: None");
    } else if peaks.len() == 1 {
        let peak_errors = peaks[0].position.start;
        println!("Peaks: Rare peak only at {}", peak_errors);
    } else if peaks.len() == 2 {
        // Can apply manual genome size estimation per https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/

        let peak_errors = peaks[0].position.start;
        let peak_genome = peaks[1].position.start;

        // Need to get the first trough. Can do this by taking the negative of the histogram and finding the first peak
        let histo_inverse = histo.iter().map(|v| -v).collect::<Vec<i32>>();
        let mut peak_finder_inverse = PeakFinder::new(&histo_inverse);
        peak_finder_inverse.with_min_prominence(200);
        let peaks_inverse = peak_finder_inverse.find_peaks();
        let trough = peaks_inverse[0].position.start;

        let n_kmers = count_kmers_in_histo(&histo[trough..]);
        let genome_size = n_kmers as f64 / (peak_genome as f64);

        println!(
            "Peaks: Rare peak at {}, trough at {}, genome peak at {}, genome size {} Mb",
            peak_errors,
            trough,
            peak_genome,
            genome_size / 1e6
        );
    } else {
        println!("Peaks: Multiple genome peaks found");
    }

    // Get the max value in the histogram
    let hist_max = histo[2..].iter().max().unwrap();
    let y_max = (*hist_max as f64 * 1.1) as i32;
    
    // Plot the histogram
    let plot_name = &format!("{out_name}_k{k}_part{chunk_i}.histo.png");
    let file_plot_path = Path::new(&directory).join(plot_name);

    let root_area = BitMapBackend::new(&file_plot_path, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("kmer histogram", ("sans-serif", 40))
        .build_cartesian_2d((0..x_max).into_segmented(), 0..y_max)
        .unwrap();

    ctx.configure_mesh().draw().unwrap();

    ctx.draw_series((0..).zip(histo.iter()).map(|(x, y)| {
        let x0 = SegmentValue::Exact(x);
        let x1 = SegmentValue::Exact(x + 2);
        let mut bar = Rectangle::new([(x0, 0), (x1, *y)], CYAN.filled());
        bar.set_margin(0, 0, 5, 5);
        bar
    }))
    .unwrap();
}

#[inline(always)]
fn string2kmers(seq: &str, k: u32) -> Vec<u64> {
    // kmers are encoded as u64, with two bits per base
    // A = 00, C = 01, G = 10, T = 11

    // Create a mask that has 1s in the last 2*k bits
    let mask: u64 = (1 << (2 * k)) - 1;

    let mut kmers: Vec<u64> = Vec::new();
    let mut frame: u64 = 0; // read the bits for each base into the least significant end of this integer
    let mut revframe: u64 = 0; // read the bits for complement into the least significant end of this integer
    let mut n_valid = 0; // number of valid bases in the frame
    for c in seq.chars() {
        let base = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            'N' => 4,
            _ => 5,
        };

        match base {
            0 | 1 | 2 | 3 => {
                frame = (frame << 2) | base;

                revframe = (revframe >> 2) | ((3 - base) << (2 * (k - 1)));

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
            }
            4 => {
                frame = 0;
                revframe = 0;
                n_valid = 0;
            }
            _ => panic!("Invalid base: {}", c),
        }
    }
    kmers
}

fn count_histogram(kmer_counts: &HashMap<u64, u64>, histo_max: u64) -> Vec<u64> {
    // Create a histogram of counts
    let mut histo: Vec<u64> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for count in kmer_counts.values() {
        if *count <= histo_max {
            histo[*count as usize] += 1;
        } else {
            histo[histo_max as usize + 1] += 1;
        }
    }
    histo
}

fn histogram_string(histo: &[u64]) -> String {
    let mut s = String::new();
    for (i, count) in histo.iter().enumerate() {
        if i > 0 {
            s.push_str(&format!("{} {}\n", i, count));
        }
    }
    s
}

fn get_fastq_stats(input_files: &Vec<String>, max_reads: u64) -> (u64, u64) {
    let mut n_records_all: u64 = 0;
    let mut n_bases_all: u64 = 0;

    'processing_files: for file_name in input_files {
        let mut n_records: u64 = 0;
        let mut n_bases: u64 = 0;

        let path = Path::new(&file_name);
        let file = File::open(path).expect("Ooops.");

        let mut line_n: u64 = 0;

        // From https://github.com/rust-lang/flate2-rs/issues/41#issuecomment-219058833
        // Handle gzip files with multiple blocks
        let mut reader = BufReader::new(file);
        loop {
            //loop over all possible gzip members
            match reader.fill_buf() {
                Ok(b) => {
                    if b.is_empty() {
                        break;
                    }
                }
                Err(e) => panic!("{}", e),
            }

            //decode the next member
            let gz = flate2::bufread::GzDecoder::new(&mut reader);
            let block_reader = BufReader::new(gz);
            for line in block_reader.lines() {
                if line_n % 4 == 1 {
                    let line = line.unwrap();
                    n_records += 1;
                    n_records_all += 1;
                    if max_reads > 0 && n_records_all == max_reads {
                        println!("Reached maximum number of reads: {}", max_reads);
                        println!(
                            "File {}: records {}, bases {}, mean read length {}",
                            file_name,
                            n_records,
                            n_bases,
                            n_bases as f64 / n_records as f64
                        );
                        break 'processing_files;
                    }
                    let line_len = line.len() as u64;
                    n_bases += line_len;
                    n_bases_all += line_len;
                }
                line_n += 1;
            }
        }

        println!(
            "File {}: records {}, bases {}, mean read length {}",
            file_name,
            n_records,
            n_bases,
            n_bases as f64 / n_records as f64
        );
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
    assert!(
        args.k < 32,
        "k must be less than 32 due to use of 64 bit integers to encode kmers"
    );
    assert!(args.k > 0, "k must be greater than 0");
    assert!(args.k % 2 == 1, "k must be odd");
    assert!(args.histo_max > 0, "histo_max must be greater than 0");
    assert!(args.n > 0, "n must be greater than 0");

    // Create an ndarray::Array2 filled with zeros to store the histogram values
    let mut histo_chunks = Array2::<u64>::zeros((args.n, args.histo_max as usize + 2));

    // Create an ndarray::Array1 filled with zeros for bases and records
    let mut n_bases_chunks = Array1::<u64>::zeros(args.n);
    let mut n_records_chunks = Array1::<u64>::zeros(args.n);

    // Get the sequence statistics
    println!("Getting input file statistics...");
    let start = std::time::Instant::now();
    let (n_records_all, n_bases_all) = get_fastq_stats(&args.input, args.max_reads);
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());
    println!(
        "Total records {}, total bases {}, mean read length {}",
        n_records_all,
        n_bases_all,
        n_bases_all as f64 / n_records_all as f64
    );

    // Ingest the data and count kmers
    let mut n_bases_processed: u64 = 0;
    let mut n_records_processed: u64 = 0;

    let chunk_size = n_records_all / args.n as u64;
    let mut kmer_counts: HashMap<u64, u64> = HashMap::new();
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
                Ok(b) => {
                    if b.is_empty() {
                        break;
                    }
                }
                Err(e) => panic!("{}", e),
            }

            //decode the next member
            let gz = flate2::bufread::GzDecoder::new(&mut reader);
            let block_reader = BufReader::new(gz);
            for line in block_reader.lines() {
                if line_n % 4 == 1 {
                    let line_text = line.unwrap();
                    n_records_processed += 1;
                    n_bases_processed += line_text.len() as u64;

                    // Count the kmers
                    let kmers = string2kmers(&line_text, args.k);
                    for kmer in kmers {
                        let count = kmer_counts.entry(kmer).or_insert(0);
                        *count += 1;
                    }

                    // If we've processed enough records, write the output
                    if (n_records_processed % chunk_size == 0) & (n_records_processed > 0) {
                        let histo = count_histogram(&kmer_counts, args.histo_max);
                        let histo_text = histogram_string(&histo);

                        println!(
                            "Cumulative chunk {} stats: {} reads, {} bases, {} unique kmers, {} kmers.", 
                            chunk_i, n_records_processed, n_bases_processed, kmer_counts.len(), kmer_counts.values().sum::<u64>());

                        // Copy the values of histo into the chunk_i row of histo_chunks
                        for (i, v) in histo.iter().enumerate() {
                            histo_chunks[[chunk_i, i]] = *v;
                        }

                        // Copy the values of n_bases_processed and n_records_processed into the chunk_i row of n_bases_chunks and n_records_chunks
                        n_bases_chunks[chunk_i] = n_bases_processed;
                        n_records_chunks[chunk_i] = n_records_processed;

                        // Write the histogram to a file
                        let file_histo_name =
                            &format!("{out_name}_k{k}_part{chunk_i}.histo", k = args.k);
                        let file_histo_path = Path::new(directory).join(file_histo_name);
                        let mut file = File::create(file_histo_path).unwrap();
                        file.write_all(histo_text.as_bytes())
                            .expect("Couldn't write output file");

                        // Write the stats to a file
                        let file_stats_name =
                            &format!("{out_name}_k{k}_part{chunk_i}.stats.tsv", k = args.k);
                        let file_stats_path = Path::new(directory).join(file_stats_name);
                        let mut file_stats = File::create(file_stats_path).unwrap();
                        writeln!(
                            &mut file_stats,
                            "sample	fastq_index	mean_readlen	num_reads	gigabases"
                        )
                        .unwrap();
                        writeln!(
                            &mut file_stats,
                            "{}   {}  {}  {}  {}",
                            out_name,
                            chunk_i,
                            n_bases_processed / n_records_processed,
                            n_records_processed,
                            n_bases_processed as f64 / 1_000_000_000.
                        )
                        .unwrap();

                        // Analyze kmers
                        infer_stats(
                            &histo_chunks,
                            &n_bases_chunks,
                            &n_records_chunks,
                            args.k,
                            chunk_i,
                            directory,
                            out_name,
                        );

                        chunk_i += 1;
                    }

                    if n_records_processed >= (args.n as u64 * chunk_size) {
                        println!(
                            "  Skipping last {} reads after the last full chunk...",
                            n_records_all - n_records_processed
                        );
                        break 'processing_files;
                    }
                }
                line_n += 1;
            }
        }
    }

    println!(
        "Total processed: records {}, bases {}",
        n_records_processed, n_bases_processed
    );
    let stop = std::time::Instant::now();
    println!("Time: {} ms", (stop - start).as_millis());

    drop(kmer_counts);

    // Write histo_chunks to a file
    let histo_file_name = &format!("{out_name}_k{k}_histo.csv", k = args.k);
    let histo_file_path = Path::new(&directory).join(histo_file_name);

    {
        let file = File::create(histo_file_path).unwrap();
        let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
        writer.serialize_array2(&histo_chunks).unwrap();
    }
}

#[cfg(test)]
mod test;
