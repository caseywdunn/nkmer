use super::*;

use std::collections::HashSet;

// Function that takes a DNA sequence as a string and returns reverse complement
fn revcomp(seq: &str) -> String {
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
fn test_tests() {
    assert_eq!(2 + 2, 4);
}

#[test]
fn test_revcomp() {
    let seq = String::from("ACGNTT");
    let rc = revcomp(&seq);
    assert_eq!(rc, "AANCGT");
}

#[test]
fn test_strings2kmer() {
    let seq = String::from("ACGNTT");
    let kmer = string2kmers(&seq, 3);

    // The only valid 3-mer in this sequence is ACG
    // A = 00, C = 01, G = 10, T = 11
    // ACG = 000110
    // CGT = 011011 , need to check reverse complement
    // It is the smallest that will be retained
    assert_eq!(
        kmer[0],
        0b0000000000000000000000000000000000000000000000000000000000000110
    );

    // Check reverse complement of the sequence, which should be the same
    let seq_rc = revcomp(&seq);
    let kmer_rc = string2kmers(&seq_rc, 3);
    assert_eq!(
        kmer_rc[0],
        0b0000000000000000000000000000000000000000000000000000000000000110
    );
}

#[test]
fn test_strings2kmer_revcomp() {
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

#[test]
fn test_count_histogram() {
    let kmer_counts: HashMap<u64, u64> = HashMap::from([
        (0, 2),
        (1, 2),
        (2, 2),
        (3, 2),
        (5, 3),
        (6, 3),
        (4, 4),
        (9, 5),
        (10, 5),
        (11, 5),
        (12, 5),
        (13, 5),
        (7, 20),
        (8, 1000),
    ]);

    let histo_counts = count_histogram(&kmer_counts, 10);
    assert_eq!(histo_counts, [0, 0, 4, 2, 1, 5, 0, 0, 0, 0, 0, 2]);
}
