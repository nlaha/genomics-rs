use log::{info, trace, warn};
use std::collections::HashMap;

mod fasta;

fn main() {
    let sequences: HashMap<String, String> = fasta::load_sequences("test_data/test1.fasta");

    // print out sequences
    for (name, seq) in sequences {
        info!("{}: {}", name, seq);
    }
}
