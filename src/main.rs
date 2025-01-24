use log::{info, trace, warn};
use ndarray::{Array, Array2, ShapeBuilder};
use num;
use pretty_env_logger;
use std::collections::HashMap;

mod fasta;

#[derive(Debug, Clone)]
struct AlignmentCell {
    insert_score: i32,
    delete_score: i32,
    sub_score: i32,
}

fn main() {
    // init logging
    pretty_env_logger::init();

    // print ascii art double helix
    // make it green with terminal escape codes
    println!(
        r#"
        \x1b[38;5;92m
        GENOME-RS
        -. .-.   .-. .-.   .-. .-.   .  
        ||\|||\ /|||\|||\ /|||\|||\ /|
        |/ \|||\|||/ \|||\|||/ \|||\||
        ~   `-~ `-`   `-~ `-`   `-~ `-
    "#
    );

    let sequences: Vec<fasta::FastaSequence> =
        fasta::load_sequences("test_data/Opsin1_colorblindness_gene.fasta");

    // if more than two sequences are found
    if sequences.len() > 2 {
        warn!("More than two sequences found. Only the first two will be used.");
    }

    // create 2D array to store dynamic programming results
    let sequence_table: Array2<AlignmentCell> = Array::from_elem(
        (
            sequences[0].sequence.len() + 1,
            sequences[1].sequence.len() + 1,
        )
            .f(),
        AlignmentCell {
            insert_score: 0,
            delete_score: 0,
            sub_score: 0,
        },
    );

    // log
    info!("Sequence table shape: {:?}", sequence_table.shape());
    info!(
        "Sequence table size (KB): {}",
        std::mem::size_of_val(&sequence_table) as f64 / 1024.0
    );
}
