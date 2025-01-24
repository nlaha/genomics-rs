use colored::Colorize;
use log::info;
use pretty_env_logger;
use sequence::{SequenceContainer, SequenceOperations};
use std::env;

mod alignment;
mod sequence;

fn main() {
    // set default log level to trace
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", "info")
    }

    // init logging
    pretty_env_logger::init();

    // print ascii art double helix
    println!(
        "{}",
        format!(
            "{}",
            r#"
        GENOME-RS
        -. .-.   .-. .-.   .-. .-.   .  
        ||\|||\ /|||\|||\ /|||\|||\ /|
        |/ \|||\|||/ \|||\|||/ \|||\||
        ~   `-~ `-`   `-~ `-`   `-~ `-
    "#
        )
        .bright_blue()
    );

    let mut sequence_container: SequenceContainer = SequenceContainer {
        sequences: Vec::new(),
    };

    sequence_container.from_fasta("test_data/Opsin1_colorblindness_gene.fasta");

    let aligned_sequences = alignment::global_alignment(sequence_container);

    info!("{}", aligned_sequences);
}
