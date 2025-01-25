#![feature(portable_simd)]

use colored::Colorize;
use config::Config;
use log::info;
use pretty_env_logger;
use sequence::{SequenceContainer, SequenceOperations};
use std::env;

mod alignment;
mod sequence;
mod config;

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

    // load config
    let config: Config = config::get_config("config.toml");

    let mut sequence_container: SequenceContainer = SequenceContainer {
        sequences: Vec::new(),
    };

    sequence_container.from_fasta(config.fasta.path.as_str());

    let aligned_sequences = alignment::global_alignment(sequence_container, config.scores);

    info!("{}", aligned_sequences);
}
