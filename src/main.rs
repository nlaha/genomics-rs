#![feature(portable_simd)]

use colored::Colorize;
use config::Config;
use log::info;
use pretty_env_logger;
use sequence::{SequenceContainer, SequenceOperations};
use std::env;

mod alignment;
mod config;
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

    // load config
    let config: Config = config::get_config("config.toml");

    let mut sequence_container: SequenceContainer = SequenceContainer {
        sequences: Vec::new(),
    };

    sequence_container.from_fasta(config.fasta.path.as_str());

    info!("{}", "Global Alignment".bright_green());
    let global_alignment = alignment::align_sequences(&sequence_container, &config.scores, false);

    info!("{}", global_alignment);

    info!("{}", "Local Alignment".bright_purple());
    let local_alignment = alignment::align_sequences(&sequence_container, &config.scores, true);

    info!("{}", local_alignment);
}
