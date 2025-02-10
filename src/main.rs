#![feature(portable_simd)]

use clap::Parser;
use colored::Colorize;
use config::Config;
use log::info;
use pretty_env_logger;
use sequence::{SequenceContainer, SequenceOperations};
use std::env;

mod alignment;
mod config;
mod sequence;

/// Tool for aligning FASTA sequences with Smith-Waterman or Needleman-Wunsch
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct CliArgs {
    /// Path to the FASTA file containing two sequences to align
    #[arg(short, long)]
    fasta_path: String,

    /// Whether to perform local or global alignment
    #[arg(short, long, default_value = "local")]
    alignment_type: String,

    /// Path to the config file
    #[arg(short, long, default_value = "config.toml")]
    config_path: String,
}

fn main() {
    // parse cli args
    let args = CliArgs::parse();

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
        GENOMICS-RS
        -. .-.   .-. .-.   .-. .-.   .  
        ||\|||\ /|||\|||\ /|||\|||\ /|
        |/ \|||\|||/ \|||\|||/ \|||\||
        ~   `-~ `-`   `-~ `-`   `-~ `-
    "#
        )
        .bright_blue()
    );

    // load config
    let config: Config = config::get_config(args.config_path.as_str());

    let mut sequence_container: SequenceContainer = SequenceContainer {
        sequences: Vec::new(),
    };

    info!("Loading sequences from {}", args.fasta_path);

    sequence_container.from_fasta(args.fasta_path.as_str());

    info!("Alignment Type: {}", args.alignment_type);

    info!("{}", "Alignment".bright_green());
    let alignment = alignment::algo::align_sequences(
        &sequence_container,
        &config.scores,
        args.alignment_type == "local" || args.alignment_type == "1",
    );

    info!("{}", alignment);
}
