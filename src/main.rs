#![feature(portable_simd)]

use clap::{Parser, Subcommand};
use colored::Colorize;
use config::Config;
use log::info;
use pretty_env_logger;
use sequence::{SequenceContainer, SequenceOperations};
use std::env;

mod alignment;
mod config;
mod sequence;
mod suffixtree;

#[derive(Debug, Subcommand)]
enum Command {
    Align {
        /// Whether to perform local or global alignment
        #[arg(short, long, default_value = "local")]
        alignment_type: String,
    },
    SuffixTree {
        #[arg(short, long)]
        alphabet_file: String,
    },
}

/// Tool for aligning FASTA sequences with Smith-Waterman or Needleman-Wunsch
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct CliArgs {
    #[clap(subcommand)]
    mode: Command,

    /// Path to the FASTA file containing two sequences to align
    #[arg(short, long)]
    fasta_path: String,

    /// Path to the config file
    #[arg(short, long, default_value = "config.toml")]
    config_path: String,
}

fn main() {
    // parse cli args
    let args = CliArgs::parse();

    // set default log level
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

    // Extract alignment_type from the Command variant
    match &args.mode {
        Command::Align { alignment_type } => {
            // log config
            info!("Using the following values for scoring:");
            info!("Match: {}", config.scores.s_match);
            info!("Mismatch: {}", config.scores.s_mismatch);
            info!("Gap: {}", config.scores.g);
            info!("Opening Gap: {}", config.scores.h);

            info!("Alignment Type: {}", alignment_type);

            info!("{}", "Alignment".bright_green());
            let alignment = alignment::algo::align_sequences(
                &sequence_container,
                &config.scores,
                alignment_type == "local" || alignment_type == "1",
            );

            info!("{}", alignment);
        }
        Command::SuffixTree { alphabet_file } => {
            info!("{}", "Suffix Tree".bright_green());

            let suffix_tree = suffixtree::tree::SuffixTree::new(
                &sequence_container.sequences[0].sequence,
                alphabet_file,
            );

            info!("{}", suffix_tree);
        }
    };
}
