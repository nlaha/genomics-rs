#![feature(portable_simd)]

use clap::{Parser, Subcommand};
use colored::Colorize;
use config::Config;
use log::info;
use pretty_env_logger;
use sequence::{SequenceContainer, SequenceOperations};
use std::{
    env,
    fs::{self, OpenOptions},
    io::{self, Write},
};

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

        /// Path to the FASTA file containing two sequences to align
        #[arg(short, long)]
        fasta_path: String,
    },
    SuffixTree {
        /// Path to the alphabet file
        #[arg(short, long)]
        alphabet_file: String,

        /// Whether to compute suffix links
        #[arg(long, default_value_t = false)]
        suffix_links: bool,

        /// Whether to compute stats
        #[arg(long, default_value_t = false)]
        stats: bool,

        /// Path to the FASTA file containing two sequences to align
        #[arg(short, long)]
        fasta_path: String,
    },
    Compare {
        /// Path to the alphabet file
        #[arg(short, long)]
        alphabet_file: String,

        /// Path to directory containing FASTA files
        #[arg(short, long)]
        fasta_dir: String,
    },
}

/// Tool for aligning FASTA sequences with Smith-Waterman or Needleman-Wunsch
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct CliArgs {
    #[clap(subcommand)]
    mode: Command,

    /// Path to the config file
    #[arg(short, long, default_value = "config.toml")]
    config_path: String,
}

fn main() -> io::Result<()> {
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

    // Extract alignment_type from the Command variant
    match &args.mode {
        Command::Align {
            alignment_type,
            fasta_path,
        } => {
            info!("MODE: {}", "Alignment".bright_yellow().bold());

            info!("Loading sequences from {}", fasta_path);
            sequence_container.from_fasta(fasta_path.as_str());

            // log config
            info!("Using the following values for scoring:");
            info!("Match: {}", config.scores.s_match);
            info!("Mismatch: {}", config.scores.s_mismatch);
            info!("Gap: {}", config.scores.g);
            info!("Opening Gap: {}", config.scores.h);

            info!("Alignment Type: {}", alignment_type);

            info!("{}", "Alignment".bright_green());
            let is_local = alignment_type == "local" || alignment_type == "1";
            let alignment_table =
                alignment::algo::alignment_table(&sequence_container, &config.scores, is_local);
            let alignment =
                alignment::algo::retrace(&sequence_container, alignment_table, is_local);

            info!("{}", alignment);
        }
        Command::SuffixTree {
            alphabet_file,
            suffix_links,
            stats,
            fasta_path,
        } => {
            info!("MODE: {}", "Suffix Tree".bright_green().bold());
            info!("Suffix links: {}", suffix_links);

            info!("Loading sequences from {}", fasta_path);
            sequence_container.from_fasta(fasta_path.as_str());

            let mut suffix_tree = suffixtree::tree::SuffixTree::new(
                alphabet_file,
                sequence_container.sequences[0].sequence.len(),
            );

            suffix_tree.insert_string(&sequence_container.sequences[0].sequence, *suffix_links);

            if *stats {
                suffix_tree.compute_stats(0);

                // delete bwt file if it exists
                let bwt_path = format!(
                    "BWT_out/{}_bwt.txt",
                    fasta_path
                        .split("/")
                        .last()
                        .unwrap()
                        .split("\\")
                        .last()
                        .unwrap()
                        .replace(".fasta", "")
                );

                info!("BWT Path: {}", bwt_path);

                if let Err(e) = std::fs::remove_file(bwt_path.clone()) {
                    eprintln!("Couldn't delete file: {}", e);
                }

                let mut out_file = OpenOptions::new()
                    .create(true)
                    .write(true)
                    .append(true)
                    .open(bwt_path)
                    .unwrap();

                // write bwt to file
                for c in suffix_tree.stats.bwt.chars() {
                    if let Err(e) = writeln!(out_file, "{}", c) {
                        eprintln!("Couldn't write to file: {}", e);
                    }
                }

                info!("{}", suffix_tree);
            }
        }
        Command::Compare {
            alphabet_file,
            fasta_dir,
        } => {
            info!("MODE: {}", "Compare".bright_green().bold());
            info!("Alphabet file: {}", alphabet_file);

            // get all .fasta files in the fasta_dir
            // and load them into the sequence container
            for file in fs::read_dir(fasta_dir)? {
                let file = file?;
                if file.path().extension().unwrap() != "fasta" {
                    continue;
                }
                sequence_container.from_fasta(file.path().to_str().unwrap());
            }

            let mut suffix_tree = suffixtree::tree::SuffixTree::new(
                alphabet_file,
                sequence_container.sequences[0].sequence.len(),
            );

            for sequence in sequence_container.sequences.iter() {
                suffix_tree.insert_string(&sequence.sequence, true);
            }

            // get lcs for each pair of sequences in the tree
            for i in 0..sequence_container.sequences.len() {
                for j in i + 1..sequence_container.sequences.len() {
                    let (start_i, start_j, length) = suffix_tree.get_lcs(i, j);
                    // get the prefixes from i and j
                    let prefix_i = &sequence_container.sequences[i].sequence[..start_i];
                    let prefix_j = &sequence_container.sequences[j].sequence[..start_j];

                    let prefix_container = SequenceContainer::from_prefixes(prefix_i, prefix_j);

                    let alignment_table =
                        alignment::algo::alignment_table(&prefix_container, &config.scores, false);
                }
            }

            info!("{}", suffix_tree);
        }
    };

    Ok(())
}
