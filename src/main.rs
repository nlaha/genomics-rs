#![feature(portable_simd)]

use clap::{Parser, Subcommand};
use colored::{control, Colorize};
use comparison::display::print_similarity_matrix;
use config::Config;
use log::info;
use ndarray::parallel::prelude::*;
use ndarray::Axis;
use ndarray::{Array2, ShapeBuilder};
use pretty_env_logger;
use rayon::iter::ParallelBridge;
use sequence::{SequenceContainer, SequenceOperations};
use std::time::Instant;
use std::{
    env,
    fs::{self, OpenOptions},
    io::{self, Write},
};

mod alignment;
mod comparison;
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

        /// Whether to compute suffix links
        #[arg(long, default_value_t = false)]
        suffix_links: bool,

        /// Number of threads to use
        #[arg(long, default_value_t = 1)]
        threads: usize,
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
            let (alignment_table, _) = alignment::algo::alignment_table(
                &sequence_container,
                &config.scores,
                is_local,
                false,
            );
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
            suffix_links,
            threads,
        } => {
            info!("MODE: {}", "Compare".bright_green().bold());
            info!("Alphabet file: {}", alphabet_file);
            info!("Suffix links: {}", suffix_links);
            info!("FASTA directory: {}", fasta_dir);

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
                suffix_tree.insert_string(&sequence.sequence, *suffix_links);
            }

            let num_sequences = sequence_container.sequences.len();
            info!("Number of sequences: {}", num_sequences);

            let mut similarity_matrix: Array2<(usize, usize, usize)> =
                Array2::from_shape_fn((num_sequences, num_sequences).f(), |_| (0, 0, 0));

            // configure thread pool for concurrent processing
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.to_owned())
                .build_global()
                .unwrap();

            // start timer
            let before_compare = Instant::now();

            // get lcs for each pair of sequences in the tree
            similarity_matrix
                .axis_iter_mut(Axis(0))
                .into_par_iter()
                .enumerate()
                .for_each(|(j, mut row)| {
                    row.indexed_iter_mut()
                        .par_bridge()
                        .for_each(|(i, similarity_score)| {
                            // skip if i > j
                            if i > j {
                                return;
                            }

                            // get initial LCS on entire sequences
                            let (start_i, start_j, lcs_length) = suffix_tree.get_lcs(i, j);

                            let get_matches = |s1: &String, s2: &String| {
                                let mut st = suffixtree::tree::SuffixTree::new(
                                    alphabet_file,
                                    s1.len() + s2.len(),
                                );

                                st.insert_string(s1, true);
                                st.insert_string(s2, true);

                                let (st_i, st_j, len) = st.get_lcs(0, 1);

                                return (len, st_i, st_j, s1.clone(), s2.clone());
                            };

                            // we're going to emulate recursion here
                            let mut stack: Vec<(usize, usize, usize, String, String)> = Vec::new();
                            stack.push((
                                lcs_length,
                                start_i,
                                start_j,
                                sequence_container.sequences[i].sequence.clone(),
                                sequence_container.sequences[j].sequence.clone(),
                            ));

                            // basically do a DFS
                            let mut score = 0;
                            while !stack.is_empty() {
                                let (lcs_length, st_i, st_j, s1, s2) = stack.pop().unwrap();

                                let prefix_i = &s1[..st_i].to_string();
                                let prefix_j = &s2[..st_j].to_string();
                                let suffix_i = &s1[st_i + lcs_length..].to_string();
                                let suffix_j = &s2[st_j + lcs_length..].to_string();

                                // push children
                                if lcs_length > 0 {
                                    stack.push(get_matches(prefix_i, prefix_j));
                                    stack.push(get_matches(suffix_i, suffix_j));
                                }

                                // process
                                score += lcs_length;
                            }

                            *similarity_score = (
                                score,
                                sequence_container.sequences[i].sequence.len(),
                                sequence_container.sequences[j].sequence.len(),
                            );
                        });
                });

            let after_compare = Instant::now();

            let elapsed = after_compare.duration_since(before_compare).as_micros();
            let elapsed_millis = after_compare.duration_since(before_compare).as_millis();
            info!(
                "[FindPath] Time taken to compare: {} us ({} ms)",
                elapsed, elapsed_millis
            );

            print_similarity_matrix(&similarity_matrix);
        }
    };

    Ok(())
}
