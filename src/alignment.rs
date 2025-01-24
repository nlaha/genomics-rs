use std::{cmp, fmt::Display};

use colored::Colorize;
use log::{debug, info, warn};
use ndarray::{Array2, ShapeBuilder};

use crate::sequence::{Sequence, SequenceContainer, SequenceOperations};

// set up parameters
const MATCH_SCORE: i64 = 1;
const MISMATCH_SCORE: i64 = -2;
const GAP_PENALTY: i64 = -2;
const H: i64 = -5;

// we need to offset our 'negative infinity' to avoid overflow errors
const NEGATIVE_INF: i64 = i64::MIN + (GAP_PENALTY * GAP_PENALTY) + (H * H);

#[derive(Debug, Clone, Copy)]
struct AlignmentCell {
    insert_score: i64,
    delete_score: i64,
    sub_score: i64,
    is_match: bool,
}

pub trait ComputeScore {
    fn substitution(&self) -> i64;
    fn score_max(&self, i_mod: i64, s_mod: i64, d_mod: i64) -> i64;
}

impl ComputeScore for AlignmentCell {
    fn substitution(&self) -> i64 {
        if self.is_match {
            MATCH_SCORE
        } else {
            MISMATCH_SCORE
        }
    }

    fn score_max(&self, i_mod: i64, s_mod: i64, d_mod: i64) -> i64 {
        cmp::max(
            self.insert_score + i_mod,
            cmp::max(self.sub_score + s_mod, self.delete_score + d_mod),
        )
    }
}

#[derive(Clone, Copy, Debug)]
#[repr(u8)]
pub enum AlignmentChoice {
    Match,
    Mismatch,
    Insert,
    Delete,
}

pub struct AlignedSequences {
    pub s1: Sequence,
    pub s2: Sequence,
    pub alignment: Vec<(AlignmentChoice, usize, usize)>,
}

impl Display for AlignedSequences {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s1_out: String = String::new();
        let mut align_out: String = String::new();
        let mut s2_out: String = String::new();

        let mut s1_idx = 0;
        let mut s2_idx = 0;
        let mut alignement_iter = self.alignment.iter().rev();

        while s1_idx < self.s1.sequence.len() || s2_idx < self.s2.sequence.len() {
            let choice = alignement_iter.next();

            let s1_base = self.s1.sequence.as_bytes()[s1_idx] as char;
            let s2_base = self.s2.sequence.as_bytes()[s2_idx] as char;

            // print s1 first
            match choice {
                Some((AlignmentChoice::Insert, _, _)) => s1_out.push('-'),
                _ => {
                    s1_out.push(s1_base);
                    s1_idx += 1;
                }
            }

            // then print the alignment row
            match choice {
                Some((AlignmentChoice::Match, _, _)) => align_out.push('|'),
                Some((AlignmentChoice::Mismatch, _, _)) => align_out.push('x'),
                Some((AlignmentChoice::Insert, _, _)) => align_out.push(' '),
                Some((AlignmentChoice::Delete, _, _)) => align_out.push(' '),
                None => align_out.push(' '),
            }

            // then print s2
            match choice {
                Some((AlignmentChoice::Delete, _, _)) => s2_out.push('-'),
                _ => {
                    s2_out.push(s2_base);
                    s2_idx += 1;
                }
            }
        }

        writeln!(f, "\n\n{}\n{}\n{}", s1_out, align_out, s2_out)
    }
}

/// Perform global alignment of two sequences
pub fn global_alignment(sequence_container: SequenceContainer) -> AlignedSequences {
    // if more than two sequences are found
    if sequence_container.sequences.len() > 2 {
        warn!("More than two sequences found. Only the first two will be used.");
    }

    let s1_len = sequence_container.sequences[0].sequence.len();
    let s2_len = sequence_container.sequences[1].sequence.len();

    // create 2D array to store dynamic programming results
    let mut sequence_table: Array2<AlignmentCell> = Array2::from_elem(
        (s1_len, s2_len).f(),
        AlignmentCell {
            insert_score: 0,
            delete_score: 0,
            sub_score: 0,
            is_match: false,
        },
    );

    // log
    info!("Sequence table shape: {:?}", sequence_table.shape());
    info!(
        "Sequence table size (KB): {}",
        std::mem::size_of_val(&sequence_table) as f64 / 1024.0
    );

    // start timer
    let start_table_init: std::time::Instant = std::time::Instant::now();

    // iterate through alignment table and fill in scores
    for i in 0..s1_len {
        for j in 0..s2_len {
            sequence_table[[i, j]] = match (i, j) {
                // origin
                (0, 0) => AlignmentCell {
                    insert_score: 0,
                    delete_score: 0,
                    sub_score: 0,
                    is_match: sequence_container.is_match(i, j),
                },
                // column
                (i, 0) => AlignmentCell {
                    insert_score: NEGATIVE_INF,
                    delete_score: H + (i as i64 * GAP_PENALTY),
                    sub_score: NEGATIVE_INF,
                    is_match: sequence_container.is_match(i, j),
                },
                // row
                (0, j) => AlignmentCell {
                    insert_score: H + (j as i64 * GAP_PENALTY),
                    delete_score: NEGATIVE_INF,
                    sub_score: NEGATIVE_INF,
                    is_match: sequence_container.is_match(i, j),
                },
                (i, j) => {
                    // recurrence relationships for S, D and I
                    let top_left = sequence_table[[i - 1, j - 1]];
                    let left = sequence_table[[i, j - 1]];
                    let top = sequence_table[[i - 1, j]];

                    AlignmentCell {
                        insert_score: left.score_max(GAP_PENALTY, H + GAP_PENALTY, H + GAP_PENALTY),
                        delete_score: top.score_max(H + GAP_PENALTY, H + GAP_PENALTY, GAP_PENALTY),
                        sub_score: top_left.substitution() + top_left.score_max(0, 0, 0),
                        is_match: sequence_container.is_match(i, j),
                    }
                }
            };
        }
    }

    let end_table_init: std::time::Instant = std::time::Instant::now();

    info!(
        "Table initialization complete, time taken: {}us",
        (end_table_init - start_table_init).as_micros()
    );

    // log score
    info!(
        "Optimal alignment found, score: {}",
        sequence_table[[s1_len - 1, s2_len - 1]].score_max(0, 0, 0)
    );

    // optimal retrace algorithm
    let start_retrace: std::time::Instant = std::time::Instant::now();

    let (mut i, mut j) = (s1_len - 1, s2_len - 1);
    let mut aligned_sequences: AlignedSequences = AlignedSequences {
        s1: sequence_container.sequences[0].clone(),
        s2: sequence_container.sequences[1].clone(),
        alignment: Vec::new(),
    };

    loop {
        let cell = sequence_table[[i, j]];
        let max = cell.score_max(0, 0, 0);

        let (i_opt, j_opt) = match max {
            // top_left move
            val if val == cell.sub_score => {
                if cell.is_match {
                    aligned_sequences
                        .alignment
                        .push((AlignmentChoice::Match, i, j));
                    debug!("Match found at ({}, {})", i, j);
                } else {
                    aligned_sequences
                        .alignment
                        .push((AlignmentChoice::Mismatch, i, j));
                    debug!("Mismatch found at ({}, {})", i, j);
                }
                (i.checked_sub(1), j.checked_sub(1))
            }
            // left move
            val if val == cell.insert_score => {
                aligned_sequences
                    .alignment
                    .push((AlignmentChoice::Insert, i, j));
                debug!("Insert found at ({}, {})", i, j);
                (Some(i), j.checked_sub(1))
            }
            // up move
            val if val == cell.delete_score => {
                aligned_sequences
                    .alignment
                    .push((AlignmentChoice::Delete, i, j));
                debug!("Delete found at ({}, {})", i, j);
                (i.checked_sub(1), Some(j))
            }
            _ => panic!("Unexpected score during retrace"),
        };

        (i, j) = match (i_opt, j_opt) {
            (Some(i), Some(j)) => (i, j),
            (None, None) => break,
            (Some(i), None) => (i, 0),
            (None, Some(j)) => (0, j),
        };
    }

    let end_retrace: std::time::Instant = std::time::Instant::now();

    info!(
        "Retrace complete, time taken: {}us",
        (end_retrace - start_retrace).as_micros()
    );

    info!(
        "Retrace alignment size: {}",
        aligned_sequences.alignment.len()
    );

    info!("Computing sequence table visualization...");
    print_sequence_table(sequence_container, sequence_table, &aligned_sequences, 100);

    return aligned_sequences;
}

/// Prints the sequence table with alignment path for visualization
fn print_sequence_table(
    sequence_container: SequenceContainer,
    sequence_table: ndarray::ArrayBase<ndarray::OwnedRepr<AlignmentCell>, ndarray::Dim<[usize; 2]>>,
    aligned_sequences: &AlignedSequences,
    downsample: usize,
) {
    let s1_len = sequence_container.sequences[0].sequence.len();
    let s2_len = sequence_container.sequences[1].sequence.len();

    println!("\nSequence Table (S1 columns, S2 rows):\n");

    // print the columns (bases of the second sequence)
    print!(" ",);
    for i in (0..s2_len).step_by(downsample) {
        print!(
            "{}",
            sequence_container.sequences[1]
                .sequence
                .chars()
                .nth(i)
                .unwrap()
        );
    }
    println!();

    // display dynamic programming table
    for i in (0..s1_len).step_by(downsample) {
        // print sequence base
        print!(
            "{}",
            sequence_container.sequences[0]
                .sequence
                .chars()
                .nth(i)
                .unwrap()
        );
        for j in (0..s2_len).step_by(downsample) {
            // pick character from BRIGHT_STRING based on the score
            let brightness_char = {
                let score: i64 = sequence_table[[i, j]].score_max(0, 0, 0);
                // average score if downsample is greater than 1
                let score = if downsample > 1 {
                    let mut sum = 0;
                    for k in 0..downsample {
                        sum += sequence_table[[i, j + k]].score_max(0, 0, 0);
                    }
                    sum / downsample as i64
                } else {
                    score
                };

                if score > 5 || sequence_table[[i, j]].is_match {
                    if sequence_table[[i, j]].is_match {
                        format!("M").dimmed()
                    } else {
                        format!("X").dimmed()
                    }
                } else {
                    format!(".").dimmed()
                }
            };

            // if this coordinate it in the alignment path, print a different character
            let alignment_choice = aligned_sequences
                .alignment
                .iter()
                .find(|(_choice, x, y)| *x == i && *y == j);

            match alignment_choice {
                Some((AlignmentChoice::Match, _, _)) => print!("{}", "M".green()),
                Some((AlignmentChoice::Mismatch, _, _)) => print!("{}", "X".red()),
                Some((AlignmentChoice::Insert, _, _)) => print!("{}", "I".blue()),
                Some((AlignmentChoice::Delete, _, _)) => print!("{}", "D".cyan()),
                None => print!("{}", brightness_char),
            }
        }
        println!();
    }
}
