use std::{cmp, fmt::Display};

use colored::Colorize;
use log::{debug, info, warn};
use ndarray::{Array2, ShapeBuilder};
use spinoff::{spinners, Color, Spinner};

use crate::sequence::{Sequence, SequenceContainer, SequenceOperations};

// set up parameters
const MATCH_SCORE: i64 = 1;
const MISMATCH_SCORE: i64 = -2;
const GAP_PENALTY: i64 = -2;
const H: i64 = -5;

// we need to offset our 'negative infinity' to avoid overflow errors
const NEGATIVE_INF: i64 = i64::MIN + (GAP_PENALTY * GAP_PENALTY) + (H * H);

const DISP_MAX_WIDTH: usize = 200;

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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum AlignmentChoice {
    Match,
    Mismatch,
    Insert,
    Delete,
    OpenInsert,
    OpenDelete,
}

pub struct AlignedSequences {
    pub s1: Sequence,
    pub s2: Sequence,
    pub alignment: Vec<(AlignmentChoice, usize, usize)>,

    // stats
    pub score: i64,
    pub matches: usize,
    pub mismatches: usize,
    pub gap_extensions: usize,
    pub opening_gaps: usize,
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
        sequence_table.len() * std::mem::size_of::<AlignmentCell>() / 1024
    );

    let mut sequence_table_spinner = Spinner::new(
        spinners::Dots,
        "Computing sequence table...",
        Color::Magenta,
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

    sequence_table_spinner.success("Sequence table computed");

    info!(
        "Table initialization complete, time taken: {}us",
        (end_table_init - start_table_init).as_micros()
    );

    // log score
    info!(
        "Optimal alignment found, score: {}",
        sequence_table[[s1_len - 1, s2_len - 1]].score_max(0, 0, 0)
    );

    let aligned_sequences = retrace(sequence_container, sequence_table);

    return aligned_sequences;
}

/// Performs a retrace of the optimal alignment path from the sequence table
fn retrace(
    sequence_container: SequenceContainer,
    sequence_table: ndarray::ArrayBase<ndarray::OwnedRepr<AlignmentCell>, ndarray::Dim<[usize; 2]>>,
) -> AlignedSequences {
    let mut retrace_spinner = Spinner::new(
        spinners::Dots,
        "Retracing optimal alignment...",
        Color::Magenta,
    );

    let s1_len = sequence_container.sequences[0].sequence.len();
    let s2_len = sequence_container.sequences[1].sequence.len();

    // optimal retrace algorithm
    let start_retrace: std::time::Instant = std::time::Instant::now();

    let (mut i, mut j) = (s1_len - 1, s2_len - 1);
    let mut aligned_sequences: AlignedSequences = AlignedSequences {
        s1: sequence_container.sequences[0].clone(),
        s2: sequence_container.sequences[1].clone(),
        alignment: Vec::new(),
        score: sequence_table[[i, j]].score_max(0, 0, 0),
        matches: 0,
        mismatches: 0,
        gap_extensions: 0,
        opening_gaps: 0,
    };

    let mut last_choice = AlignmentChoice::Match;
    loop {
        let cell = sequence_table[[i, j]];
        let max = cell.score_max(0, 0, 0);

        let (i_opt, j_opt) = match max {
            // top_left move
            val if val == cell.sub_score => {
                if cell.is_match {
                    last_choice = AlignmentChoice::Match;
                    aligned_sequences.matches += 1;
                    aligned_sequences
                        .alignment
                        .push((AlignmentChoice::Match, i, j));
                    debug!("Match found at ({}, {})", i, j);
                } else {
                    last_choice = AlignmentChoice::Mismatch;
                    aligned_sequences.mismatches += 1;
                    aligned_sequences
                        .alignment
                        .push((AlignmentChoice::Mismatch, i, j));
                    debug!("Mismatch found at ({}, {})", i, j);
                }
                (i.checked_sub(1), j.checked_sub(1))
            }
            // left move
            val if val == cell.insert_score => {
                let choice = if last_choice == AlignmentChoice::Insert {
                    aligned_sequences.gap_extensions += 1;
                    AlignmentChoice::Insert
                } else {
                    aligned_sequences.opening_gaps += 1;
                    AlignmentChoice::OpenInsert
                };

                aligned_sequences.alignment.push((choice, i, j));
                debug!("Insert found at ({}, {})", i, j);
                last_choice = AlignmentChoice::Insert;
                (Some(i), j.checked_sub(1))
            }
            // up move
            val if val == cell.delete_score => {
                let choice = if last_choice == AlignmentChoice::Delete {
                    aligned_sequences.gap_extensions += 1;
                    AlignmentChoice::Delete
                } else {
                    aligned_sequences.opening_gaps += 1;
                    AlignmentChoice::OpenDelete
                };

                aligned_sequences.alignment.push((choice, i, j));
                debug!("Delete found at ({}, {})", i, j);
                last_choice = AlignmentChoice::Delete;
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

    retrace_spinner.success("Retrace complete");

    info!(
        "Retrace complete, time taken: {}us",
        (end_retrace - start_retrace).as_micros()
    );

    info!(
        "Retrace alignment size: {}",
        aligned_sequences.alignment.len()
    );

    // skip visualization if the table is too large
    if sequence_table.shape()[0] < DISP_MAX_WIDTH && sequence_table.shape()[1] < DISP_MAX_WIDTH * 10
    {
        info!("Computing sequence table visualization...");
        print_sequence_table(sequence_container, sequence_table, &aligned_sequences);
    } else {
        warn!("Sequence table too large to visualize");
    }

    aligned_sequences
}

impl Display for AlignedSequences {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // print original and aligned sequences if they're not too long
        if self.s1.sequence.len() <= DISP_MAX_WIDTH && self.s2.sequence.len() <= DISP_MAX_WIDTH {
            info!("Original Sequences:");
            info!("{}", self.s1);
            info!("{}", self.s2);
        } else {
            warn!("Sequences are too long to display.");
        }

        let mut s1_out: String = String::new();
        let mut align_out: String = String::new();
        let mut s2_out: String = String::new();

        let mut s1_idx = 0;
        let mut s2_idx = 0;
        let mut alignement_iter = self.alignment.iter().rev();

        let mut horizontal_len = 0;
        let mut align_idx = 0;

        while align_idx < self.alignment.len() {
            let choice = alignement_iter.next();

            if horizontal_len > DISP_MAX_WIDTH {
                // print the start/end index we're on
                writeln!(f, "\n\n{}-{}:\n", align_idx - DISP_MAX_WIDTH, align_idx)?;
                // print the current chunk
                writeln!(f, "{}\n{}\n{}", s1_out, align_out, s2_out)?;
                // reset the output strings
                s1_out.clear();
                align_out.clear();
                s2_out.clear();
                horizontal_len = 0;
            }

            // print s1 first
            match choice {
                Some((AlignmentChoice::Insert | AlignmentChoice::OpenInsert, _, _)) => {
                    s1_out.push('-')
                }
                _ => {
                    s1_out.push(self.s1.sequence.as_bytes()[s1_idx] as char);
                    s1_idx += 1;
                }
            }

            // then print the alignment row
            match choice {
                Some((AlignmentChoice::Match, _, _)) => align_out.push('|'),
                Some((AlignmentChoice::Mismatch, _, _)) => align_out.push('x'),
                Some((AlignmentChoice::Insert, _, _)) => align_out.push(' '),
                Some((AlignmentChoice::Delete, _, _)) => align_out.push(' '),
                Some((AlignmentChoice::OpenInsert | AlignmentChoice::OpenDelete, _, _)) => {
                    align_out.push('%')
                }
                None => align_out.push(' '),
            }

            // then print s2
            match choice {
                Some((AlignmentChoice::Delete | AlignmentChoice::OpenDelete, _, _)) => {
                    s2_out.push('-')
                }
                _ => {
                    s2_out.push(self.s2.sequence.as_bytes()[s2_idx] as char);
                    s2_idx += 1;
                }
            }

            horizontal_len += 1;
            align_idx += 1;
        }

        writeln!(f, "\n\n{}-{}:\n", align_idx - s1_out.len(), align_idx)?;
        writeln!(f, "{}\n{}\n{}", s1_out, align_out, s2_out)?;

        // print stats report
        writeln!(f, "\n\nAlignment Score: {}", self.score)?;
        writeln!(
            f,
            "Matches: {}/{} ({:.2}%)",
            self.matches,
            align_idx,
            (self.matches as f64 / align_idx as f64) * 100.0
        )?;
        writeln!(
            f,
            "Mismatches: {}/{} ({:.2}%)",
            self.mismatches,
            align_idx,
            (self.mismatches as f64 / align_idx as f64) * 100.0
        )?;
        writeln!(
            f,
            "Gap Extensions: {}/{} ({:.2}%)",
            self.gap_extensions,
            align_idx,
            (self.gap_extensions as f64 / align_idx as f64) * 100.0
        )?;
        writeln!(
            f,
            "Opening Gaps: {}/{} ({:.2}%)",
            self.opening_gaps,
            align_idx,
            (self.opening_gaps as f64 / align_idx as f64) * 100.0
        )
    }
}

/// Prints the sequence table with alignment path for visualization
fn print_sequence_table(
    sequence_container: SequenceContainer,
    sequence_table: ndarray::ArrayBase<ndarray::OwnedRepr<AlignmentCell>, ndarray::Dim<[usize; 2]>>,
    aligned_sequences: &AlignedSequences,
) {
    let s1_len = sequence_container.sequences[0].sequence.len();
    let s2_len = sequence_container.sequences[1].sequence.len();

    println!("\nSequence Table (S1 columns, S2 rows):\n");

    // print the columns (bases of the second sequence)
    print!(" ",);
    for i in 0..s2_len {
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
    for i in 0..s1_len {
        // print sequence base
        print!(
            "{}",
            sequence_container.sequences[0]
                .sequence
                .chars()
                .nth(i)
                .unwrap()
        );
        for j in 0..s2_len {
            // pick character from BRIGHT_STRING based on the score
            let brightness_char = {
                let score: i64 = sequence_table[[i, j]].score_max(0, 0, 0);

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
                Some((AlignmentChoice::OpenInsert, _, _)) => print!("{}", "I".blue().bold()),
                Some((AlignmentChoice::OpenDelete, _, _)) => print!("{}", "D".cyan().bold()),
                None => print!("{}", brightness_char),
            }
        }
        println!();
    }
}
