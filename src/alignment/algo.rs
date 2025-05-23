use std::{
    ops::Add,
    simd::{
        i64x4,
        num::{SimdInt, SimdUint},
        usizex4,
    },
};

use log::{debug, error, info, warn};
use ndarray::{Array2, ShapeBuilder};
use num::Zero;
use spinoff::{spinners, Color, Spinner};

use crate::{
    alignment::display::print_alignment_table,
    config::Scores,
    sequence::{Sequence, SequenceContainer, SequenceOperations},
};

/// Cell in the dynamic programming table
/// * `insert_score` - score for an insertion
/// * `delete_score` - score for a deletion
/// * `sub_score` - score for a substitution
#[derive(PartialEq, Debug, Clone, Copy, Default)]
#[repr(C)]
pub struct AlignmentCell {
    pub insert_score: i64,
    pub delete_score: i64,
    pub sub_score: i64,

    pub insert_matches: usize,
    pub delete_matches: usize,
    pub sub_matches: usize,
}

impl Add for AlignmentCell {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        AlignmentCell {
            insert_score: self.insert_score + other.insert_score,
            delete_score: self.delete_score + other.delete_score,
            sub_score: self.sub_score + other.sub_score,

            insert_matches: self.insert_matches + other.insert_matches,
            delete_matches: self.delete_matches + other.delete_matches,
            sub_matches: self.sub_matches + other.sub_matches,
        }
    }
}

impl Zero for AlignmentCell {
    fn zero() -> Self {
        AlignmentCell {
            insert_score: 0,
            delete_score: 0,
            sub_score: 0,

            insert_matches: 0,
            delete_matches: 0,
            sub_matches: 0,
        }
    }

    fn is_zero(&self) -> bool {
        self.insert_score == 0
            && self.delete_score == 0
            && self.sub_score == 0
            && self.insert_matches == 0
            && self.delete_matches == 0
            && self.sub_matches == 0
    }

    fn set_zero(&mut self) {
        self.insert_score = 0;
        self.delete_score = 0;
        self.sub_score = 0;

        self.insert_matches = 0;
        self.delete_matches = 0;
        self.sub_matches = 0;
    }
}

/// Handles score calculation for alignment cells
pub trait ComputeScore {
    fn score_max(&self, i_mod: i64, s_mod: i64, d_mod: i64, is_local: bool) -> i64;
    fn max_matches(&self) -> usize;
}

impl ComputeScore for AlignmentCell {
    /// Computes the maximum score for the cell
    /// * `i_mod` - insertion score modifier
    /// * `s_mod` - substitution score modifier
    /// * `d_mod` - deletion score modifier
    /// * Returns the maximum score (where the modifiers are applied to each)
    fn score_max(&self, i_mod: i64, s_mod: i64, d_mod: i64, is_local: bool) -> i64 {
        let v: std::simd::Simd<i64, 4> = i64x4::from_array([
            self.insert_score + i_mod,
            self.sub_score + s_mod,
            self.delete_score + d_mod,
            if is_local { 0 } else { i64::MIN },
        ]);

        v.reduce_max()
    }

    /**
     * Finds the maximum number of matches in the cell from each direction
     */
    fn max_matches(&self) -> usize {
        let v: std::simd::Simd<usize, 4> = usizex4::from_array([
            self.insert_matches,
            self.sub_matches,
            self.delete_matches,
            0,
        ]);

        v.reduce_max() as usize
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
/// * `sequence_container` - container struct containing two sequences
/// * `scores` - scoring parameters
pub fn alignment_table(
    sequence_container: &SequenceContainer,
    scores: &Scores,
    is_local: bool,
    reverse_sequences: bool,
) -> (Array2<AlignmentCell>, usize) {
    let mut maximum_score = i64::MIN;
    let mut max_cell = (0, 0);

    // if more than two sequences are found
    if sequence_container.sequences.len() > 2 {
        warn!("More than two sequences found. Only the first two will be used.");
    }

    // we need to offset our 'negative infinity' to avoid overflow errors
    let negative_inf: i64 = i64::MIN + num::abs(scores.g + scores.h);

    let s1_len = sequence_container.sequences[0].sequence.len();
    let s2_len = sequence_container.sequences[1].sequence.len();

    // create 2D array to store dynamic programming results
    let mut alignment_table: Array2<AlignmentCell> = Array2::zeros((s1_len + 1, s2_len + 1).f());

    // log
    info!("Sequence table shape: {:?}", alignment_table.shape());
    info!(
        "Sequence table size (KB): {}",
        alignment_table.len() * std::mem::size_of::<AlignmentCell>() / 1024
    );

    let mut alignment_table_spinner = Spinner::new(
        spinners::Dots,
        "Computing sequence table...",
        Color::Magenta,
    );

    // start timer
    let start_table_init: std::time::Instant = std::time::Instant::now();

    // iterate through alignment table and fill in scores
    for i in 0..s1_len + 1 {
        for j in 0..s2_len + 1 {
            alignment_table[[i, j]] = match (i, j) {
                // origin
                (0, 0) => AlignmentCell {
                    insert_score: 0,
                    delete_score: 0,
                    sub_score: 0,
                    insert_matches: 0,
                    delete_matches: 0,
                    sub_matches: 0,
                },
                // column
                (i, 0) => AlignmentCell {
                    insert_score: negative_inf,
                    delete_score: scores.h + (i as i64 * scores.g),
                    sub_score: negative_inf,
                    insert_matches: 0,
                    delete_matches: 0,
                    sub_matches: 0,
                },
                // row
                (0, j) => AlignmentCell {
                    insert_score: scores.h + (j as i64 * scores.g),
                    delete_score: negative_inf,
                    sub_score: negative_inf,
                    insert_matches: 0,
                    delete_matches: 0,
                    sub_matches: 0,
                },
                (i, j) => {
                    // recurrence relationships for S, D and I
                    let top_left = alignment_table[[i - 1, j - 1]];
                    let left = alignment_table[[i - 1, j]];
                    let top = alignment_table[[i, j - 1]];

                    let is_match = sequence_container.is_match(i - 1, j - 1, reverse_sequences);

                    // build alignment cell struct
                    let cell = AlignmentCell {
                        insert_score: top.score_max(
                            scores.g,
                            scores.h + scores.g,
                            scores.h + scores.g,
                            is_local,
                        ),

                        delete_score: left.score_max(
                            scores.h + scores.g,
                            scores.h + scores.g,
                            scores.g,
                            is_local,
                        ),

                        sub_score: match is_match {
                            true => scores.s_match + top_left.score_max(0, 0, 0, is_local),
                            false => scores.s_mismatch + top_left.score_max(0, 0, 0, is_local),
                        },

                        insert_matches: top.max_matches(),
                        delete_matches: left.max_matches(),
                        sub_matches: match is_match {
                            true => top_left.max_matches() + 1,
                            false => top_left.max_matches(),
                        },
                    };

                    let max_cell_score = cell.score_max(0, 0, 0, is_local);
                    if maximum_score < max_cell_score {
                        max_cell = (i, j);
                        maximum_score = max_cell_score;
                    }

                    cell
                }
            };
        }
    }

    let end_table_init: std::time::Instant = std::time::Instant::now();

    alignment_table_spinner.success("Sequence table computed");

    info!(
        "Table initialization complete, time taken: {}us",
        (end_table_init - start_table_init).as_micros()
    );

    let matches_at_max = alignment_table[max_cell].max_matches();

    return (alignment_table, matches_at_max);
}

/// Performs a retrace of the optimal alignment path from the sequence table
/// * `sequence_container` - container struct containing two sequences
/// * `alignment_table` - the dynamic programming table
pub fn retrace(
    sequence_container: &SequenceContainer,
    alignment_table: ndarray::Array2<AlignmentCell>,
    is_local: bool,
) -> AlignedSequences {
    let mut retrace_spinner = Spinner::new(
        spinners::Dots,
        "Retracing optimal alignment...",
        Color::Magenta,
    );

    let s1_len = sequence_container.sequences[0].sequence.len();
    let s2_len: usize = sequence_container.sequences[1].sequence.len();

    // optimal retrace algorithm
    let start_retrace: std::time::Instant = std::time::Instant::now();

    // start at end if performing global alignment
    // or at the highest scoring cell if performing local alignment
    let (mut i, mut j) = match is_local {
        // our sequence table is 1-indexed
        false => (s1_len, s2_len),
        // get coordinates of the highest scoring cell
        true => {
            alignment_table
                // get indexed iter that gives us (x,y), cell
                .indexed_iter()
                // max by lambda that compares the score of the cell
                .max_by(|a, b| {
                    a.1.score_max(0, 0, 0, is_local)
                        .cmp(&b.1.score_max(0, 0, 0, is_local))
                })
                .unwrap()
                // get the index of the max cell
                .0
        }
    };

    info!("Starting at ({}, {})", i, j);

    let mut aligned_sequences: AlignedSequences = AlignedSequences {
        s1: sequence_container.sequences[0].clone(),
        s2: sequence_container.sequences[1].clone(),
        alignment: Vec::new(),
        score: alignment_table[[i, j]].score_max(0, 0, 0, is_local),
        matches: 0,
        mismatches: 0,
        gap_extensions: 0,
        opening_gaps: 0,
    };

    let mut last_choice = AlignmentChoice::Match;
    loop {
        let cell = alignment_table[[i, j]];

        // get the score of the current cell
        let max = cell.score_max(0, 0, 0, is_local);

        // debug log cell
        debug!(
            "({}, {}) -> max:{}, i:{}, s:{}, d:{}",
            i, j, max, cell.insert_score, cell.sub_score, cell.delete_score
        );

        let (i_opt, j_opt) = match max {
            // top_left move
            val if val == cell.sub_score => {
                if sequence_container.is_match(i, j, false) {
                    last_choice = AlignmentChoice::Match;
                    aligned_sequences.matches += 1;
                    aligned_sequences
                        .alignment
                        .push((AlignmentChoice::Match, i, j));
                    debug!("Match found at ({}, {}) -> {}", i, j, val);
                } else {
                    last_choice = AlignmentChoice::Mismatch;
                    aligned_sequences.mismatches += 1;
                    aligned_sequences
                        .alignment
                        .push((AlignmentChoice::Mismatch, i, j));
                    debug!("Mismatch found at ({}, {}) -> {}", i, j, val);
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
                debug!("Insert found at ({}, {}) -> {}", i, j, val);
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
                debug!("Delete found at ({}, {}) -> {}", i, j, val);
                last_choice = AlignmentChoice::Delete;
                (i.checked_sub(1), Some(j))
            }
            _ => {
                if is_local && max == 0 {
                    info!("Ending local alignment at ({}, {})", i, j);
                    break;
                }

                error!("No valid move found at ({}, {}) -> {}", i, j, max);
                panic!("Unexpected score during retrace: {}", max);
            }
        };

        (i, j) = match (i_opt, j_opt) {
            (Some(i), Some(j)) => (i, j),
            (None, None) => break,
            (Some(i), None) => (i, 0),
            (None, Some(j)) => (0, j),
        };

        if i == 0 && j == 0 {
            break;
        }
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

    print_alignment_table(&aligned_sequences, &alignment_table);

    aligned_sequences
}
