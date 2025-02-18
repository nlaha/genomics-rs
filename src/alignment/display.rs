use colored::Colorize;
use log::{info, warn};
use std::fmt::Display;

use super::algo::{AlignedSequences, AlignmentChoice};

const DISP_MAX_WIDTH: usize = 200;

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
                    if s1_idx < self.s1.sequence.len() {
                        s1_out.push(self.s1.sequence.as_bytes()[s1_idx] as char);
                        s1_idx += 1;
                    }
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
                    if s2_idx < self.s2.sequence.len() {
                        s2_out.push(self.s2.sequence.as_bytes()[s2_idx] as char);
                        s2_idx += 1;
                    }
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
        )?;
        writeln!(
            f,
            "Percent Identity {}%",
            (self.matches as f64 / align_idx as f64) * 100.0
        )
    }
}

/// Prints the sequence table with alignment path for visualization
/// * `aligned_sequences` - the aligned sequences
pub fn print_sequence_table(aligned_sequences: &AlignedSequences) {
    let s1_len = aligned_sequences.s1.sequence.len();
    let s2_len = aligned_sequences.s2.sequence.len();

    // skip visualization if the table is too large
    if s1_len < DISP_MAX_WIDTH && s2_len < DISP_MAX_WIDTH * 10 {
        info!("Computing sequence table visualization...");
    } else {
        warn!("Sequence table too large to visualize");
        return;
    }

    println!("\nSequence Table (S1 columns, S2 rows):\n");

    // print the columns (bases of the second sequence)
    print!(" ",);
    for i in 0..s2_len {
        print!("{}", aligned_sequences.s2.sequence.chars().nth(i).unwrap());
    }
    println!();

    // display dynamic programming table
    for i in 0..s1_len {
        // print sequence base
        print!("{}", aligned_sequences.s1.sequence.chars().nth(i).unwrap());
        for j in 0..s2_len {
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
                None => print!("."),
            }
        }
        println!();
    }
}
