use log::{debug, error, info, trace, warn};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

/// A struct to hold a FASTA sequence and its name.
#[derive(Debug, Clone)]
pub struct FastaSequence {
    pub name: String,
    pub sequence: String,
}

/// Load one or more sequences from a FASTA file.
/// @param filepath: Path to the FASTA file
pub fn load_sequences(filepath: &str) -> Vec<FastaSequence> {
    let mut sequences: Vec<FastaSequence> = Vec::new();
    let mut sequence_name: Option<String> = None;

    // read file line by line
    // look for ">" to identify sequence headers
    if let Ok(lines) = read_lines(filepath) {
        for line in lines.map_while(Result::ok) {
            // skip if empty
            if line.is_empty() {
                continue;
            }

            if line.starts_with(">") {
                // get the name minus the '>'
                let name: String = line[1..].trim().to_string();
                info!("Sequence Found: {}", name);
                sequences.push(FastaSequence {
                    name: name.clone(),
                    sequence: String::new(),
                });

                sequence_name = Some(name);
            } else {
                if let Some(sequence_name) = sequence_name.as_mut() {
                    // trim any whitespace in the line and append to the current sequence
                    let sequence = sequences.last_mut().unwrap();
                    sequence.sequence.push_str(&line.trim());
                } else {
                    warn!("Sequence data found without a header");
                }
            }
        }
    } else {
        error!("Could not open file: {}", filepath);
    }

    // print some debug info
    debug!("Loaded {} sequences", sequences.len());
    for (name, sequence) in sequences.iter().map(|s| (&s.name, &s.sequence)) {
        debug!("  {}: {} bases", name, sequence.len());
    }

    return sequences;
}

// Source: https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html
// The output is wrapped in a Result to allow matching on errors.
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
