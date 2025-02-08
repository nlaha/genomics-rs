use genomics_rs::config;
use genomics_rs::config::Config;

const TEST_CONFIG: Config = Config {
    scores: config::Scores {
        s_match: 1,
        s_mismatch: -2,
        g: -2,
        h: -5,
    },
    fasta: config::FastaConfig {
        path: String::new(),
    },
};

#[cfg(test)]
mod test_alignment {
    use genomics_rs::{
        alignment::{self, AlignmentChoice},
        sequence::{Sequence, SequenceContainer},
    };

    use crate::TEST_CONFIG;

    #[test]
    fn test_simple_matches() {
        let s1: Sequence = Sequence {
            name: "s1".to_string(),
            sequence: "ACGT".to_string(),
        };

        let s2: Sequence = Sequence {
            name: "s2".to_string(),
            sequence: "ACGT".to_string(),
        };

        let sc: SequenceContainer = SequenceContainer {
            sequences: vec![s1, s2],
        };

        let aligned_sequences = alignment::global_alignment(sc, TEST_CONFIG.scores);

        assert_eq!(aligned_sequences.score, 4);
        assert_eq!(aligned_sequences.matches, 4);
        assert_eq!(aligned_sequences.mismatches, 0);
        assert_eq!(aligned_sequences.opening_gaps, 0);
        assert_eq!(aligned_sequences.gap_extensions, 0);
        assert_eq!(
            aligned_sequences.alignment,
            vec![
                (AlignmentChoice::Match, 3, 3),
                (AlignmentChoice::Match, 2, 2),
                (AlignmentChoice::Match, 1, 1),
                (AlignmentChoice::Match, 0, 0),
            ]
        );
    }

    #[test]
    fn test_gaps() {
        let s1: Sequence = Sequence {
            name: "s1".to_string(),
            sequence: "ACGT".to_string(),
        };

        let s2: Sequence = Sequence {
            name: "s2".to_string(),
            sequence: "AACGT".to_string(),
        };

        let sc: SequenceContainer = SequenceContainer {
            sequences: vec![s1, s2],
        };

        let aligned_sequences = alignment::global_alignment(sc, TEST_CONFIG.scores);

        assert_eq!(aligned_sequences.score, -1);
        assert_eq!(aligned_sequences.matches, 4);
        assert_eq!(aligned_sequences.mismatches, 0);
        assert_eq!(aligned_sequences.opening_gaps, 1);
        assert_eq!(aligned_sequences.gap_extensions, 0);
        assert_eq!(
            aligned_sequences.alignment,
            vec![
                (AlignmentChoice::Match, 3, 4),
                (AlignmentChoice::Match, 2, 3),
                (AlignmentChoice::Match, 1, 2),
                (AlignmentChoice::OpenInsert, 0, 1),
                (AlignmentChoice::Match, 0, 0)
            ]
        );
    }

    #[test]
    fn test_affine_gap() {
        let s1: Sequence = Sequence {
            name: "s1".to_string(),
            sequence: "ACGGATAAAAAAAATC".to_string(),
        };

        let s2: Sequence = Sequence {
            name: "s2".to_string(),
            sequence: "ACGGATAAAATC".to_string(),
        };

        let sc: SequenceContainer = SequenceContainer {
            sequences: vec![s1, s2],
        };

        let aligned_sequences = alignment::global_alignment(sc, TEST_CONFIG.scores);

        assert_eq!(aligned_sequences.score, -1);
        assert_eq!(aligned_sequences.matches, 10);
        assert_eq!(aligned_sequences.mismatches, 0);
        assert_eq!(aligned_sequences.opening_gaps, 1);
        assert_eq!(aligned_sequences.gap_extensions, 3);
    }
}
