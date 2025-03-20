use genomics_rs::config;
use genomics_rs::config::Config;

const TEST_CONFIG: Config = Config {
    scores: config::Scores {
        s_match: 1,
        s_mismatch: -2,
        g: -2,
        h: -5,
    },
};

#[cfg(test)]
mod test_alignment {
    use genomics_rs::alignment::algo::{AlignedSequences, AlignmentChoice::*};
    use genomics_rs::{
        alignment::{self},
        sequence::{Sequence, SequenceContainer},
    }; // Import enum variants

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

        let aligned_sequences: AlignedSequences =
            alignment::algo::align_sequences(&sc, &TEST_CONFIG.scores, false);

        println!("{}", aligned_sequences);
        assert_eq!(aligned_sequences.score, 4);
        assert_eq!(aligned_sequences.matches, 4);
        assert_eq!(aligned_sequences.mismatches, 0);
        assert_eq!(aligned_sequences.opening_gaps, 0);
        assert_eq!(aligned_sequences.gap_extensions, 0);
        assert_eq!(
            aligned_sequences.alignment,
            vec![(Match, 4, 4), (Match, 3, 3), (Match, 2, 2), (Match, 1, 1)]
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
            sequence: "AGCGT".to_string(),
        };

        let sc: SequenceContainer = SequenceContainer {
            sequences: vec![s1, s2],
        };

        let aligned_sequences: AlignedSequences =
            alignment::algo::align_sequences(&sc, &TEST_CONFIG.scores, false);

        println!("{}", aligned_sequences);
        assert_eq!(aligned_sequences.matches, 3);
        assert_eq!(aligned_sequences.mismatches, 1);
        assert_eq!(aligned_sequences.opening_gaps, 1);
        assert_eq!(aligned_sequences.gap_extensions, 0);
        assert_eq!(
            aligned_sequences.alignment,
            vec![
                (Match, 4, 5),
                (Match, 3, 4),
                (Match, 2, 3),
                (OpenInsert, 1, 2),
                (Mismatch, 1, 1)
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

        let aligned_sequences: AlignedSequences =
            alignment::algo::align_sequences(&sc, &TEST_CONFIG.scores, false);

        println!("{}", aligned_sequences);
        assert_eq!(aligned_sequences.matches, 12);
        assert_eq!(aligned_sequences.mismatches, 0);
        assert_eq!(aligned_sequences.opening_gaps, 1);
        assert_eq!(aligned_sequences.gap_extensions, 3);

        assert_eq!(
            aligned_sequences.alignment,
            vec![
                (Match, 16, 12),
                (Match, 15, 11),
                (Match, 14, 10),
                (Match, 13, 9),
                (Match, 12, 8),
                (Match, 11, 7),
                (OpenDelete, 10, 6),
                (Delete, 9, 6),
                (Delete, 8, 6),
                (Delete, 7, 6),
                (Match, 6, 6),
                (Match, 5, 5),
                (Match, 4, 4),
                (Match, 3, 3),
                (Match, 2, 2),
                (Match, 1, 1)
            ]
        );
    }
}
