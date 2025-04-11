#[cfg(test)]
mod test_suffixtree {

    use genomics_rs::sequence::{SequenceContainer, SequenceOperations};
    use genomics_rs::suffixtree::tree::SuffixTree;

    #[test]
    fn test_tree_simple2() {
        let mut tree = SuffixTree::new("alphabets/dna.txt", 10);
        tree.insert_string("ACA", true);
        tree.compute_stats(0);

        assert_eq!(tree.stats.num_nodes, 6);
    }

    #[test]
    fn test_tree_simple3() {
        let mut tree = SuffixTree::new("alphabets/banana.txt", 10);
        tree.insert_string("BANANA", true);
        tree.compute_stats(0);

        println!("{}", tree);

        assert_eq!(tree.stats.num_internal, 3);
        assert_eq!(tree.stats.num_leaves, 7);
        assert_eq!(tree.stats.num_nodes, 11);

        assert_eq!(tree.stats.average_string_depth, 2.0);
        assert_eq!(tree.stats.max_string_depth, 3);
        assert_eq!(tree.stats.bwt, "ANNB$AA".to_string());
    }

    #[test]
    fn test_tree_simple4() {
        let mut tree = SuffixTree::new("alphabets/english.txt", 20);
        tree.insert_string("MISSISSIPPI", true);
        tree.compute_stats(0);

        println!("{}", tree);

        assert_eq!(tree.stats.num_internal, 6);
        assert_eq!(tree.stats.num_leaves, 12);
        assert_eq!(tree.stats.num_nodes, 19);

        assert_eq!(tree.stats.average_string_depth, 2.0);
        assert_eq!(tree.stats.max_string_depth, 4);
        assert_eq!(tree.stats.bwt, "IPSSM$PISSII".to_string());
    }

    #[test]
    fn test_tree_covid_wuhan() {
        let mut sequence_container: SequenceContainer = SequenceContainer {
            sequences: Vec::new(),
        };

        sequence_container.from_fasta("test_data/Covid_Wuhan.fasta");

        let mut suffix_tree = SuffixTree::new("alphabets/dna.txt", 200000);
        suffix_tree.insert_string(&sequence_container.sequences[0].sequence, true);
        suffix_tree.compute_stats(0);

        assert_eq!(suffix_tree.stats.num_internal, 19098);
        assert_eq!(suffix_tree.stats.num_leaves, 29904);
        assert_eq!(suffix_tree.stats.num_nodes, 49003);

        // load BWT from file and compare to the computed BWT line by line
        let bwt = std::fs::read_to_string("BWTs/Covid_Wuhan.fasta.BWT.out")
            .unwrap()
            .replace("\n", "");

        let mut idx = 0;
        for (computed, expected) in suffix_tree.stats.bwt.chars().zip(bwt.chars()) {
            assert_eq!(
                computed, expected,
                "[idx {}] Computed: {}, Expected: {}",
                idx, computed, expected
            );
            idx += 1;
        }
    }

    #[test]
    fn test_tree_human_brca2() {
        let mut sequence_container: SequenceContainer = SequenceContainer {
            sequences: Vec::new(),
        };

        sequence_container.from_fasta("test_data/Human-BRCA2-cds.fasta");

        let mut suffix_tree = SuffixTree::new("alphabets/dna.txt", 200000);
        suffix_tree.insert_string(&sequence_container.sequences[0].sequence, true);
        suffix_tree.compute_stats(0);

        assert_eq!(suffix_tree.stats.num_internal, 7299);
        assert_eq!(suffix_tree.stats.num_leaves, 11383);
        assert_eq!(suffix_tree.stats.num_nodes, 18683);

        // load BWT from file and compare to the computed BWT line by line
        let bwt = std::fs::read_to_string("BWTs/Human-BRCA2-cds.fasta.BWT.txt")
            .unwrap()
            .replace("\n", "");

        for (computed, expected) in suffix_tree.stats.bwt.chars().zip(bwt.chars()) {
            assert_eq!(computed, expected);
        }
    }

    #[test]
    fn test_tree_slyco() {
        let mut sequence_container: SequenceContainer = SequenceContainer {
            sequences: Vec::new(),
        };

        sequence_container.from_fasta("test_data/Slyco.fasta");

        let mut suffix_tree = SuffixTree::new("alphabets/dna.txt", 200000);
        suffix_tree.insert_string(&sequence_container.sequences[0].sequence, true);
        suffix_tree.compute_stats(0);

        assert_eq!(suffix_tree.stats.num_internal, 98972);
        assert_eq!(suffix_tree.stats.num_leaves, 155462);
        assert_eq!(suffix_tree.stats.num_nodes, 254435);

        // load BWT from file and compare to the computed BWT line by line
        let bwt = std::fs::read_to_string("BWTs/Slyco.fas.BWT.out")
            .unwrap()
            .replace("\n", "");

        for (computed, expected) in suffix_tree.stats.bwt.chars().zip(bwt.chars()) {
            assert_eq!(computed, expected);
        }
    }

    #[test]
    fn test_generalized_suffix_tree() {
        let mut tree = SuffixTree::new("alphabets/banana.txt", 10);
        tree.insert_string("BANANA", true);
        tree.insert_string("ABANANA", true);
        tree.compute_stats(0);

        let (s1, s2, length) = tree.get_lcs(0, 1);
        println!("LCS: S1 -> {}, S2 -> {}, Length -> {}", s1, s2, length);

        assert_eq!(s1, 0);
        assert_eq!(s2, 1);
        assert_eq!(length, 6);

        println!("{}", tree);
    }

    #[test]
    fn test_generalized_suffix_tree2() {
        // // set default log level
        // env::set_var("RUST_LOG", "debug");

        // // init logging
        // pretty_env_logger::init();

        let mut tree = SuffixTree::new("alphabets/banana.txt", 10);
        tree.insert_string("BANANA", true);
        tree.insert_string("BANANAB", true);
        tree.insert_string("ABABABA", true);
        tree.compute_stats(0);

        let (s1, s2, length) = tree.get_lcs(1, 2);
        println!("LCS: S1 -> {}, S2 -> {}, Length -> {}", s1, s2, length);

        assert_eq!(s1, 5);
        assert_eq!(s2, 4);
        assert_eq!(length, 2);

        println!("{}", tree);
    }

    // #[test]
    // fn test_tree_chr12() {
    //     let mut sequence_container: SequenceContainer = SequenceContainer {
    //         sequences: Vec::new(),
    //     };

    //     sequence_container.from_fasta("test_data/chr12.fasta");

    //     let mut suffix_tree = SuffixTree::new(
    //         &sequence_container.sequences[0].sequence,
    //         "alphabets/dna.txt",
    //         true,
    //     );
    //     suffix_tree.compute_stats(false);

    //     assert_eq!(suffix_tree.stats.num_internal, 699519);
    //     assert_eq!(suffix_tree.stats.num_leaves, 1078176);
    //     assert_eq!(suffix_tree.stats.num_nodes, 1777696);
    // }
}
