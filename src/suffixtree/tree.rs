use std::panic;
use std::time::Instant;

use log::{debug, error, info};

/**
 * Represents a node in the suffix tree
 */
#[derive(Debug, Clone, PartialEq)]
pub struct TreeNode {
    pub id: usize,
    pub string_depth: usize,
    pub string_idx: usize,
    pub edge_start: usize,
    pub edge_end: usize,
    pub parent: Option<usize>,
    pub suffix_link: Option<usize>,
    pub children: Vec<Option<usize>>,
}

pub struct TreeStats {
    pub num_internal: usize,
    pub num_leaves: usize,
    pub num_nodes: usize,
    pub average_string_depth: f64,
    pub max_string_depth: usize,
    pub bwt: String,
    pub longest_repeat_len: usize,
    pub longest_repeat_start: usize,
}

pub struct SuffixTree {
    // nodes are stored in a flat vector
    pub strings: Vec<String>,
    last_internal_id: usize,
    last_leaf_id: usize,
    pub alphabet: Vec<char>,
    pub nodes: Vec<Option<TreeNode>>,
    pub stats: TreeStats,
}

/**
 * Gets the index of a character in the alphabet
 */
pub fn get_child_index(alphabet: &Vec<char>, c: char) -> usize {
    return alphabet
        .iter()
        .enumerate()
        .find(|(_, &x)| x == c)
        .expect(format!("Character {} not found in alphabet", c).as_str())
        .0;
}

impl SuffixTree {
    /**
     * Inserts a new string into the suffix tree
     */
    pub fn insert_string(&mut self, new_string: &str, enable_suffix_links: bool) {
        let string_length: usize = new_string.len();

        // offset all indices for the next string so we have a layout like
        // S1_leaves S1_internal S2_leaves S2_internal
        if self.last_internal_id != 0 {
            self.last_leaf_id = self.last_internal_id;
        }
        self.last_internal_id += string_length + 2;

        // resize nodes array to fit new string
        self.nodes
            .resize(self.nodes.len() + (string_length * 2 + 1), None);

        let before_tree = Instant::now();

        self.strings.push(new_string.to_string() + "$");
        let string_idx = self.strings.len() - 1;

        // build a set of suffixes with '$' appended to the end
        for i in 0..string_length + 1 {
            if enable_suffix_links {
                self.suffix_link_traversal(i, string_idx);
            } else {
                self.find_path(i, 0, string_idx);
            }
        }

        let after_tree = Instant::now();

        let elapsed = after_tree.duration_since(before_tree).as_micros();
        let elapsed_millis = after_tree.duration_since(before_tree).as_millis();
        info!(
            "[FindPath] Time taken to build suffix tree: {} us ({} ms)",
            elapsed, elapsed_millis
        );
    }

    /**
     * Create a new suffix tree with a given number of children
     */
    pub fn new(alphabet_file: &str, initial_allocation: usize) -> SuffixTree {
        // load alphabet from file
        let mut alphabet = match std::fs::read_to_string(alphabet_file) {
            Ok(a) => a.replace(" ", "").chars().collect::<Vec<char>>(),
            Err(_) => {
                panic!("Could not read alphabet file: {}", alphabet_file);
            }
        };

        alphabet.insert(0, '$');

        let mut alphabet_sorted = alphabet.clone();
        alphabet_sorted.sort();

        let mut tree = SuffixTree {
            strings: vec![],
            last_internal_id: 0,
            last_leaf_id: 1,
            alphabet: alphabet_sorted,
            nodes: vec![None; initial_allocation * 2 + 1],
            stats: TreeStats {
                num_internal: 0,
                num_leaves: 0,
                num_nodes: 0,
                average_string_depth: 0.0,
                max_string_depth: 0,
                bwt: "".to_string(),
                longest_repeat_len: 0,
                longest_repeat_start: 0,
            },
        };

        tree.nodes[0] = Some(TreeNode {
            id: 0,
            string_depth: 0,
            parent: None,
            children: vec![None; tree.alphabet.len()],
            edge_start: 0,
            edge_end: 0,
            suffix_link: Some(0),
            string_idx: 0,
        });

        // print the size of a single node in RAM
        info!(
            "Size of a single node in RAM: {} bytes",
            std::mem::size_of_val(&tree.nodes[0])
        );

        return tree;
    }

    /**
     * Follows the suffix link
     */
    fn suffix_link_traversal(&mut self, suffix_idx: usize, string_idx: usize) {
        // update suffix link
        let u_idx = self.nodes[self.last_leaf_id - 1]
            .as_ref()
            .unwrap()
            .parent
            .unwrap_or(0);

        let v_idx = self.nodes[u_idx].as_mut().unwrap().suffix_link;

        let mut should_establish_link = false;

        debug!("u_idx: {}", u_idx);

        // find the next node to go to
        let v = match v_idx {
            Some(v_idx) => {
                // CASE 1
                // suffix link is known (u is not the last node inserted)

                debug!("CASE 1: u = {}", u_idx);
                debug!("CASE 1: v = {}", v_idx);
                if u_idx == 0 {
                    // CASE 1A - u is the root node
                    // go to v (the root)
                    v_idx
                } else {
                    // CASE 1B - u is not the root node
                    // go to v
                    v_idx
                }
            }
            None => {
                // CASE 2
                // suffix link is not known (u is the last node inserted)

                // get u' (the parent of u)
                let u_prime = self.nodes[u_idx].as_ref().unwrap().parent.unwrap();

                // get v' (the suffix link of u')
                let u_prime_ref = self.nodes[u_prime].as_ref().unwrap();
                let v_prime = u_prime_ref
                    .suffix_link
                    .expect(format!("u' {} has no suffix link", u_prime).as_str());

                // get u so we can find the edge leading off of u' (beta')
                let u_ref = self.nodes[u_idx].as_ref().unwrap();

                debug!("CASE 2: u' = {}", u_prime);
                debug!("CASE 2: v' = {}", v_prime);

                let v = match u_prime == 0 {
                    // CASE 2B
                    true => self.node_hops(
                        v_prime,
                        u_ref.edge_start + 1,
                        u_ref.edge_end,
                        u_ref.string_idx,
                    ),

                    // CASE 2A
                    false => {
                        self.node_hops(v_prime, u_ref.edge_start, u_ref.edge_end, u_ref.string_idx)
                    }
                };

                debug!("After NodeHops: v is {}", v);

                should_establish_link = true;

                v
            }
        };

        // establish suffix link from u to v
        if u_idx != 0 && should_establish_link {
            // if the string depth of u is 1, link to the root node
            if self.nodes[u_idx].as_ref().unwrap().string_depth == 1 {
                debug!("Linking u {} to root", u_idx);
                self.nodes[u_idx].as_mut().unwrap().suffix_link = Some(0);
            } else {
                debug!("Linking u {} to v {}", u_idx, v);
                self.nodes[u_idx].as_mut().unwrap().suffix_link = Some(v);
            }
        }

        self.find_path(suffix_idx, v, string_idx);

        // print graphviz at this point
        debug!("\n{}", self.write_graphviz());
    }

    /**
     * Prints the string depth of each node from left to right
     * in the suffix tree
     */
    pub fn display_string_depth(&self, f: &mut std::fmt::Formatter<'_>) {
        writeln!(f, "String Depth: depth (node ID)").unwrap();
        self.dfs(&mut |node: &TreeNode| {
            write!(f, "{} (n{}), ", node.string_depth, node.id).unwrap();
        });
        writeln!(f, "\n").unwrap();
    }

    /**
     * Adds a child node to the given parent node
     */
    pub fn add_child(&mut self, parent: usize, mut child: TreeNode, string_idx: usize) {
        child.parent = Some(parent);

        // figure out where we should insert it
        let child_idx = match self.strings[string_idx].bytes().nth(child.edge_start) {
            Some(c) => get_child_index(&self.alphabet, c as char),
            None => {
                panic!("Child node has no edge label");
            }
        };

        debug!(
            "Adding child onto {} at index {} with first character {}",
            parent,
            child_idx,
            self.strings[string_idx]
                .bytes()
                .nth(child.edge_start)
                .unwrap() as char
        );

        // add to parent's children
        let parent_ref = self.nodes[parent]
            .as_mut()
            .expect("Parent node not found in add_child");
        parent_ref.children[child_idx] = Some(child.id);

        // add to nodes array
        let child_id = child.id;
        self.nodes[child_id] = Some(child);
    }

    /**
     * Performs a Depth First Search (DFS) on the suffix tree
     * executing a callback on each node
     */
    pub fn dfs(&self, callback: &mut dyn FnMut(&TreeNode)) {
        let mut stack = vec![self.nodes[0].as_ref().expect("Root node not found")];

        while let Some(node) = stack.pop() {
            //debug!("DFS: {}", node.id);
            callback(node);
            for child in node.children.iter().rev().flatten() {
                stack.push(self.nodes[*child].as_ref().expect("Child node not found"));
            }
        }
    }

    /**
     * Breaks the edge at the given index and inserts a new node
     * also populates the leaf node with the given label
     */
    pub fn break_edge(
        &mut self,
        node: usize,
        break_idx: usize,
        leaf_start: usize,
        leaf_end: usize,
        create_leaf: bool,
        string_idx: usize,
    ) -> usize {
        // if it doesn't have a parent, panic
        let node_ref = &mut self.nodes[node]
            .as_mut()
            .expect("Node not found in break_edge");

        // Get the current label
        let original_label_start = node_ref.edge_start;
        let original_label_end = node_ref.edge_end;

        // update the current node's edge
        // so that it's the last part of the current label
        node_ref.edge_start = break_idx;
        node_ref.edge_end = original_label_end;

        // create a new internal node
        // with edge = the first part of the current label
        let parent = node_ref
            .parent
            .expect("Node has no parent - cannot break edge");
        let new_internal_node =
            self.create_internal_node(parent, node, original_label_start, break_idx, string_idx);

        if create_leaf {
            // add the leaf on the new internal node
            return self.create_leaf(new_internal_node, leaf_start, leaf_end, string_idx);
        } else {
            return new_internal_node;
        }
    }

    /**
     * Creates a new internal node with the given parent and edge label
     */
    pub fn create_internal_node(
        &mut self,
        parent: usize,
        original_node: usize,
        edge_start: usize,
        edge_end: usize,
        string_idx: usize,
    ) -> usize {
        let parent_ref = self.nodes[parent].as_ref().expect("Parent node not found");

        let internal_id = self.last_internal_id;

        let string_depth = parent_ref.string_depth + (edge_end - edge_start);
        let internal_node = TreeNode {
            id: internal_id,
            string_depth: string_depth,
            parent: Some(parent),
            children: vec![None; self.alphabet.len()],
            edge_start: edge_start,
            edge_end: edge_end,
            // suffix link to last internal node
            suffix_link: None,
            string_idx: string_idx,
        };

        self.last_internal_id += 1;

        // add new internal node as a child of the parent of the original node
        self.add_child(parent, internal_node, string_idx);

        // add original node as a child of the new internal node
        self.add_child(
            internal_id,
            self.nodes[original_node].as_ref().unwrap().clone(),
            string_idx,
        );

        return internal_id;
    }

    /**
     * Creates a new leaf node with the given parent and edge label
     */
    pub fn create_leaf(
        &mut self,
        parent: usize,
        edge_start: usize,
        edge_end: usize,
        string_idx: usize,
    ) -> usize {
        let parent_ref = self.nodes[parent].as_ref().expect("Parent node not found");
        let leaf_id: usize = self.last_leaf_id;

        let leaf = TreeNode {
            id: leaf_id,
            string_depth: parent_ref.string_depth + (edge_end - edge_start),
            parent: Some(parent),
            children: vec![None; self.alphabet.len()],
            edge_start: edge_start,
            edge_end: edge_end,
            suffix_link: None,
            string_idx: string_idx,
        };

        self.last_leaf_id += 1;

        self.add_child(parent, leaf, string_idx);

        return leaf_id;
    }

    /**
     * Hops down the tree until we reach a node where the edge label is
     * greater than our suffix progress
     */
    pub fn node_hops(
        &mut self,
        current_node: usize,
        beta_start: usize,
        beta_end: usize,
        beta_string_idx: usize,
    ) -> usize {
        debug!("[NodeHops] Hop from {}", current_node);

        let beta_length: usize = beta_end
            .checked_sub(beta_start)
            .expect(format!("beta (start): {}, beta (end): {}", beta_start, beta_end).as_str());

        let mut current_node_idx = current_node;

        // if beta is of length 0, return the current node
        if beta_length == 0 {
            debug!(
                "[NodeHops] Can't hop from {} anymore (beta length is 0)",
                current_node_idx
            );
            return current_node_idx;
        }

        let mut remaining_beta = beta_length;
        while remaining_beta > 0 {
            let current_node_ref = self.nodes[current_node_idx].as_ref().unwrap();
            debug!("[NodeHops] Current node: {}", current_node_ref.id);

            let c = self.strings[beta_string_idx].as_bytes()[beta_end - remaining_beta] as char;

            let child_idx: usize = get_child_index(&self.alphabet, c as char);
            let child_node_idx = current_node_ref.children[child_idx];

            match child_node_idx {
                Some(child_node_idx) => {
                    let child_node_ref = self.nodes[child_node_idx].as_ref().unwrap();
                    let child_edge_length = child_node_ref.edge_end - child_node_ref.edge_start;

                    if child_edge_length <= remaining_beta {
                        // hop to the next node
                        remaining_beta -= child_edge_length;
                        current_node_idx = child_node_ref.id;
                        debug!(
                            "[NodeHops] Hopping to child node {} (remaining beta: {})",
                            child_node_ref.id, remaining_beta
                        );
                    } else {
                        debug!(
                            "[NodeHops] Found child node with partial match {}",
                            child_node_ref.id
                        );

                        let mut break_idx = 0;
                        for i in 0..remaining_beta {
                            let c1 = self.strings[child_node_ref.string_idx].as_bytes()
                                [child_node_ref.edge_start + i]
                                as char;
                            let c2 = self.strings[beta_string_idx].as_bytes()
                                [beta_end - remaining_beta + i]
                                as char;

                            if c1 != c2 {
                                debug!("[NodeHops] Found mismatch at {}: {} != {}", i, c1, c2);
                                break;
                            } else {
                                debug!("[NodeHops] Found match at {}: {} == {}", i, c1, c2);
                            }

                            break_idx = child_node_ref.edge_start + i + 1;
                        }

                        // break edge and create internal node (v)
                        current_node_idx = self.break_edge(
                            child_node_idx,
                            break_idx,
                            0,
                            0,
                            false,
                            child_node_ref.string_idx,
                        );
                        break;
                    }
                }
                None => {
                    // return current node
                    debug!(
                        "[NodeHops] Can't hop from {} anymore (no child node)",
                        current_node_ref.id
                    );
                    break;
                }
            }
        }

        debug!("[NodeHops] Returning node {}", current_node_idx);

        return current_node_idx;
    }

    /**
     * Computes the statistics of the suffix tree
     * as well as the Burrows-Wheeler Transform
     */
    pub fn compute_stats(&mut self, string_idx: usize) {
        // compute burrows-wheeler transform
        let mut bwt: Vec<char> = vec![' '; self.strings[string_idx].len()];
        let mut idx: usize = 0;
        let mut num_leaves = 0;
        let mut num_internal = 0;
        let mut max_string_depth = 0;
        let mut string_depth_sum = 0;
        let mut longest_repeat_len = 0;
        let mut longest_repeat_start = 0;
        let mut longest_repeat = false;

        self.dfs(&mut |node: &TreeNode| {
            // if it's a leaf
            if node.id > 0
                && node.id < self.strings[string_idx].len() + 1
                && idx < self.strings[string_idx].len()
            {
                // if we previously found a longest repeat (internal node)
                // set the longest repeat start to the node id
                if longest_repeat {
                    longest_repeat_start = node.id;
                    longest_repeat = false;
                }

                num_leaves += 1;
                if node.id == 1 {
                    bwt[idx] = '$';
                } else {
                    bwt[idx] = self.strings[string_idx].as_bytes()[node.id - 2] as char;
                }
                idx += 1;
            } else {
                // don't count the root node
                if node.id == 0 {
                    return;
                }

                num_internal += 1;
                string_depth_sum += node.string_depth;
                if node.string_depth > max_string_depth {
                    // this is also the longest matching repeat
                    longest_repeat_len = node.string_depth;
                    longest_repeat = true;
                    max_string_depth = node.string_depth;
                }
            }
        });

        self.stats.longest_repeat_len = longest_repeat_len;
        self.stats.longest_repeat_start = longest_repeat_start;

        self.stats.num_leaves = num_leaves;
        self.stats.num_internal = num_internal;
        self.stats.num_nodes = self.stats.num_internal + self.stats.num_leaves + 1;
        self.stats.bwt = bwt.iter().collect::<String>().trim().to_string();

        self.stats.average_string_depth = string_depth_sum as f64 / self.stats.num_internal as f64;
        self.stats.max_string_depth = max_string_depth;
    }

    /**
     * Walks down the tree and inserts new leaf for the given suffix
     */
    pub fn find_path(&mut self, suffix_idx: usize, start_node: usize, string_idx: usize) {
        let mut current_node = self.nodes[start_node]
            .as_ref()
            .expect(format!("Start node {} not found", start_node).as_str());
        let suffix_len = self.strings[string_idx].len() - suffix_idx;

        debug!("Current node: {}", current_node.id);

        // how far we've already walked down the suffix
        let mut suffix_sub_idx =
            current_node.string_depth - (current_node.edge_end - current_node.edge_start);

        debug!(
            "Suffix sub-index: {}, suffix index: {}",
            suffix_sub_idx, suffix_idx
        );

        loop {
            // walk down label on current node's edge
            for label_idx in current_node.edge_start..current_node.edge_end {
                if suffix_sub_idx > suffix_len {
                    debug!("Reached end of suffix");
                    break;
                }

                let suffix_char = self.strings[string_idx].as_bytes()[suffix_idx + suffix_sub_idx];
                let c = self.strings[current_node.string_idx].as_bytes()[label_idx];

                // if the suffix character is not equal to the edge character
                // break a new edge and insert a new leaf
                debug!(
                    "Comparing [{}] {} to {} [{}]",
                    suffix_sub_idx, suffix_char as char, c as char, label_idx
                );

                if suffix_char != c {
                    self.break_edge(
                        current_node.id,
                        label_idx,
                        // leaf label is the position we're at in the suffix
                        // and the index of the end of the original string (since it's a suffix and we're at the end)
                        suffix_idx + suffix_sub_idx,
                        self.strings[string_idx].len(),
                        true,
                        string_idx,
                    );
                    return;
                }

                suffix_sub_idx += 1;
            }

            debug!("Done with edge of current node {}", current_node.id);

            if suffix_sub_idx >= suffix_len {
                error!(
                    "ERR: Suffix sub-idx {} is greater than suffix length {}",
                    suffix_idx, suffix_len
                );
                return;
            }

            // if we've reached the end of the edge label, move to the next node
            // compare the first character of the edge label of each child node
            // to the next character in the suffix
            let c = self.strings[string_idx].as_bytes()[suffix_idx + suffix_sub_idx] as char;

            debug!(
                "Checking for child with {}, suffix index: {}",
                c, suffix_sub_idx
            );

            let child_idx: usize = get_child_index(&self.alphabet, c);
            let child_node = current_node.children[child_idx];

            match child_node {
                // there's already a child node, so we need to keep moving down the tree
                Some(n) => {
                    debug!("Found next child node for {}, child idx: {}", c, child_idx);
                    current_node = self.nodes[n]
                        .as_ref()
                        .expect("Child node not found when walking");
                }
                // there's no child node, so we need to add a new leaf
                None => {
                    debug!("No child node found, adding new leaf");
                    self.create_leaf(
                        current_node.id,
                        suffix_idx + suffix_sub_idx,
                        self.strings[string_idx].len(),
                        string_idx,
                    );
                    return;
                }
            }
        }
    }
}

#[cfg(test)]
mod test {

    use std::env;

    use crate::sequence::{SequenceContainer, SequenceOperations};

    use super::SuffixTree;

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
        tree.compute_stats(0);

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
