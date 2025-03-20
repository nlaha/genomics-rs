use std::cell::RefCell;
use std::panic;
use std::rc::Rc;

use log::{debug, info};

/**
 * Represents a node in the suffix tree
 */
#[derive(Debug, Clone, PartialEq)]
pub struct TreeNode {
    pub id: usize,
    pub string_depth: usize,
    pub parent: Option<usize>,
    pub children: Vec<Option<usize>>,
    pub edge_start: usize,
    pub edge_end: usize,
}

pub struct TreeStats {
    pub num_internal: usize,
    pub num_leaves: usize,
    pub num_nodes: usize,
    pub average_string_depth: f64,
    pub max_string_depth: usize,
    pub bwt: String,
}

pub struct SuffixTree {
    // nodes are stored in a flat vector
    pub original_string: String,
    last_internal_id: usize,
    last_leaf_id: usize,
    pub suffixes: Vec<String>,
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
     * Create a new suffix tree with a given number of children
     */
    pub fn new(original_string: &str, alphabet_file: &str) -> SuffixTree {
        let string_length = original_string.len();

        // load alphabet from file
        let alphabet = match std::fs::read_to_string(alphabet_file) {
            Ok(a) => a.chars().collect::<Vec<char>>(),
            Err(_) => {
                panic!("Could not read alphabet file: {}", alphabet_file);
            }
        };

        let mut alphabet_sorted = alphabet.clone();
        alphabet_sorted.sort();

        let mut tree = SuffixTree {
            original_string: original_string.to_string() + "$",
            last_internal_id: string_length + 2,
            last_leaf_id: 1,
            suffixes: Vec::with_capacity(string_length),
            alphabet: alphabet_sorted,
            nodes: vec![None; string_length * 2 + 2],
            stats: TreeStats {
                num_internal: 0,
                num_leaves: 0,
                num_nodes: 0,
                average_string_depth: 0.0,
                max_string_depth: 0,
                bwt: "".to_string(),
            },
        };

        tree.nodes[0] = Some(TreeNode {
            id: 0,
            string_depth: 0,
            parent: None,
            children: vec![None; tree.alphabet.len()],
            edge_start: 0,
            edge_end: 0,
        });

        // build a set of suffixes with '$' appended to the end
        for i in 0..string_length + 1 {
            let suffix = tree.original_string[i..].to_string();
            // if suffix is longer than 100 characters, truncate it
            if suffix.len() > 100 {
                info!("[FindPath] {}/{} {}...", i, string_length, &suffix[..100]);
            } else {
                info!("[FindPath] {}/{} {}", i, string_length, suffix);
            }
            tree.suffixes.push(suffix.to_string());
            tree.find_path(i);
        }

        tree.stats.average_string_depth = tree
            .nodes
            .iter()
            .flatten()
            .map(|n| n.string_depth)
            .sum::<usize>() as f64
            / tree.nodes.len() as f64;
        tree.stats.max_string_depth = tree
            .nodes
            .iter()
            .flatten()
            .map(|n| n.string_depth)
            .max()
            .unwrap_or(0);

        // compute burrows-wheeler transform
        let mut bwt: Vec<char> = vec![' '; tree.original_string.len()];
        let mut idx = 0;
        let mut num_leaves = 0;
        let mut num_internal = 0;

        tree.dfs(&mut |node: TreeNode| {
            // if it's a leaf
            if node.id > 0 && node.id < tree.suffixes.len() + 1 && idx < tree.original_string.len()
            {
                num_leaves += 1;
                if node.id == 1 {
                    bwt[idx] = '$';
                } else {
                    bwt[idx] = tree.original_string.as_bytes()[node.id - 2] as char;
                }
                idx += 1;
            } else {
                num_internal += 1;
            }
        });

        tree.stats.num_leaves = num_leaves;
        tree.stats.num_internal = num_internal;
        tree.stats.num_nodes = tree.stats.num_internal + tree.stats.num_leaves;
        tree.stats.bwt = bwt.iter().collect::<String>().trim().to_string();

        return tree;
    }

    pub fn add_child(&mut self, parent: usize, child: TreeNode) {
        // figure out where we should insert it
        let child_idx = match self.original_string.bytes().nth(child.edge_start) {
            Some(c) => get_child_index(&self.alphabet, c as char),
            None => {
                panic!("Child node has no edge label");
            }
        };

        // add to parent's children
        let parent_ref = self.nodes[parent]
            .as_mut()
            .expect("Parent node not found in add_child");
        parent_ref.children[child_idx] = Some(child.id);

        // add to nodes array
        self.nodes[child.id] = Some(child.clone());

        // update the child's parent
        let child_ref = self.nodes[child.id].as_mut().unwrap();
        child_ref.parent = Some(parent);
    }

    /**
     * Performs a Depth First Search (DFS) on the suffix tree
     * executing a callback on each node
     */
    pub fn dfs(&self, callback: &mut dyn FnMut(TreeNode)) {
        let mut stack = vec![self.nodes[0].as_ref().expect("Root node not found")];

        while let Some(node) = stack.pop() {
            debug!("DFS: {}", node.id);
            callback(node.clone());
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
    ) {
        // if it doesn't have a parent, panic
        let node_ref = &mut self.nodes[node]
            .as_mut()
            .expect("Node not found in break_edge");

        debug!("Breaking edge at index {}", break_idx);

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
            self.create_internal_node(parent, node, original_label_start, break_idx);

        // add the leaf on the new internal node
        self.create_leaf(new_internal_node, leaf_start, leaf_end);
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
    ) -> usize {
        let parent_ref = self.nodes[parent].as_ref().expect("Parent node not found");

        let internal_id = self.last_internal_id;
        debug!(
            "Creating internal node with id: {} and label: {}",
            internal_id,
            &self.original_string[edge_start..edge_end]
        );

        let internal_node = TreeNode {
            id: internal_id,
            string_depth: parent_ref.string_depth + (edge_end - edge_start),
            parent: Some(parent),
            children: vec![None; self.alphabet.len()],
            edge_start: edge_start,
            edge_end: edge_end,
        };

        self.last_internal_id += 1;

        // add new internal node as a child of the parent of the original node
        self.add_child(parent, internal_node);

        // add original node as a child of the new internal node
        self.add_child(
            internal_id,
            self.nodes[original_node].as_ref().unwrap().clone(),
        );

        return internal_id;
    }

    /**
     * Creates a new leaf node with the given parent and edge label
     */
    pub fn create_leaf(&mut self, parent: usize, edge_start: usize, edge_end: usize) -> usize {
        let parent_ref = self.nodes[parent].as_ref().expect("Parent node not found");
        let leaf_id: usize = self.last_leaf_id;
        debug!(
            "Creating leaf node with id: {} and label: {}",
            leaf_id,
            &self.original_string[edge_start..edge_end]
        );

        let leaf = TreeNode {
            id: leaf_id,
            string_depth: parent_ref.string_depth + (edge_end - edge_start),
            parent: Some(parent),
            children: vec![None; self.alphabet.len()],
            edge_start: edge_start,
            edge_end: edge_end,
        };

        self.last_leaf_id += 1;

        self.add_child(parent, leaf);

        return leaf_id;
    }

    /**
     * Walks down the tree and inserts new leaf for the given suffix
     */
    pub fn find_path(&mut self, suffix_idx: usize) {
        let mut current_node = self.nodes[0].as_ref().expect("Root node not found");
        let suffix_len = self.original_string.len() - suffix_idx;

        // how far we've already walked down the suffix
        let mut suffix_sub_idx = 0;

        loop {
            // walk down label on current node's edge
            debug!(
                "Edge label: {}",
                &self.original_string[current_node.edge_start..current_node.edge_end]
            );
            for label_idx in current_node.edge_start..current_node.edge_end {
                if suffix_sub_idx > suffix_len {
                    debug!("Reached end of suffix");
                    break;
                }

                let suffix_char = self.original_string.as_bytes()[suffix_idx + suffix_sub_idx];
                let c = self.original_string.as_bytes()[label_idx];

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
                        self.original_string.len(),
                    );
                    return;
                }

                suffix_sub_idx += 1;
            }

            debug!("Done with edge of current node {}", current_node.id);

            if suffix_sub_idx >= suffix_len {
                return;
            }

            // if we've reached the end of the edge label, move to the next node
            // compare the first character of the edge label of each child node
            // to the next character in the suffix
            let c = self.original_string.as_bytes()[suffix_idx + suffix_sub_idx] as char;

            debug!(
                "Checking for child with {}, suffix index: {}",
                c, suffix_sub_idx
            );

            let child_idx = get_child_index(&self.alphabet, c);
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
                        self.original_string.len(),
                    );
                    return;
                }
            }
        }
    }
}

#[cfg(test)]
mod test {

    use crate::sequence::{SequenceContainer, SequenceOperations};

    use super::SuffixTree;

    // #[test]
    // fn test_tree_simple() {
    //     let tree = SuffixTree::new("A", "alphabets/dna.txt");

    //     assert_eq!(tree.suffixes.len(), 2);
    // }

    // #[test]
    // fn test_tree_simple2() {
    //     let tree = SuffixTree::new("ACA", "alphabets/dna.txt");

    //     assert_eq!(tree.suffixes.len(), 4);
    // }

    #[test]
    fn test_tree_simple3() {
        let tree = SuffixTree::new("BANANA", "alphabets/banana.txt");

        println!("{}", tree);

        assert_eq!(tree.suffixes.len(), 7);
        assert_eq!(tree.stats.num_internal, 3);
        assert_eq!(tree.stats.num_leaves, 7);
        assert_eq!(tree.stats.num_nodes, 10);
        assert_eq!(tree.stats.average_string_depth, 2.0);
        assert_eq!(tree.stats.max_string_depth, 3);
        assert_eq!(tree.stats.bwt, "ANNB$AA".to_string());
    }

    // #[test]
    // fn test_tree_complex() {
    //     let mut sequence_container: SequenceContainer = SequenceContainer {
    //         sequences: Vec::new(),
    //     };

    //     sequence_container.from_fasta("test_data/Covid_Wuhan.fasta");

    //     let suffix_tree = SuffixTree::new(
    //         &sequence_container.sequences[0].sequence,
    //         "alphabets/dna.txt",
    //     );

    //     println!("{}", suffix_tree);
    // }
}
