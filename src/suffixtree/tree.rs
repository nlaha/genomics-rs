use std::cell::RefCell;
use std::panic;
use std::rc::Rc;

use log::{debug, info};

/**
 * Indicates the start and end of an edge in the suffix tree
 */
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edge {
    pub label: String,
}

/**
 * Represents a node in the suffix tree
 */
#[derive(Debug, Clone, PartialEq)]
pub struct TreeNode {
    pub id: usize,
    pub string_depth: usize,
    pub parent: RefCell<Option<Rc<TreeNode>>>,
    pub children: RefCell<Vec<Option<Rc<TreeNode>>>>,
    pub edge: RefCell<Edge>,
}

impl TreeNode {
    pub fn set_parent(&self, parent: Rc<TreeNode>) {
        *self.parent.borrow_mut() = Some(parent);
    }

    pub fn add_child(self: &Rc<Self>, alphabet: &Vec<char>, child: Rc<TreeNode>) {
        let child_idx = match child.edge.borrow().label.chars().next() {
            Some(c) => get_child_index(alphabet, c),
            None => {
                panic!("Child node has no edge label");
            }
        };

        self.children.borrow_mut()[child_idx] = Some(child.clone());
        child.set_parent(self.clone());
    }
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
    pub root: RefCell<Rc<TreeNode>>,
    pub alphabet: Vec<char>,
    pub internal_nodes: Vec<Rc<TreeNode>>,
    pub leaf_nodes: Vec<Rc<TreeNode>>,
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
            root: RefCell::new(Rc::new(TreeNode {
                id: 0,
                string_depth: 0,
                parent: RefCell::new(None),
                children: RefCell::new(vec![None; alphabet.len()]),
                edge: RefCell::new(Edge {
                    label: "".to_string(),
                }),
            })),
            alphabet: alphabet_sorted,
            internal_nodes: Vec::new(),
            leaf_nodes: Vec::new(),
            stats: TreeStats {
                num_internal: 0,
                num_leaves: 0,
                num_nodes: 0,
                average_string_depth: 0.0,
                max_string_depth: 0,
                bwt: "".to_string(),
            },
        };

        // build a set of suffixes with '$' appended to the end
        for i in 0..string_length + 1 {
            let suffix = tree.original_string[i..].to_string();
            // if suffix is longer than 100 characters, truncate it
            if suffix.len() > 100 {
                info!("[FindPath] {}/{} {}...", i, string_length, &suffix[..100]);
            } else {
                info!("[FindPath] {}/{} {}", i, string_length, suffix);
            }
            tree.find_path(suffix.as_str());
        }

        tree.stats.num_nodes = tree.leaf_nodes.len() + tree.internal_nodes.len();
        tree.stats.num_internal = tree.internal_nodes.len();
        tree.stats.num_leaves = tree.leaf_nodes.len();
        tree.stats.average_string_depth = tree
            .internal_nodes
            .iter()
            .map(|n| n.string_depth)
            .sum::<usize>() as f64
            / tree.internal_nodes.len() as f64;
        tree.stats.max_string_depth = tree
            .internal_nodes
            .iter()
            .map(|n| n.string_depth)
            .max()
            .unwrap_or(0);

        // compute burrows-wheeler transform
        let mut bwt: Vec<char> = vec![' '; tree.original_string.len()];
        let mut idx = 0;
        tree.dfs(&mut |node: Rc<TreeNode>| {
            // if it's a leaf
            if node.id > 0 && node.id < tree.suffixes.len() + 1 && idx < tree.original_string.len()
            {
                if node.id == 1 {
                    bwt[idx] = '$';
                } else {
                    bwt[idx] = tree.original_string.as_bytes()[node.id - 2] as char;
                }
                idx += 1;
            }
        });
        tree.stats.bwt = bwt.iter().collect::<String>().trim().to_string();

        return tree;
    }

    /**
     * Performs a Depth First Search (DFS) on the suffix tree
     * executing a callback on each node
     */
    pub fn dfs(&self, callback: &mut dyn FnMut(Rc<TreeNode>)) {
        let mut stack = vec![self.root.borrow().clone()];

        while let Some(node) = stack.pop() {
            callback(node.clone());
            for child in node.children.borrow().iter().rev().flatten() {
                stack.push(child.clone());
            }
        }
    }

    /**
     * Breaks the edge at the given index and inserts a new node
     * also populates the leaf node with the given label
     */
    pub fn break_edge(&mut self, node: Rc<TreeNode>, break_idx: usize, leaf_label: &str) {
        // if it doesn't have a parent, panic
        let parent = node
            .parent
            .borrow()
            .clone()
            .expect("Node has no parent - cannot break edge");

        debug!("Breaking edge at {}", break_idx);

        // Get the current label
        let current_label = node.edge.borrow().label.to_owned();

        // update the current node's edge
        *node.edge.borrow_mut() = Edge {
            label: current_label[break_idx..].to_string().to_owned(),
        };

        // create a new internal node
        let new_internal_node =
            self.create_internal_node(parent.clone(), &current_label[..break_idx]);

        // update the current node's parent
        new_internal_node.add_child(&self.alphabet, node.clone());

        // add the leaf on the new internal node
        self.create_leaf(new_internal_node.clone(), leaf_label);
    }

    /**
     * Creates a new internal node with the given parent and edge label
     */
    pub fn create_internal_node(&mut self, parent: Rc<TreeNode>, label: &str) -> Rc<TreeNode> {
        let internal_id = self.last_internal_id;
        debug!(
            "Creating internal node with id: {} and label: {}",
            internal_id, label
        );

        let internal_node = Rc::new(TreeNode {
            id: internal_id,
            string_depth: parent.string_depth + label.bytes().len(),
            parent: RefCell::new(Some(parent.clone())),
            children: RefCell::new(vec![None; self.alphabet.len()]),
            edge: RefCell::new(Edge {
                label: label.to_string(),
            }),
        });

        self.internal_nodes.push(internal_node.clone());
        self.last_internal_id = internal_id + 1;

        parent.add_child(&self.alphabet, internal_node.clone());

        return internal_node;
    }

    /**
     * Creates a new leaf node with the given parent and edge label
     */
    pub fn create_leaf(&mut self, parent: Rc<TreeNode>, label: &str) -> Rc<TreeNode> {
        let leaf_id: usize = self.last_leaf_id;
        debug!(
            "Creating leaf node with id: {} and label: {}",
            leaf_id, label
        );

        let leaf = Rc::new(TreeNode {
            id: leaf_id,
            string_depth: parent.string_depth + label.bytes().len(),
            parent: RefCell::new(Some(parent.clone())),
            children: RefCell::new(vec![None; self.alphabet.len()]),
            edge: RefCell::new(Edge {
                label: label.to_string(),
            }),
        });

        self.leaf_nodes.push(leaf.clone());
        self.last_leaf_id = leaf_id + 1;

        parent.add_child(&self.alphabet, leaf.clone());

        return leaf;
    }

    /**
     * Walks down the tree and inserts new leaf for the given suffix
     */
    pub fn find_path(&mut self, suffix: &str) {
        self.suffixes.push(suffix.to_string());

        let mut current_node = self.root.borrow().clone();

        let mut suffix_idx = 0;
        loop {
            let edge_label: &str = &current_node.edge.borrow().label.clone();

            // enumerate characters in edge label
            if edge_label.len() > 100 {
                debug!("Edge label: {}...", edge_label[0..100].to_string());
            } else {
                debug!("Edge label: {}", edge_label);
            }

            for (i, c) in edge_label.bytes().enumerate() {
                if (suffix_idx + i) >= suffix.len() {
                    break;
                }

                let suffix_char = suffix.as_bytes()[i + suffix_idx];

                // if the suffix character is not equal to the edge character
                // break a new edge and insert a new leaf
                debug!("Comparing {} to {}", suffix_char as char, c as char);
                if suffix_char != c {
                    self.break_edge(current_node.clone(), i, &suffix[i..]);
                    return;
                }
            }

            if suffix_idx >= suffix.len() {
                return;
            }

            // if we've reached the end of the edge label, move to the next node
            // compare the first character of the edge label of each child node
            // to the next character in the suffix
            if edge_label.len() > 0 {
                suffix_idx += edge_label.len();
            }
            let c = suffix.as_bytes()[suffix_idx] as char;

            debug!("Suffix index: {}", suffix_idx);

            let child_idx = get_child_index(&self.alphabet, c);
            let child_node: Option<Rc<TreeNode>> =
                current_node.children.borrow()[child_idx].clone();

            debug!(
                "Finding next child node for {}, child idx: {}",
                c, child_idx
            );

            match child_node {
                // there's already a child node, so we need to keep moving down the tree
                Some(n) => {
                    if n.edge.borrow().label.len() > 100 {
                        debug!(
                            "Found child node: {}, label: {}...",
                            n.id,
                            n.edge.borrow().label[0..100].to_string()
                        );
                    } else {
                        debug!(
                            "Found child node: {}, label: {}",
                            n.id,
                            n.edge.borrow().label
                        );
                    }
                    current_node = n;
                }
                // there's no child node, so we need to add a new leaf
                None => {
                    debug!("No child node found, adding new leaf");
                    current_node.add_child(
                        &self.alphabet.clone(),
                        self.create_leaf(
                            current_node.clone(),
                            &suffix[edge_label.len()..].to_string(),
                        ),
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

    // #[test]
    // fn test_tree_simple3() {
    //     let tree = SuffixTree::new("BANANA", "alphabets/banana.txt");

    //     println!("{}", tree);

    //     assert_eq!(tree.suffixes.len(), 7);
    //     assert_eq!(tree.stats.num_internal, 3);
    //     assert_eq!(tree.stats.num_leaves, 7);
    //     assert_eq!(tree.stats.num_nodes, 10);
    //     assert_eq!(tree.stats.average_string_depth, 2.0);
    //     assert_eq!(tree.stats.max_string_depth, 3);
    //     assert_eq!(tree.stats.bwt, "ANNB$AA".to_string());
    // }

    #[test]
    fn test_tree_complex() {
        let mut sequence_container: SequenceContainer = SequenceContainer {
            sequences: Vec::new(),
        };

        sequence_container.from_fasta("test_data/Covid_Wuhan.fasta");

        let suffix_tree = SuffixTree::new(
            &sequence_container.sequences[0].sequence,
            "alphabets/dna.txt",
        );

        println!("{}", suffix_tree);
    }
}
