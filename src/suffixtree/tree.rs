use slab::Slab;
use std::panic;

/**
 * Indicates the start and end of an edge in the suffix tree
 */
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edge {
    label: String,
}

/**
 * Represents a node in the suffix tree
 */
#[derive(Debug, Clone, PartialEq)]
pub struct TreeNode {
    id: usize,
    string_depth: usize,
    parent: Option<usize>,
    children: [Option<usize>; 5],
    edge: Edge,
}

pub struct SuffixTree {
    // nodes are stored in a flat vector
    original_string: String,
    last_id: usize,
    nodes: Vec<Option<TreeNode>>,
    suffixes: Vec<String>,
}

impl SuffixTree {
    /**
     * Create a new suffix tree with a given number of children
     */
    pub fn new(original_string: &str) -> SuffixTree {
        let string_length = original_string.len();

        let mut tree = SuffixTree {
            original_string: original_string.to_string(),
            last_id: string_length,
            nodes: vec![None; string_length * 6],
            suffixes: Vec::with_capacity(string_length),
        };

        // insert root
        let root = TreeNode {
            id: 0,
            string_depth: 0,
            parent: None,
            children: [None, None, None, None, None],
            edge: Edge {
                label: "".to_string(),
            },
        };

        tree.nodes.push(Some(root));

        // build a set of suffixes with '$' appended to the end
        for i in 0..string_length {
            let suffix = tree.original_string[i..].to_string() + "$";
            tree.find_path(suffix.as_str());
        }

        return tree;
    }

    /**
     * Increments the ID counter and returns the ID to be
     * used for an internal node
     */
    pub fn get_next_node_id(&mut self) -> usize {
        self.last_id += 6;

        return self.last_id;
    }

    /**
     * Breaks the edge at the given index and inserts a new node
     * also populates the leaf node with the given label
     */
    pub fn break_edge(&mut self, node_idx: usize, break_idx: usize, leaf_label: &str) {
        let node = &mut self.nodes[node_idx].clone().expect("Node is null");
        // if it doesn't have a parent, panic
        if node.parent.is_none() {
            panic!("Node has no parent - cannot break edge");
        }

        let new_internal_node = TreeNode {
            id: 0,
            string_depth: node.parent.unwrap().clone() + break_idx,
            parent: node.parent.clone(),
            children: [None, None, None, None, None],
            edge: Edge {
                label: node.edge.label[..break_idx].to_string(),
            },
        };

        // update the current node's edge
        node.edge.label = node.edge.label[break_idx..].to_string();

        // update the current node's parent
        node.parent = Some(new_internal_node.id);

        self.add_leaf(leaf_label, new_internal_node);
    }

    /**
     * Adds a leaf node to the internal node
     */
    fn add_leaf(&mut self, leaf_label: &str, internal_node: TreeNode) {
        // create a leaf node
        let new_leaf_node = TreeNode {
            id: 0,
            string_depth: internal_node.string_depth + leaf_label.len(),
            parent: Some(internal_node.id),
            children: [None, None, None, None, None],
            edge: Edge {
                label: leaf_label.to_string(),
            },
        };

        // add the new leaf node to the internal node
        self.add_child(&internal_node, new_leaf_node);
    }

    /**
     * Adds a child with a given edge to the parent node
     * at the correct child index based on the alphabetic order of the first character
     */
    pub fn add_child(&mut self, parent: &TreeNode, mut child: TreeNode) {
        let next_id: usize = self.get_next_node_id();
        child.id = next_id;
        match child.edge.label.chars().next() {
            Some(c) => match c {
                '$' => {
                    self.nodes[parent.id + 1] = Some(child);
                }
                'a' => {
                    self.nodes[parent.id + 2] = Some(child);
                }
                'c' => {
                    self.nodes[parent.id + 3] = Some(child);
                }
                'g' => {
                    self.nodes[parent.id + 4] = Some(child);
                }
                't' => {
                    self.nodes[parent.id + 5] = Some(child);
                }
                _ => {
                    panic!("Unsupported character '{}' in edge label", c);
                }
            },
            None => {
                panic!("Child node has no edge label");
            }
        }
    }

    /**
     * Walks down the tree and inserts new leaf for the given suffix
     */
    pub fn find_path(&mut self, suffix: &str) {
        let mut node_index = 0; // Start with the root node index

        self.suffixes.push(suffix.to_string());

        loop {
            let current_node = self.nodes[node_index]
                .clone()
                .expect(format!("Node is null: {}", node_index).as_str());
            let edge_label: &str = &current_node.edge.label;

            // enumerate characters in edge label
            for (i, c) in edge_label.chars().enumerate() {
                let suffix_char = suffix.as_bytes()[i] as char;

                // if the suffix character is not equal to the edge character
                // break a new edge and insert a new leaf
                if suffix_char != c {
                    self.break_edge(node_index, i, &suffix[i..]);
                    return;
                }
            }

            // if we've reached the end of the edge label, move to the next node
            // compare the first character of the edge label of each child node
            // to the next character in the suffix
            let c = suffix.as_bytes()[edge_label.len()] as char;
            println!("{}", suffix);
            let child_offset = match c {
                '$' => 0,
                'a' => 1,
                'c' => 2,
                'g' => 3,
                't' => 4,
                _ => {
                    panic!("Unsupported character '{}' in edge label", c);
                }
            };

            let child_node: Option<usize> = current_node.children[child_offset];
            match child_node {
                // there's already a child node, so we need to keep moving down the tree
                Some(n) => {
                    node_index = n;
                }
                // there's no child node, so we need to add a new leaf
                None => {
                    self.add_leaf(&suffix[edge_label.len()..], current_node.clone());
                    return;
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::SuffixTree;

    #[test]
    fn test_tree_simple() {
        let tree = SuffixTree::new("a");

        assert!(tree.nodes.len() == 2);
        assert!(tree.suffixes.len() == 1);
        assert!(tree.nodes[1].clone().unwrap().edge.label == "a");
        assert!(tree.nodes[1].clone().unwrap().string_depth == 1);
        assert!(tree.nodes[1].clone().unwrap().parent == Some(0));
    }
}
