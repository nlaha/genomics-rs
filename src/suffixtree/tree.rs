use std::cell::{Ref, RefCell};
use std::panic;
use std::rc::Rc;
use std::thread::current;

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
    parent: RefCell<Option<Rc<TreeNode>>>,
    children: RefCell<[Option<Rc<TreeNode>>; 5]>,
    edge: RefCell<Edge>,
}

impl TreeNode {
    pub fn set_parent(&self, parent: Rc<TreeNode>) {
        *self.parent.borrow_mut() = Some(parent);
    }

    pub fn add_child(self: &Rc<Self>, child: Rc<TreeNode>) {
        let child_idx = match child.edge.borrow().label.chars().next() {
            Some(c) => get_child_index(c),
            None => {
                panic!("Child node has no edge label");
            }
        };

        self.children.borrow_mut()[child_idx] = Some(child.clone());
        child.set_parent(self.clone());
    }
}

pub struct SuffixTree {
    // nodes are stored in a flat vector
    original_string: String,
    last_id: usize,
    suffixes: Vec<String>,
    root: RefCell<Rc<TreeNode>>,
}

pub fn get_child_index(c: char) -> usize {
    match c {
        '$' => 0,
        'a' => 1,
        'c' => 2,
        'g' => 3,
        't' => 4,
        _ => {
            panic!("Unsupported character '{}' in edge label", c);
        }
    }
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
            suffixes: Vec::with_capacity(string_length),
            root: RefCell::new(Rc::new(TreeNode {
                id: 0,
                string_depth: 0,
                parent: RefCell::new(None),
                children: RefCell::new([None, None, None, None, None]),
                edge: RefCell::new(Edge {
                    label: "".to_string(),
                }),
            })),
        };

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
    pub fn break_edge(&mut self, node: Rc<TreeNode>, break_idx: usize, leaf_label: &str) {
        // if it doesn't have a parent, panic
        let parent = node
            .parent
            .borrow()
            .clone()
            .expect("Node has no parent - cannot break edge");

        // Get the current label
        let current_label = node.edge.borrow().label.clone();

        // update the current node's edge
        *node.edge.borrow_mut() = Edge {
            label: current_label[break_idx..].to_string(),
        };

        // create a new internal node
        let next_id: usize = self.get_next_node_id();
        let new_internal_node = Rc::new(TreeNode {
            id: next_id,
            string_depth: parent.string_depth + break_idx,
            parent: node.parent.clone(),
            children: RefCell::new([None, None, None, None, None]),
            edge: RefCell::new(Edge {
                label: node.edge.borrow().label[..break_idx].to_string(),
            }),
        });

        // update the current node's parent
        new_internal_node.add_child(node.clone());

        // add the leaf on the new internal node
        let leaf_node = Rc::new(TreeNode {
            id: next_id + 1,
            string_depth: new_internal_node.string_depth + leaf_label.len(),
            parent: RefCell::new(Some(new_internal_node.clone())),
            children: RefCell::new([None, None, None, None, None]),
            edge: RefCell::new(Edge {
                label: leaf_label.to_string(),
            }),
        });

        new_internal_node.add_child(leaf_node.clone());
    }

    /**
     * Walks down the tree and inserts new leaf for the given suffix
     */
    pub fn find_path(&mut self, suffix: &str) {
        self.suffixes.push(suffix.to_string());

        let mut current_node = self.root.borrow().clone();

        loop {
            let edge_label: &str = &current_node.edge.borrow().label.clone();

            // enumerate characters in edge label
            for (i, c) in edge_label.chars().enumerate() {
                let suffix_char = suffix.as_bytes()[i] as char;

                // if the suffix character is not equal to the edge character
                // break a new edge and insert a new leaf
                if suffix_char != c {
                    self.break_edge(current_node.clone(), i, &suffix[i..]);
                    return;
                }
            }

            // if we've reached the end of the edge label, move to the next node
            // compare the first character of the edge label of each child node
            // to the next character in the suffix
            let c = suffix.as_bytes()[edge_label.len()] as char;
            println!("{}", suffix);

            let child_node: Option<Rc<TreeNode>> =
                current_node.children.borrow()[get_child_index(c)].clone();
            match child_node {
                // there's already a child node, so we need to keep moving down the tree
                Some(n) => {
                    current_node = n;
                }
                // there's no child node, so we need to add a new leaf
                None => {
                    current_node.add_child(Rc::new(TreeNode {
                        id: self.get_next_node_id(),
                        string_depth: current_node.string_depth + suffix.len(),
                        parent: RefCell::new(Some(current_node.clone())),
                        children: RefCell::new([None, None, None, None, None]),
                        edge: RefCell::new(Edge {
                            label: suffix[edge_label.len()..].to_string(),
                        }),
                    }));
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
    }

    #[test]
    fn test_tree_simple2() {
        let tree = SuffixTree::new("aca");

        assert_eq!(tree.suffixes.len(), 3);
    }
}
