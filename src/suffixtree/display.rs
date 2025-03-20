use std::{fmt::Display, rc::Rc};

use petgraph::{dot::Dot, Graph};

use super::tree::{SuffixTree, TreeNode, TreeStats};

impl Display for TreeStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "
            Internal nodes: {}
            Leaves: {}
            Nodes: {}
            Average string depth: {}
            Max string depth: {}
            BWT: {}
            ",
            self.num_internal,
            self.num_leaves,
            self.num_nodes,
            self.average_string_depth,
            self.max_string_depth,
            self.bwt
        )
    }
}

impl Display for SuffixTree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // dfs traversal of the tree
        let mut graph = Graph::<String, String>::new();
        let mut parent: Option<Rc<TreeNode>> = None;
        self.dfs(&mut |node: Rc<TreeNode>| {
            let node_idx = graph.add_node(node.id.to_string());
            parent = node.parent.borrow().clone();

            match &parent {
                Some(parent) => {
                    let parent_idx = graph
                        .node_indices()
                        .find(|idx| {
                            let node: &String = graph.node_weight(*idx).unwrap();
                            node == &parent.id.to_string()
                        })
                        .unwrap();
                    graph.add_edge(parent_idx, node_idx, node.edge.borrow().label.to_string());
                }
                None => {
                    // root node case
                }
            }
        });

        writeln!(
            f,
            "Graphviz:\n{:?}",
            Dot::with_config(&graph, &[petgraph::dot::Config::NodeNoLabel])
        )?;

        writeln!(f, "Stats: {}", self.stats)
    }
}
