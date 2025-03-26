use std::{env, fmt::Display};

use log::warn;
use petgraph::{dot::Dot, Graph};

use super::tree::{SuffixTree, TreeNode, TreeStats};

impl Display for TreeStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "
            BWT: {}
            BWT Length: {}
            Internal nodes: {}
            Leaves: {}
            Nodes: {}
            Average string depth: {}
            Max string depth: {}
            ",
            if self.bwt.len() > 100 {
                self.bwt[..100].to_string() + "... (truncated)"
            } else {
                self.bwt.to_string()
            },
            self.bwt.len(),
            self.num_internal,
            self.num_leaves,
            self.num_nodes,
            self.average_string_depth,
            self.max_string_depth
        )
    }
}

impl Display for SuffixTree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // if the number of nodes is < 50, compute graphviz
        if self.nodes.len() < 50 {
            let dot = self.write_graphviz();
            writeln!(f, "Graphviz:\n {}", dot)?;
        } else {
            warn!("Graphviz output is too large to display.",);
        }

        if env::var("RUST_LOG") == Ok("debug".to_string()) {
            self.display_string_depth(f);
        }

        writeln!(f, "\nStats: {}", self.stats)
    }
}

impl SuffixTree {
    pub fn write_graphviz(&self) -> String {
        let mut graph = Graph::<String, String>::new();
        self.dfs(&mut |node: TreeNode| {
            let node_idx = graph.add_node(node.id.to_string());

            match node.parent {
                Some(parent) => {
                    let parent_ref = self.nodes[parent].as_ref().unwrap();
                    let parent_idx = graph
                        .node_indices()
                        .find(|idx| {
                            let node: &String = graph.node_weight(*idx).unwrap();
                            node == &parent_ref.id.to_string()
                        })
                        .unwrap();
                    graph.add_edge(
                        parent_idx,
                        node_idx,
                        self.original_string[node.edge_start..node.edge_end].to_string(),
                    );
                }
                None => {
                    // root node case
                }
            }
        });

        // iterate through nodes and create suffix links
        for c_node in self.nodes.iter() {
            if let Some(c_node) = c_node {
                if let Some(suffix_link) = c_node.suffix_link {
                    let suffix_link_ref = self.nodes[suffix_link].as_ref().unwrap();
                    let suffix_link_idx = graph
                        .node_indices()
                        .find(|idx| {
                            let node: &String = graph.node_weight(*idx).unwrap();
                            node == &suffix_link_ref.id.to_string()
                        })
                        .unwrap();
                    let node_idx = graph
                        .node_indices()
                        .find(|idx| {
                            let node: &String = graph.node_weight(*idx).unwrap();
                            node == &c_node.id.to_string()
                        })
                        .unwrap();
                    graph.add_edge(node_idx, suffix_link_idx, "*".to_string());
                }
            }
        }

        format!("{}", Dot::new(&graph))
    }
}
