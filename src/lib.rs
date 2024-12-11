pub mod io;
pub mod fasta;

// Reexports
pub use fasta::read_fasta;

use blart::visitor::{TreeStats, TreeStatsCollector};
use blart::{ConcreteNodePtr, InnerNode, LeafNode, Node, NodePtr, OpaqueNodePtr, TreeMap};
use std::ffi::CString;
use std::ops::Index;

// Types
type PeptideId = u64;

#[derive(Debug)]
pub struct PeptideTrie {
    _tree: TreeMap<CString, PeptideId>,
    _next_id: PeptideId,
}

impl PeptideTrie {
    pub fn new() -> PeptideTrie {
        PeptideTrie {
            _tree: TreeMap::new(),
            _next_id: 0,
        }
    }

    pub fn add(&mut self, peptide: &[u8]) {
        let id = self._next_id;
        self._next_id = id + 1;
        self._tree.insert(CString::new(peptide).expect("Invalid string!"), id);
    }

    pub fn len(&self) -> usize {
        self._tree.len()
    }

    pub fn stats(&self) -> Option<TreeStats> {
        TreeStatsCollector::collect(&self._tree)
    }
}

fn annotate_sequence<const N: usize>(node: &OpaqueNodePtr<CString, PeptideId, N>, seq: &[u8]) -> Vec<PeptideId> {
    let mut res = Vec::new();

    for i in 0..seq.len() {
        _annotate_sequence(node, seq, i, &mut res);
    }

    res
}

/// Search in the given tree for any peptides that are substrings of `seq` starting at
/// `start`.
///
/// This is essentially equivalent to `search_unchecked` but permits recording multiple
/// matches while advancing through the sequence. We use the same adaptive prefix
/// matching approach, but unconditionally compare full key strings (optimistic strategy),
/// which will be slightly less efficient for sufficiently short peptides.
fn _annotate_sequence<const N: usize>(node: &OpaqueNodePtr<CString, PeptideId, N>, seq: &[u8], start: usize, res: &mut Vec<PeptideId>) {
    let mut depth = 0;
    let mut current_node = *node;
    loop {
        let step = match current_node.to_node_ptr() {
            ConcreteNodePtr::Node4(t) => {
                do_inner_lookup(seq, start, depth, res, t)
            }
            ConcreteNodePtr::Node16(t) => {
                do_inner_lookup(seq, start, depth, res, t)
            }
            ConcreteNodePtr::Node48(t) => {
                do_inner_lookup(seq, start, depth, res, t)
            }
            ConcreteNodePtr::Node256(t) => {
                do_inner_lookup(seq, start, depth, res, t)
            }
            ConcreteNodePtr::LeafNode(_) => {
                panic!("Encountered leaf unexpectedly!")
            }
        };

        if step.is_none() {
            break
        }

        (current_node, depth) = step.unwrap();
    }
}

/// Handle this inner node as a potential match for the given start and depth.
/// Will perform the following:
/// - Check for explicit or implicit prefix match to the node, and return None on mismatch
/// - Check for keys terminating after the node's prefix (child for byte \0) and handle the leaf
/// - Locate any next node further in the sequence
fn do_inner_lookup<T: Node<N, Key=CString, Value=PeptideId> + InnerNode<N>, const N: usize>(seq: &[u8], start: usize, depth: usize, res: &mut Vec<PeptideId>, t: NodePtr<N, T>)
    -> Option<(OpaqueNodePtr<CString, PeptideId, N>, usize)>
{
    let inner_node = t.read();

    // println!("attempting match at depth {} ({})", depth, &seq.as_str()[start + depth..]);

    let prefix_match = inner_node.attempt_pessimistic_match_prefix(&seq[start + depth..]);

    if prefix_match.is_err() {
        // println!("failed to match at depth {} ({})", depth, &seq.as_str()[start + depth..]);
        return None
    }

    let m = prefix_match.unwrap();

    let new_depth = depth + m.matched_bytes;

    if start + new_depth >= seq.len() {
        // println!("too close to end at depth {}", depth);
        return None
    }

    // println!("matched {} bytes (prefix: {})", m.matched_bytes, &seq[start..start + new_depth]);

    // Two possible paths from here -- either a string terminating zero, or the next char
    inner_node.lookup_child(b'\0').map(
        |l| match l.to_node_ptr() {
            ConcreteNodePtr::LeafNode(n) => {
                handle_leaf(&seq, start, res, n);
            }
            _ => { panic!("Found non-leaf for null byte!") }
        }
    );

    let next = inner_node.lookup_child(seq[start + new_depth]);

    if next.is_none() {
        return None
    }

    match next.unwrap().to_node_ptr() {
        ConcreteNodePtr::LeafNode(n) => {
            handle_leaf(&seq, start, res, n);
            None
        }
        n => {
            Some((n.to_opaque(), new_depth + 1))
        }
    }
}

fn handle_leaf<const N: usize>(seq: &[u8], start: usize, res: &mut Vec<PeptideId>, t: NodePtr<N, LeafNode<CString, PeptideId, N>>) {
    // Due to the potential for optimistic matching, we need to check the full
    // key matches. In the future, we can track pessimistic/optimistic match
    // to elide this check if all matches to the prefix were explicit.

    let n = t.read();
    let key = n.key_ref();

    // println!("checking leaf {:?}", key);
    // TODO: change to check equality as we walk through looking for null terminator
    let len = key.as_bytes().len();

    if start + len < seq.len() - 1 && key.as_bytes()[0..len].eq(&seq[start..start + len]) {
        // println!("key matches! adding {:?} to result", key);
        res.push(*t.read().value_ref());
    }
}

#[cfg(test)]
mod tests {

    // #[test]
    // fn it_works() {
    //     let result = add(2, 2);
    //     assert_eq!(result, 4);
    // }
}
