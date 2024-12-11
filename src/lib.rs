pub mod io;
pub mod fasta;

// Reexports
pub use fasta::read_fasta;

use blart::visitor::{TreeStats, TreeStatsCollector};
use blart::{ConcreteNodePtr, InnerNode, Node, NodePtr, OpaqueNodePtr, TreeMap};
// Imports
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
        _annotate_sequence(node, seq, i, 0, &mut res);
    }

    res
}

/// Search in the given tree for any peptides that are substrings of `seq` starting at
/// `start`, having already traversed `depth` characters within the tree.
///
/// This is essentially equivalent to `search_unchecked` but permits recording multiple
/// matches while advancing through the sequence. We use the same adaptive prefix
/// matching approach, but unconditionally compare full key strings (optimistic strategy),
/// which will be slightly less efficient for sufficiently short peptides.
fn _annotate_sequence<const N: usize>(node: &OpaqueNodePtr<CString, PeptideId, N>, seq: &[u8], start: usize, depth: usize, res: &mut Vec<PeptideId>) {
    match node.to_node_ptr() {
        ConcreteNodePtr::Node4(t) => {
            do_inner_lookup(seq, start, depth, res, t);
        }
        ConcreteNodePtr::Node16(t) => {
            do_inner_lookup(seq, start, depth, res, t);
        }
        ConcreteNodePtr::Node48(t) => {
            do_inner_lookup(seq, start, depth, res, t);
        }
        ConcreteNodePtr::Node256(t) => {
            do_inner_lookup(seq, start, depth, res, t);
        }
        ConcreteNodePtr::LeafNode(t) => {
            // Due to the potential for optimistic matching, we need to check the full
            // key matches. In the future, we can track pessimistic/optimistic match
            // to elide this check if all matches to the prefix were explicit.

            let n = t.read();
            let key = n.key_ref();

            // println!("checking leaf {:?}", key);
            // TODO: change to check equality as we walk through looking for null terminator
            let len = key.as_bytes().len();

            if start + len < seq.len()-1 && key.as_bytes()[0..len].eq(&seq[start..start+len]) {
                // println!("key matches! adding {:?} to result", key);
                res.push(*t.read().value_ref());
            }
        }
    }
}

fn do_inner_lookup<T: Node<N, Key=CString, Value=PeptideId> + InnerNode<N>, const N: usize>(seq: &[u8], start: usize, depth: usize, res: &mut Vec<PeptideId>, t: NodePtr<N, T>) -> () {
    let inner_node = t.read();

    // println!("attempting match at depth {} ({})", depth, &seq.as_str()[start + depth..]);

    let prefix_match = inner_node.attempt_pessimistic_match_prefix(&seq[start + depth..]);

    if prefix_match.is_err() {
        // println!("failed to match at depth {} ({})", depth, &seq.as_str()[start + depth..]);
        return
    }

    let m = prefix_match.unwrap();

    let new_depth = depth + m.matched_bytes;

    if start + new_depth >= seq.len() - 1 {
        // println!("too close to end at depth {}", depth);
        return
    }

    // println!("matched {} bytes (prefix: {})", m.matched_bytes, &seq[start..start + new_depth]);

    // Two possible paths from here -- either a string terminating zero, or the next char
    for byte in vec![
        b'\0',
        seq[start + new_depth]
    ] {
        // let brep = if byte != 0 {
        //     format!("{}", byte as char)
        // } else {
        //     "\\0".into()
        // };
        // println!("searching for next byte {} at depth {}", brep, new_depth);

        let next = inner_node.lookup_child(byte);

        if next.is_none() {
            // println!("no match!")
        } else {
            let n = next.unwrap();
            // println!("matched a child node; descending into {:?}", n.to_node_ptr());

            _annotate_sequence(
                &n,
                seq,
                start,
                new_depth + 1, // we matched one additional byte with `lookup_child`
                res
            );
        }
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
