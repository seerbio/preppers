pub mod io;
pub mod fasta;

use std::cmp::min;
// Reexports
pub use fasta::read_fasta;

use blart::visitor::{TreeStats, TreeStatsCollector};
use blart::{AttemptOptimisticPrefixMatch, ConcreteNodePtr, InnerNode, LeafNode, Node, NodePtr, NodeType, OpaqueNodePtr, OptimisticMismatch, PessimisticMismatch, PrefixMatch, TreeMap};
use std::ffi::CString;

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

fn annotate_sequence<const N: usize>(root: &OpaqueNodePtr<CString, PeptideId, N>, seq: &[u8]) -> Vec<PeptideId> {
    let mut res = Vec::new();

    for i in 0..seq.len() {
        // SAFETY: Since we have an immutable reference to the root node, that
        // means there can only exist other immutable references aside from this one,
        // and no mutable references. That means that no mutating operations can occur
        // on the root node or any child of the root node.
        unsafe {
            _annotate_sequence(root, seq, i, &mut res);
        }
    }

    res
}

/// Search in the given tree for any peptides that are substrings of `seq` starting at
/// `start`.
///
/// This is essentially equivalent to `search_unchecked` but permits recording multiple
/// matches while advancing through the sequence.
///
/// # Safety
///  - This function cannot be called concurrently with any mutating operation
///    on `root` or any child node of `root`. This function will arbitrarily
///    read to any child in the given tree.
unsafe fn _annotate_sequence<const N: usize>(root: &OpaqueNodePtr<CString, PeptideId, N>, seq: &[u8], start: usize, res: &mut Vec<PeptideId>) {
    let mut depth = 0;
    let mut pessimistic_depth = 0;
    let mut current_node = *root;
    loop {
        let step = match current_node.to_node_ptr() {
            ConcreteNodePtr::Node4(t) => unsafe {
                // SAFETY: The safety requirement is covered by the safety requirement on the
                // containing function
                do_inner_lookup(seq, start, depth, pessimistic_depth, res, t)
            }
            ConcreteNodePtr::Node16(t) => unsafe {
                // SAFETY: The safety requirement is covered by the safety requirement on the
                // containing function
                do_inner_lookup(seq, start, depth, pessimistic_depth, res, t)
            }
            ConcreteNodePtr::Node48(t) => unsafe {
                // SAFETY: The safety requirement is covered by the safety requirement on the
                // containing function
                do_inner_lookup(seq, start, depth, pessimistic_depth, res, t)
            }
            ConcreteNodePtr::Node256(t) => unsafe {
                // SAFETY: The safety requirement is covered by the safety requirement on the
                // containing function
                do_inner_lookup(seq, start, depth, pessimistic_depth, res, t)
            }
            ConcreteNodePtr::LeafNode(_) => {
                panic!("Encountered leaf unexpectedly!")
            }
        };

        if step.is_none() {
            break
        }

        let any_implicit: bool;
        (current_node, depth, any_implicit) = step.unwrap();

        if !any_implicit {
            pessimistic_depth = depth;
        }
    }
}

/// Handle this inner node as a potential match for the given start and depth.
/// Will perform the following:
/// - Check for explicit or implicit prefix match to the node, and return None on mismatch
/// - Check for keys terminating after the node's prefix (child for byte \0) and handle the leaf
/// - Locate any next node further in the sequence
///
/// # Safety
///  - No other access or mutation to the `t` Node can happen while this function runs.
unsafe fn do_inner_lookup<T: Node<N, Key=CString, Value=PeptideId> + InnerNode<N>, const N: usize>(seq: &[u8], start: usize, depth: usize, pessimistic_depth: usize, res: &mut Vec<PeptideId>, t: NodePtr<N, T>)
    -> Option<(OpaqueNodePtr<CString, PeptideId, N>, usize, bool)>
{
    // SAFETY: The lifetime produced from this is bounded to this scope and does not
    // escape. Further, no other code mutates the node referenced, which is further
    // enforced the "no concurrent reads or writes" requirement on the
    // `_annotate_sequence` function.
    let inner_node = unsafe { t.as_ref() };

    // println!("attempting match at depth {} ({})", depth, &seq.as_str()[start + depth..]);

    let trunc_seq = &seq[start + depth..];

    let matched_bytes: usize;
    let was_optimistic: bool;
    if pessimistic_depth < depth {
        // We must switch to optimistic matching
        match inner_node.optimistic_match_prefix(trunc_seq) {
            Ok(m) => {
                matched_bytes = m.matched_bytes;
                was_optimistic = true;
            }
            Err(_) => {
                return None
            }
        }
    } else {
        match inner_node.attempt_pessimistic_match_prefix(trunc_seq) {
            Ok(m) => {
                matched_bytes = m.matched_bytes;
                was_optimistic = m.any_implicit_bytes;
            }
            Err(_) => {
                return None
            }
        }
    }

    let new_depth = depth + matched_bytes;

    if start + new_depth >= seq.len() {
        // println!("too close to end at depth {}", depth);
        return None
    }

    // println!("matched {} bytes (prefix: {})", m.matched_bytes, &seq[start..start + new_depth]);

    // Two possible paths from here -- either a string terminating zero, or the next char
    let check_start = if was_optimistic {
        pessimistic_depth
    } else {
        new_depth
    };
    inner_node.lookup_child(b'\0').map(
        |l| match l.to_node_ptr() {
            ConcreteNodePtr::LeafNode(n) => {
                // SAFETY: The safety requirement is covered by the safety requirement on the
                // containing function
                unsafe { handle_leaf(&seq[start..], res, n, check_start); }
            }
            _ => { panic!("Found non-leaf for null byte!") }
        }
    );

    let next = inner_node.lookup_child(seq[start + new_depth]);

    if next.is_none() {
        return None
    }

    let n = next.unwrap();

    // Check node type before calling `to_node_ptr`, as converting back
    // to opaque (to return) is expensive!
    match n.node_type() {
        NodeType::Leaf => {
            match n.to_node_ptr() {
                ConcreteNodePtr::LeafNode(l) => {
                    // SAFETY: The safety requirement is covered by the safety requirement on the
                    // containing function
                    unsafe { handle_leaf(&seq[start..], res, l, check_start); }
                    None
                }
                _ => { panic!("Unwrapped unexpected node type!") }
            }
        }
        _ => {
            // Plus one, as we matched one additional byte with lookup_child
            Some((n, new_depth + 1, was_optimistic))
        }
    }
}

///
/// # Safety
///  - No other access or mutation to the `t` Node can happen while this function runs.
unsafe fn handle_leaf<const N: usize>(seq: &[u8], res: &mut Vec<PeptideId>, t: NodePtr<N, LeafNode<CString, PeptideId, N>>, check_start: usize) {
    // Due to the potential for optimistic matching, we need to check the full
    // key matches. In the future, we can track pessimistic/optimistic match
    // to elide this check if all matches to the prefix were explicit.

    // SAFETY: The lifetime produced from this is bounded to this scope and does not
    // escape. Further, no other code mutates the node referenced, which is further
    // enforced the "no concurrent reads or writes" requirement on the
    // `_annotate_sequence` function.
    let inner_node = unsafe { t.as_ref() };

    let key = inner_node.key_ref();
    let key_bytes = key.as_bytes_with_nul();

    let mut i = check_start;
    while i < seq.len() {
        if key_bytes[i] == b'\0' {
            // All bytes matched
            res.push(*inner_node.value_ref());
        }

        if key_bytes[i] != seq[i] {
            // Mismatch
            return
        }

        i += 1;
    }

    // No more bytes in the sequence before exhausting the key; no match
}

#[cfg(test)]
mod tests {

    // #[test]
    // fn it_works() {
    //     let result = add(2, 2);
    //     assert_eq!(result, 4);
    // }
}
