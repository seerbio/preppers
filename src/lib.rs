pub mod io;
pub mod fasta;

// Reexports
pub use fasta::read_fasta;

use blart::visitor::{TreeStats, TreeStatsCollector};
use blart::{ConcreteNodePtr, InnerNode, LeafNode, Node, NodePtr, NodeType, OpaqueNodePtr, TreeMap};
use std::ffi::CString;
use blart::map::EntryRef;

// Types
pub type PeptideId = u64;

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

    pub fn insert(&mut self, peptide: &[u8]) -> PeptideId {
        // Build reversed key to store in the trie
        let key = CString::new(peptide.iter().rev().copied().collect::<Vec<_>>()).expect(&format!("Invalid peptide sequence {:?}", &peptide));

        let entry = self._tree.try_entry_ref(&key);

        match entry {
            Ok(v) => {
                match v {
                    EntryRef::Occupied(o) => {
                        // Return existing ID
                        *o.get()
                    }
                    EntryRef::Vacant(v) => {
                        // Insert a new ID
                        let id = self._next_id;
                        self._next_id = id + 1;

                        v.insert(id);

                        id
                    }
                }
            }
            Err(_) => {
                panic!("Attempted insertion of illegal key. This may be a result of previous corruption!");
            }
        }
    }

    pub fn len(&self) -> usize {
        self._tree.len()
    }

    pub fn stats(&self) -> Option<TreeStats> {
        TreeStatsCollector::collect(&self._tree)
    }
}

/// Given a protein sequence `seq`, traverse a trie of peptides and return a `Vec` of peptide IDs
/// whose sequences are found within the protein's sequence. Peptides are found within the sequence
/// regardless of the location of any cleavage sites; as a result this function provides no
/// guarantee that peptides have any number of enzymatic termini within the given sequence!
///
/// The given slice `seq` will be filtered to remove any newline characters before processing.
fn annotate_sequence<const N: usize>(root: &OpaqueNodePtr<CString, PeptideId, N>, seq: &[u8]) -> Vec<PeptideId> {
    let mut res = Vec::new();

    // Filter the sequence
    let mut filtseq = seq.iter().filter(|&c| !b"\n\r".contains(c)).copied().collect::<Vec<_>>();

    // Iterate backwards to match key storage
    filtseq.reverse();

    for i in 0..filtseq.len() {
        // SAFETY: Since we have an immutable reference to the root node, that
        // means there can only exist other immutable references aside from this one,
        // and no mutable references. That means that no mutating operations can occur
        // on the root node or any child of the root node.
        unsafe {
            _annotate_sequence(&root, &filtseq, i, &mut res);
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

    // The point where we will start checking the key, if we find a leaf node
    let check_start = if was_optimistic {
        pessimistic_depth
    } else {
        new_depth
    };

    // Two possible paths from here -- either a string terminating zero, or the next char
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

    let next = inner_node.lookup_child(seq[start + new_depth])?;

    // Check node type before calling `to_node_ptr`, as converting back
    // to opaque (to return) is expensive!
    match next.node_type() {
        NodeType::Leaf => {
            match next.to_node_ptr() {
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
            Some((next, new_depth + 1, was_optimistic))
        }
    }
}

/// Handle a leaf node as a potential match for the given start. Handles verify any
/// key bytes that were optimistically matched and adding the result to the result vector
/// when a match is confirmed.
///
/// # Safety
///  - No other access or mutation to the `t` Node can happen while this function runs.
unsafe fn handle_leaf<const N: usize>(seq: &[u8], res: &mut Vec<PeptideId>, t: NodePtr<N, LeafNode<CString, PeptideId, N>>, check_start: usize) {
    // SAFETY: The lifetime produced from this is bounded to this scope and does not
    // escape. Further, no other code mutates the node referenced, which is further
    // enforced the "no concurrent reads or writes" requirement on the
    // `_annotate_sequence` function.
    let inner_node = unsafe { t.as_ref() };

    let key = inner_node.key_ref();
    let key_bytes = key.as_bytes();

    if seq.len() < key_bytes.len() {
        // Not enough bytes to match
        return
    }

    if seq[check_start..].starts_with(&key_bytes[check_start..]) {
        res.push(*inner_node.value_ref());
    }
}

#[cfg(test)]
mod tests {
    use blart::TreeMap;
    use crate::{annotate_sequence, PeptideTrie};
    use crate::fasta::annotate_fasta;

    #[test]
    fn test_empty_tree() {
        let tree = PeptideTrie::new();

        let fasta = crate::fasta::Fasta::new(">HEADER\nANYSEQUENCE".as_bytes().into());

        let res = annotate_fasta(
            &fasta,
            tree,
        );

        assert!(res.is_none());

        // If we start parsing the fasta anyway:
        //
        // assert!(res.is_some());
        //
        // let coll_res: Vec<_> = res.unwrap().collect();
        //
        // assert_eq!(coll_res.len(), 1);
        //
        // let prepped_entry = coll_res.iter().next().unwrap();
        //
        // assert!(prepped_entry.peptides().is_empty());
    }

    #[test]
    fn test_singleton_tree() {
        let mut tree = PeptideTrie::new();

        tree.insert("APEPTIDEK".as_bytes());

        let root = TreeMap::into_raw(tree._tree).unwrap();

        let res = annotate_sequence(
            &root,
            "ANYSEQUENCE".as_bytes(),
        );
    }

    #[test]
    fn test_match_at_start() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        tree.insert("APEPTIDER".as_bytes());

        let root = TreeMap::into_raw(tree._tree).unwrap();

        let res = annotate_sequence(
            &root,
            "APEPTIDEKANOTHER".as_bytes(),
        );

        assert!(res.contains(&pep_id));
    }

    #[test]
    fn test_match_at_end() {
        let mut tree = PeptideTrie::new();

        tree.insert("APEPTIDEK".as_bytes());
        let pep_id = tree.insert("ANOTHER".as_bytes());

        let root = TreeMap::into_raw(tree._tree).unwrap();

        let res = annotate_sequence(
            &root,
            "APEPTIDEKANOTHER".as_bytes(),
        );

        assert!(res.contains(&pep_id));
    }
}
