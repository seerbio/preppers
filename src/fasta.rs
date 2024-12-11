use super::io::slurp_file;
use super::{PeptideId, PeptideTrie};
use blart::{AsBytes, AttemptOptimisticPrefixMatch, ConcreteNodePtr, InnerNode, NodePtr, OpaqueNodePtr, PessimisticMismatch, TreeMap};
use std::ffi::CString;
use std::path::PathBuf;

pub fn read_fasta(fasta_path: PathBuf) -> impl Iterator<Item: FastaEntry> {
    FastaIterator {
        file_bytes: slurp_file(fasta_path),
        byte_index: 0,
    }
}

pub fn annotate_fasta(fasta_path: PathBuf, peptides: PeptideTrie) -> impl Iterator<Item=PreppedFastaEntry> {
    annotate_iter(read_fasta(fasta_path), TreeMap::into_raw(peptides._tree).unwrap())
}

fn annotate_iter<T: Iterator<Item: FastaEntry>, const N: usize>(iter: T, peptides: OpaqueNodePtr<CString, PeptideId, N>) -> impl Iterator<Item=PreppedFastaEntry> {
    iter.map(
        move |entry| annotate(&entry, &peptides)
    )
}

fn annotate<const N: usize>(entry: &impl FastaEntry, peptides: &OpaqueNodePtr<CString, PeptideId, N>) -> PreppedFastaEntry {
    let seq = entry.sequence().to_owned();
    let peps = get_peptides_for_sequence(&peptides, &seq);

    PreppedFastaEntry{
        header: entry.header().to_owned(),
        sequence: seq,
        peptides: peps,
    }
}

fn get_peptides_for_sequence<const N: usize>(peptides: &OpaqueNodePtr<CString, PeptideId, N>, seq: &String) -> Vec<PeptideId> {
    let mut res = Vec::new();

    for i in 0..seq.len() {
        match_peptides(peptides, &seq, i, 0, &mut res);
    }

    res
}

struct FastaIterator {
    file_bytes: Vec<u8>,
    byte_index: usize,
}

impl FastaIterator {
    fn peek(&self) -> &u8 {
        &self.file_bytes[self.byte_index]
    }

    fn eof(&self) -> bool {
        self.byte_index >= self.file_bytes.len()
    }
}

impl Iterator for FastaIterator {
    type Item = PlainFastaEntry;

    fn next(&mut self) -> Option<Self::Item> {
        if self.eof() {
            return None
        }

        while !self.eof() && b"\n\r".contains(self.peek()) {
            self.byte_index += 1
        }

        // Read header
        if !self.eof() && *self.peek() != b'>' {
            panic!("Did not find FASTA header at index {}", self.byte_index)
        }
        let h_start = self.byte_index;
        while !self.eof() && !b"\n\r".contains(self.peek()) {
            self.byte_index += 1
        }
        let h_end = self.byte_index;
        let header = &self.file_bytes[h_start..h_end];

        // Read sequence
        let s_start = self.byte_index;
        while !self.eof() && *self.peek() != b'>' {
            self.byte_index += 1
        }
        let s_end = self.byte_index;
        let sequence = &self.file_bytes[s_start..s_end];

        Some(
            PlainFastaEntry {
                header: String::from_utf8(header.to_vec()).expect("Invalid UTF8 in header!"),
                sequence: String::from_utf8(
                        sequence.to_vec()
                            .into_iter()
                            .filter(
                                |b| !b"\r\n".contains(b)
                            )
                            .collect::<Vec<_>>()
                    ).expect("Invalid UTF8 in sequence!"),
            }
        )
    }
}

pub trait FastaEntry {
    fn header(&self) -> &String;
    fn sequence(&self) -> &String;
}

pub struct PlainFastaEntry {
    header: String,
    sequence: String,
}

impl FastaEntry for PlainFastaEntry {
    fn header(&self) -> &String {
        &self.header
    }

    fn sequence(&self) -> &String {
        &self.sequence
    }
}

pub struct PreppedFastaEntry {
    header: String,
    sequence: String,
    peptides: Vec<PeptideId>,
}

impl PreppedFastaEntry {
    pub fn peptides(&self) -> &Vec<PeptideId> {
        &self.peptides
    }
}

impl FastaEntry for PreppedFastaEntry {
    fn header(&self) -> &String {
        &self.header
    }

    fn sequence(&self) -> &String {
        &self.sequence
    }
}

/// Search in the given tree for any peptides that are substrings of `seq` starting at
/// `start`, having already traversed `depth` characters within the tree.
///
/// This is essentially equivalent to `search_unchecked` but permits recording multiple
/// matches while advancing through the sequence. However, due to implementation details
/// within `blart` we skip the optimization of using optimistic/pessimistic matching and
/// simply compare the full key at each leaf we find (the pessimisitic approach).
///
/// The challenge here is that path compression of the tree means that a single node
/// might represent many bytes in the key, and those bytes might be implicit. This
/// means that we don't know what length of key to check!
fn match_peptides<const N: usize>(node: &OpaqueNodePtr<CString, PeptideId, { N }>, seq: &String, start: usize, depth: usize, res: &mut Vec<PeptideId>) {
    match node.to_node_ptr() {
        ConcreteNodePtr::Node4(t) => {
            let inner_node = t.read();

            let prefix_match = inner_node.attempt_pessimistic_match_prefix(&seq.as_bytes()[start + depth..]);

            if prefix_match.is_err() {
                println!("failed to match at depth {} ({})", depth, &seq.as_str()[start + depth..]);
                return
            }

            let m = prefix_match.unwrap();

            let new_depth = depth + m.matched_bytes + 1;

            if new_depth >= seq.len() {
                println!("too close to end at depth {}", depth);
                return
            }

            println!("matched {} bytes; searching for next byte {} at depth {}", m.matched_bytes, &seq.as_str()[start + new_depth..start + new_depth + 1], new_depth);

            let next = inner_node.lookup_child(seq.as_bytes()[start + new_depth]);

            next.as_ref().map(|n|
                match_peptides(
                    n,
                    seq,
                    start,
                    new_depth + 1, // we matched one additional byte with `lookup_child`
                    res
                )
            );
        }
        ConcreteNodePtr::Node16(t) => {
            todo!()
        }
        ConcreteNodePtr::Node48(t) => {
            todo!()
        }
        ConcreteNodePtr::Node256(t) => {
            todo!()
        }
        ConcreteNodePtr::LeafNode(t) => {
            // We might reach a leaf early, if it has a unique prefix shorter than the key.
            // Here we just check the length.
            let n = t.read();
            let key = n.key_ref();
            println!("checking leaf {:?}", key);
            let len = key.as_bytes().len();
            if start+len < seq.len()-1 && key.as_bytes().eq(&seq.as_bytes()[start..start+len]) {
                res.push(*t.read().value_ref());
            }
        }
    }
}
