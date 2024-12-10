use std::collections::BTreeSet;
use std::ffi::CString;
use std::path::PathBuf;
use blart::{InnerNode, InnerNode16, InnerNode256, InnerNode4, InnerNode48, LeafNode, NodePtr, NodeType, OpaqueNodePtr, TreeMap};
use blart::visitor::{Visitable, Visitor};
use super::{PeptideId, PeptideTrie};
use super::io::slurp_file;

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

fn annotate<'a, const N: usize>(entry: &impl FastaEntry, peptides: &'a OpaqueNodePtr<CString, PeptideId, N>) -> PreppedFastaEntry {
    let seq = entry.sequence().to_owned();
    let peps = get_peptides_for_sequence(&seq, &peptides);

    PreppedFastaEntry{
        header: entry.header().to_owned(),
        sequence: seq,
        peptides: peps,
    }
}

fn get_peptides_for_sequence<const N: usize>(seq: &String, peptides: &OpaqueNodePtr<CString, PeptideId, N>) -> Vec<PeptideId> {
    let mut res = Vec::new();

    for i in 0..seq.len() {
        peptides.visit_with(&mut ProteinSequenceVisitor::new(seq.as_bytes(), i, &mut res));
    }

    res
}

struct ProteinSequenceVisitor<'a> {
    seq: &'a [u8],
    start: usize,
    res: &'a mut Vec<PeptideId>,
    idx: usize,
}

impl<'a> ProteinSequenceVisitor<'a> {
    fn new(seq: &'a [u8], start: usize, res: &'a mut Vec<PeptideId>) -> ProteinSequenceVisitor<'a> {
        ProteinSequenceVisitor { seq, start, res, idx: 0, }
    }
}

impl<const N: usize> Visitor<CString, PeptideId, N> for ProteinSequenceVisitor<'_> {
    type Output=();

    fn default_output(&self) -> Self::Output { }

    fn combine_output(&self, o1: Self::Output, o2: Self::Output) -> Self::Output { }

    fn visit_node4(&mut self, t: &InnerNode4<CString, PeptideId, N>) -> Self::Output {
        if self.start + self.idx >= self.seq.len() {
            return ()
        }

        let next = t.lookup_child(self.seq[self.idx]);

        next.map(|n|
            n.visit_with(
                &mut ProteinSequenceVisitor {
                    seq: self.seq,
                    start: self.start,
                    res: self.res,
                    idx: self.idx + 1,
                }
            )
        );
    }

    fn visit_node16(&mut self, t: &InnerNode16<CString, PeptideId, N>) -> Self::Output {
        if self.start + self.idx >= self.seq.len() {
            return ()
        }

        let next = t.lookup_child(self.seq[self.idx]);

        next.map(|n|
            n.visit_with(
                &mut ProteinSequenceVisitor {
                    seq: self.seq,
                    start: self.start,
                    res: self.res,
                    idx: self.idx + 1,
                }
            )
        );
    }

    fn visit_node48(&mut self, t: &InnerNode48<CString, PeptideId, N>) -> Self::Output {
        if self.start + self.idx >= self.seq.len() {
            return ()
        }

        let next = t.lookup_child(self.seq[self.idx]);

        next.map(|n|
            n.visit_with(
                &mut ProteinSequenceVisitor {
                    seq: self.seq,
                    start: self.start,
                    res: self.res,
                    idx: self.idx + 1,
                }
            )
        );
    }

    fn visit_node256(&mut self, t: &InnerNode256<CString, PeptideId, N>) -> Self::Output {
        if self.start + self.idx >= self.seq.len() {
            return ()
        }

        let next = t.lookup_child(self.seq[self.idx]);

        next.map(|n|
            n.visit_with(
                &mut ProteinSequenceVisitor {
                    seq: self.seq,
                    start: self.start,
                    res: self.res,
                    idx: self.idx + 1,
                }
            )
        );

    }

    fn visit_leaf(&mut self, t: &LeafNode<CString, PeptideId, N>) -> Self::Output {
        let seq = &self.seq[self.start..self.start + self.idx - 1];

        if t.matches_full_key(seq) {
            self.res.push(*t.value_ref());
        } else {
            println!("Leaf mismatch! Reached {:?} with {:?}", t.key_ref(), String::from_utf8(seq.to_vec()).unwrap());
        }
    }
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