use super::io::slurp_file;
use super::{annotate_sequence, PeptideId, PeptideTrie};
use blart::{OpaqueNodePtr, TreeMap};
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
    let peps = annotate_sequence(&peptides, seq.as_bytes());

    PreppedFastaEntry{
        header: entry.header().to_owned(),
        sequence: seq,
        peptides: peps,
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
