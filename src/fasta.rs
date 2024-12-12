use super::io::slurp_file;
use super::{annotate_sequence, PeptideId, PeptideTrie};
use blart::{OpaqueNodePtr, TreeMap};
use std::ffi::CString;
use std::path::PathBuf;
use std::time::Instant;

pub fn read_fasta(fasta_path: PathBuf) -> impl Iterator<Item=PlainFastaEntry> {
    let read_start = Instant::now();

    let fasta_bytes = slurp_file(fasta_path);

    let read_duration = read_start.elapsed();
    println!("Read FASTA in {:.4} sec", read_duration.as_secs_f64());

    FastaIterator {
        file_bytes: fasta_bytes,
        byte_index: 0,
    }
}

pub fn annotate_fasta(fasta_path: PathBuf, peptides: PeptideTrie) -> impl Iterator<Item=PreppedFastaEntry> {
    annotate_iter(read_fasta(fasta_path), TreeMap::into_raw(peptides._tree).unwrap())
}

fn annotate_iter<T: Iterator<Item=PlainFastaEntry>, const N: usize>(iter: T, peptides: OpaqueNodePtr<CString, PeptideId, N>) -> impl Iterator<Item=PreppedFastaEntry> {
    iter.map(
        move |entry| annotate(entry, &peptides)
    )
}

fn annotate<const N: usize>(entry: PlainFastaEntry, peptides: &OpaqueNodePtr<CString, PeptideId, N>) -> PreppedFastaEntry {
    let peps = annotate_sequence(peptides, entry.sequence());

    PreppedFastaEntry{
        entry: entry,
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
                // TODO: implicit copy!
                header: header.to_vec(),
                sequence: sequence
                    .to_vec()
                    .into_iter()
                    // TODO: filtering and collecting is an implicit copy!
                    .filter(
                        |b| !b"\r\n".contains(b)
                    )
                    .collect()
            }
        )
    }
}

pub trait FastaEntry {
    fn header(&self) -> &[u8];
    fn sequence(&self) -> &[u8];
}

pub struct PlainFastaEntry {
    header: Vec<u8>,
    sequence: Vec<u8>,
}

impl FastaEntry for PlainFastaEntry {
    fn header(&self) -> &[u8] {
        &self.header
    }

    fn sequence(&self) -> &[u8] {
        &self.sequence
    }
}

pub struct PreppedFastaEntry {
    entry: PlainFastaEntry,
    peptides: Vec<PeptideId>,
}

impl PreppedFastaEntry {
    pub fn peptides(&self) -> &Vec<PeptideId> {
        &self.peptides
    }
}

impl FastaEntry for PreppedFastaEntry {
    fn header(&self) -> &[u8] {
        self.entry.header()
    }

    fn sequence(&self) -> &[u8] {
        self.entry.sequence()
    }
}
