use std::collections::BTreeSet;
use std::ffi::CString;
use std::path::PathBuf;
use super::{PeptideId, PeptideTrie};
use super::io::slurp_file;

pub fn read_fasta(fasta_path: PathBuf) -> impl Iterator<Item: FastaEntry> {
    FastaIterator {
        file_bytes: slurp_file(fasta_path),
        byte_index: 0,
    }
}

pub fn annotate_fasta<'a>(fasta_path: PathBuf, peptides: &'a PeptideTrie) -> impl Iterator<Item=PreppedFastaEntry> + use <'a> {
    annotate_iter(read_fasta(fasta_path), peptides)
}

fn annotate_iter<'a, T: Iterator<Item: FastaEntry>>(iter: T, peptides: &'a PeptideTrie) -> impl Iterator<Item=PreppedFastaEntry> + use <'a, T> {
    iter.map(
        |entry| annotate(&entry, peptides)
    )
}

fn annotate<'a>(entry: &impl FastaEntry, peptides: &'a PeptideTrie) -> PreppedFastaEntry {
    let seq = entry.sequence().to_owned();
    let peps = get_peptides_for_sequence(&seq, peptides);
    PreppedFastaEntry{
        header: entry.header().to_owned(),
        sequence: seq,
        peptides: peps,
    }
}

const MIN_PFX: usize = 5;

fn get_peptides_for_sequence(seq: &String, peptides: &PeptideTrie) -> Vec<PeptideId> {
    let (first_key, _) = peptides._tree.first_key_value().unwrap();
    let (last_key, _) = peptides._tree.last_key_value().unwrap();

    let mut state = BTreeSet::<usize>::new();
    let mut res = Vec::new();

    for i in 0..seq.len() - MIN_PFX {
        let pfx = &seq[i..i + MIN_PFX];

        // let res: Vec<_> = map.prefix(pfx.as_bytes()).collect();
        // println!("{i}: {pfx} -- {res:?}");

        // Due to the implementation of prefix() we must check
        // the prefix is within the range of values in the tree;
        // if it's outside we instead get the full set of keys!!
        let is_inrange = pfx.as_bytes() >= first_key.as_bytes() && pfx.as_bytes() <= last_key.as_bytes();

        // let peps: Vec<_> = match is_inrange {
        //     true => peptides._tree.prefix(pfx.as_bytes())
        //         // This filter is required if the prefix is outside the tree's bounds?
        //         //.filter(|(k, _)| k.to_bytes()[0..MIN_PFX] == *pfx.as_bytes())
        //         .collect(),
        //     _ => Vec::new(),
        // };
        // println!("{i}: {pfx} {is_inrange} -- {peps:?}");
        //
        // if peps.iter().any(|(k, _)| k.to_bytes()[0..MIN_PFX] != *pfx.as_bytes()) {
        //     panic!("Some returned values had wrong prefix! Prefix: {pfx}; Result: {peps:?}")
        // }

        if is_inrange && peptides._tree.prefix(pfx.as_bytes()).next().is_some() {
            // This prefix has at least one peptide, so we will keep checking it
            state.insert(i);
        }

        // Check for peptides ending at this index
        let mut to_rm = Vec::new();
        for start in &state {
            let putseq = &seq[*start..i];

            let prefix = peptides._tree.prefix(putseq.as_bytes());
            let mut n = 0;
            for (pep, id) in prefix {
                n += 1;

                if pep.as_bytes() == putseq.as_bytes() {
                    // We found a peptide; add it to the result
                    res.push(*id);
                }
            }
            if n == 0 {
                to_rm.push(*start);
            }
        }

        for start in to_rm {
            state.remove(&start);
        }
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