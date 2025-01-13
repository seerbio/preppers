use super::io::slurp_file;
use super::{annotate_sequence, PeptideId, PeptideTrie};
use blart::{OpaqueNodePtr, TreeMap};
use std::ffi::CString;
use std::path::PathBuf;

pub fn read_fasta(fasta_path: PathBuf) -> Fasta {
    let fasta_bytes = slurp_file(fasta_path);

    Fasta {
        file_bytes: fasta_bytes,
    }
}

pub fn annotate_fasta<'a>(fasta: &'a Fasta, peptides: PeptideTrie) -> Option<impl Iterator<Item=PreppedFastaEntry<'a>>> {
    Some(annotate_iter(fasta.iter(), TreeMap::into_raw(peptides._tree)?))
}

fn annotate_iter<'a, T: Iterator<Item=PlainFastaEntry<'a>>, const N: usize>(iter: T, peptides: OpaqueNodePtr<CString, PeptideId, N>) -> impl Iterator<Item=PreppedFastaEntry<'a>> {
    iter.map(
        move |entry| annotate(entry, &peptides)
    )
}

fn annotate<'a, const N: usize>(entry: PlainFastaEntry<'a>, peptides: &OpaqueNodePtr<CString, PeptideId, N>) -> PreppedFastaEntry<'a> {
    let peps = annotate_sequence(peptides, entry.sequence());

    PreppedFastaEntry{
        entry: entry,
        peptides: peps,
    }
}

pub struct Fasta {
    file_bytes: Vec<u8>,
}

impl Fasta {
    pub fn new(file_bytes: Vec<u8>) -> Fasta {
        Fasta {
            file_bytes,
        }
    }
}


impl<'a> Fasta {
    pub fn iter(&'a self) -> FastaIterator<'a> {
        FastaIterator {
            fasta: self,
            byte_index: 0,
        }
    }
}


pub struct FastaIterator<'a> {
    fasta: &'a Fasta,
    byte_index: usize,
}

impl<'a> FastaIterator<'a> {
    fn peek(&'a self) -> &'a u8 {
        &self.fasta.file_bytes[self.byte_index]
    }

    fn eof(&'a self) -> bool {
        self.byte_index >= self.fasta.file_bytes.len()
    }
}

impl<'a> Iterator for FastaIterator<'a> {
    type Item = PlainFastaEntry<'a>;

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
        let h_start = self.byte_index; // Increment to omit ">" from header
        while !self.eof() && !b"\n\r".contains(self.peek()) {
            self.byte_index += 1
        }
        let h_end = self.byte_index;
        let header = &self.fasta.file_bytes[h_start..h_end];

        // Read sequence
        let s_start = self.byte_index + 1;
        while !self.eof() && *self.peek() != b'>' {
            self.byte_index += 1
        }
        let s_end = self.byte_index - 1;
        let sequence = &self.fasta.file_bytes[s_start..s_end];

        Some(
            PlainFastaEntry {
                header,
                sequence,  // TODO: must handle filtering newline bytes!!!
            }
        )
    }
}

pub trait FastaEntry<'a> {
    fn header(&self) -> &'a [u8];
    fn sequence(&self) -> &'a [u8];
}

pub struct PlainFastaEntry<'a> {
    header: &'a [u8],
    sequence: &'a [u8]
}

impl<'a> FastaEntry<'a> for PlainFastaEntry<'a> {
    fn header(&self) -> &'a [u8] {
        &self.header
    }

    fn sequence(&self) -> &'a [u8] {
        &self.sequence
    }
}

pub struct PreppedFastaEntry<'a> {
    entry: PlainFastaEntry<'a>,
    peptides: Vec<PeptideId>,
}

impl PreppedFastaEntry<'_> {
    pub fn peptides(&self) -> &Vec<PeptideId> {
        &self.peptides
    }
}

impl<'a> FastaEntry<'a> for PreppedFastaEntry<'a> {
    fn header(&self) -> &'a [u8] {
        self.entry.header()
    }

    fn sequence(&self) -> &'a [u8] {
        self.entry.sequence()
    }
}

#[cfg(test)]
mod tests {
    use super::{Fasta, FastaEntry};

    #[test]
    fn test_parse_fasta() {
        /// Test parsing a basic fasta returns the correct result
        let fasta = Fasta::new(b">header1\nAAA\n>header2\nBBB\n".to_vec());
        let mut iter = fasta.iter();

        let entry1 = iter.next().unwrap();
        assert_eq!(entry1.header(), b">header1");
        assert_eq!(entry1.sequence(), b"AAA");

        let entry2 = iter.next().unwrap();
        assert_eq!(entry2.header(), b">header2");
        assert_eq!(entry2.sequence(), b"BBB");

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_parse_fasta_end_no_newline() {
        /// Test parsing a basic fasta returns the correct result when the file
        /// is missing a newline at the end
        let fasta = Fasta::new(b">header1\nAAA\n>header2\nBBB".to_vec());
        let mut iter = fasta.iter();

        // We only care about the last entry
        let entry = iter.last().expect("FASTA should not be empty!");
        assert_eq!(entry.header(), b">header2");
        assert_eq!(entry.sequence(), b"BBB");
    }
}
