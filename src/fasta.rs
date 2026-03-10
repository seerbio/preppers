use super::io::slurp_file;
use super::{annotate_sequence, termini, PeptideId, PeptideTrie};
use blart::{OpaqueNodePtr, TreeMap};
use std::ffi::CString;
use std::path::PathBuf;
use std::slice::Iter;
use fancy_regex::Regex;

pub fn read_fasta(fasta_path: PathBuf) -> Fasta {
    let fasta_bytes = slurp_file(fasta_path);

    Fasta {
        file_bytes: fasta_bytes,
    }
}

/// Annotate the fasta using peptides from the given trie, returning all substring matches.
/// Returns `None` if the trie is empty.
pub fn annotate_fasta(fasta: &Fasta, peptides: PeptideTrie) -> Option<impl Iterator<Item=PreppedFastaEntry<'_>>> {
    Some(annotate_iter(fasta.iter(), TreeMap::into_raw(peptides._tree)?))
}

/// Annotate the fasta using peptides from the given trie, filtering results using the given pattern
/// string. Returns `None` if the trie is empty.
pub fn annotate_fasta_filtered(fasta: &Fasta, peptides: PeptideTrie, enzyme_patt: Regex, n_req_termini: u8) -> Option<impl Iterator<Item=PreppedFastaEntry<'_>>> {
    let iter = match annotate_fasta(fasta, peptides) {
        Some(iter) => iter,
        None => return None,
    };

    Some(termini::filter_entry_termini(iter, enzyme_patt, n_req_termini))
}

fn annotate_iter<'a, T: Iterator<Item=PlainFastaEntry<'a>>, const N: usize>(iter: T, peptides: OpaqueNodePtr<CString, PeptideId, N>) -> impl Iterator<Item=PreppedFastaEntry<'a>> {
    iter.map(
        move |entry| annotate(entry, &peptides)
    )
}

fn annotate<'a, const N: usize>(entry: PlainFastaEntry<'a>, peptides: &OpaqueNodePtr<CString, PeptideId, N>) -> PreppedFastaEntry<'a> {
    let (seq, idxs) = annotate_sequence(peptides, entry.sequence());

    PreppedFastaEntry{
        entry: entry,
        sequence: seq,
        peptide_indices: idxs,
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
    /// Iterate over the entries in the FASTA, returning entries
    /// as slices into the FASTA's contents; the returned sequences
    /// thus may contain newline characters.
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
        let h_start = self.byte_index;  // This will include the header's '>'

        // Read until the next newline (or EOF)
        while !self.eof() && !b"\n\r".contains(self.peek()) {
            self.byte_index += 1
        }
        let h_end = self.byte_index;
        let header = &self.fasta.file_bytes[h_start..h_end];

        // Read sequence
        let s_start = self.byte_index + 1;

        // Read until the next header (or EOF); don't worry about newlines
        while !self.eof() && *self.peek() != b'>' {
            self.byte_index += 1
        }
        let s_end = self.byte_index;

        let sequence = &self.fasta.file_bytes[s_start..s_end];

        Some(
            PlainFastaEntry {
                header,
                sequence
            }
        )
    }
}

pub trait FastaEntry {
    fn header(&self) -> &[u8];
    fn sequence(&self) -> &[u8];
}

#[derive(Debug)]
pub struct PlainFastaEntry<'a> {
    header: &'a [u8],
    sequence: &'a [u8]
}

impl<'a> FastaEntry for PlainFastaEntry<'a> {
    fn header(&self) -> &'a [u8] {
        &self.header
    }

    fn sequence(&self) -> &'a [u8] {
        &self.sequence
    }
}

/// A peptide match represented as `(peptide_id, start, stop)`.
///
/// Index semantics are **closed**: `start` is the index of the first residue in
/// the peptide match, and `stop` is the index of the last residue in the peptide
/// match (inclusive), not one-past-the-end.
pub type PeptideHit = (PeptideId, usize, usize);

#[derive(Debug)]
pub struct PreppedFastaEntry<'a> {
    /// Zero-copy reference to the original FASTA entry
    pub (crate) entry: PlainFastaEntry<'a>,

    /// Owned copy of the sequence, post-normalization
    pub (crate) sequence: Vec<u8>,

    pub (crate) peptide_indices: Vec<PeptideHit>,
}

impl PreppedFastaEntry<'_> {
    pub fn peptides(&self) -> impl ExactSizeIterator<Item=&PeptideId> {
        self.peptide_indices.iter().map(|(id, _, _)| id)
    }

    pub fn peptide_indices(&self) -> Iter<'_, PeptideHit> {
        self.peptide_indices.iter()
    }
}

impl<'a> IntoIterator for &'a PreppedFastaEntry<'_> {
    type Item = &'a PeptideHit;
    type IntoIter = Iter<'a, PeptideHit>;

    fn into_iter(self) -> Self::IntoIter {
        self.peptide_indices()
    }
}

impl FastaEntry for PreppedFastaEntry<'_> {
    fn header(&self) -> &[u8] {
        self.entry.header()
    }

    fn sequence(&self) -> &[u8] {
        self.sequence.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::{annotate_fasta, annotate_fasta_filtered, Fasta, FastaEntry};
    use crate::PeptideTrie;
    use blart::AsBytes;
    use fancy_regex::Regex;

    /// Test parsing a basic fasta returns the correct result
    #[test]
    fn test_parse_fasta() {
        let fasta = Fasta::new(b">header1\nAAA\nAAA\n>header2\nBBBBBB\n".to_vec());
        let mut iter = fasta.iter();

        let entry1 = iter.next().unwrap();
        assert_eq!(entry1.header(), b">header1");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry1.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"AAAAAA");

        let entry2 = iter.next().unwrap();
        assert_eq!(entry2.header(), b">header2");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry2.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"BBBBBB");

        assert!(iter.next().is_none());
    }

    /// Test parsing a basic fasta returns the correct result when the file
    /// is missing a newline at the end
    #[test]
    fn test_parse_fasta_no_end_newline() {
        let fasta = Fasta::new(b">header1\nAAA\n>header2\nBBB".to_vec());
        let iter = fasta.iter();

        // We only care about the last entry
        let entry = iter.last().expect("FASTA should not be empty!");
        assert_eq!(entry.header(), b">header2");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"BBB");
    }

    /// Test parsing a fasta with windows line endings returns the correct result
    #[test]
    fn test_parse_fasta_windows() {
        let fasta = Fasta::new(b">header1\r\nAAA\r\nAAA\r\n>header2\r\nBBBBBB\r\n".to_vec());
        let mut iter = fasta.iter();

        let entry1 = iter.next().unwrap();
        assert_eq!(entry1.header(), b">header1");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry1.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"AAAAAA");

        let entry2 = iter.next().unwrap();
        assert_eq!(entry2.header(), b">header2");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry2.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"BBBBBB");

        assert!(iter.next().is_none());
    }

    /// Test parsing a fasta with widnows line endings returns the correct result
    /// when the file is missing a newline at the end
    #[test]
    fn test_parse_fasta_no_end_newline_windows() {
        let fasta = Fasta::new(b">header1\r\nAAA\r\n>header2\r\nBBB".to_vec());
        let iter = fasta.iter();

        // We only care about the last entry
        let entry = iter.last().expect("FASTA should not be empty!");
        assert_eq!(entry.header(), b">header2");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"BBB");
    }

    /// Test parsing a fasta with mixed line endings returns the correct result
    #[test]
    fn test_parse_fasta_mixed() {
        let fasta = Fasta::new(b">header1\nAAA\nAAA\n>header2\r\nBBBBBB\r\n>header3\rCCCCCC\r".to_vec());
        let mut iter = fasta.iter();

        let entry1 = iter.next().unwrap();
        assert_eq!(entry1.header(), b">header1");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry1.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"AAAAAA");

        let entry2 = iter.next().unwrap();
        assert_eq!(entry2.header(), b">header2");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry2.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"BBBBBB");

        let entry3 = iter.next().unwrap();
        assert_eq!(entry3.header(), b">header3");

        // ignore newline characters in sequence; these are stripped downstream
        assert_eq!(entry3.sequence().iter().filter(|b| !b"\r\n".contains(b)).copied().collect::<Vec<_>>().as_bytes(), b"CCCCCC");

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_empty_tree() {
        let tree = PeptideTrie::new();

        let fasta = Fasta::new(">HEADER\nANYSEQUENCE".as_bytes().into());

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

        let pep_id = tree.insert("APEPTIDEK".as_bytes());

        let fasta = Fasta::new(">HEADER\nAPEPTIDEKANOTHER".as_bytes().into());

        let res = annotate_fasta(
            &fasta,
            tree,
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.sequence(), "APEPTIDEKANOTHER".as_bytes());
        assert_eq!(prepped_entry.peptides().len(), 1);
        assert_eq!(*prepped_entry.peptides().next().unwrap(), pep_id);
        assert_eq!(*prepped_entry.peptide_indices().next().unwrap(), (pep_id, 0, 8));
    }

    #[test]
    fn test_match_post_normalization_0() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        tree.insert("APEPTIDER".as_bytes());

        let fasta = Fasta::new(">HEADER\nAPEPTIDEK\nANOTHER".as_bytes().into());

        let res = annotate_fasta(
            &fasta,
            tree,
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 1);
        assert_eq!(*prepped_entry.peptides().next().unwrap(), pep_id);
        assert_eq!(*prepped_entry.peptide_indices().next().unwrap(), (pep_id, 0, 8));
    }

    #[test]
    fn test_match_post_normalization_1() {
        let mut tree = PeptideTrie::new();

        tree.insert("APEPTIDER".as_bytes());
        let pep_id = tree.insert("ANOTHER".as_bytes());

        let fasta = Fasta::new(">HEADER\nAPEPTIDEK\nANOTHER".as_bytes().into());

        let res = annotate_fasta(
            &fasta,
            tree,
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 1);
        assert_eq!(*prepped_entry.peptides().next().unwrap(), pep_id);
        assert_eq!(*prepped_entry.peptide_indices().next().unwrap(), (pep_id, 9, 15));
    }

    #[test]
    fn test_match_post_normalization_2() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        tree.insert("APEPTIDER".as_bytes());

        let fasta = Fasta::new(">HEADER\nAPEPT\nIDEK\nANOTHER".as_bytes().into());

        let res = annotate_fasta(
            &fasta,
            tree,
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 1);
        assert_eq!(*prepped_entry.peptides().next().unwrap(), pep_id);
        assert_eq!(*prepped_entry.peptide_indices().next().unwrap(), (pep_id, 0, 8));
    }

    #[test]
    fn test_match_post_normalization_3() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        tree.insert("APEPTIDER".as_bytes());

        let fasta = Fasta::new(">HEADER\r\nAPEPT\r\nIDEK\r\nANOTHER".as_bytes().into());

        let res = annotate_fasta(
            &fasta,
            tree,
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 1);
        assert_eq!(*prepped_entry.peptides().next().unwrap(), pep_id);
        assert_eq!(*prepped_entry.peptide_indices().next().unwrap(), (pep_id, 0, 8));
    }

    #[test]
    fn test_annotate_filtered() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        tree.insert("PEPTIDER".as_bytes());

        let fasta = Fasta::new(">HEADER\nPEPTIDEKPEPTIDER\nAPEPT\nIDEKANOTHER".as_bytes().into());

        let res = annotate_fasta_filtered(
            &fasta,
            tree,
            Regex::new("(?<=[KR])(?!P)").unwrap(),
            2
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 1);
        assert_eq!(*prepped_entry.peptide_indices().next().unwrap(), (pep_id, 16, 24));
    }

    #[test]
    fn test_annotate_filtered_semi() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        let pep_id_2 = tree.insert("PEPTIDER".as_bytes());

        let fasta = Fasta::new(">HEADER\nPEPTIDEKPEPTIDER\nAPEPT\nIDEKANOTHER".as_bytes().into());

        let res = annotate_fasta_filtered(
            &fasta,
            tree,
            Regex::new("(?<=[KR])(?!P)").unwrap(),
            1
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 2);
        assert!(prepped_entry.peptide_indices().any(|h| *h == (pep_id, 16, 24)));
        assert!(prepped_entry.peptide_indices().any(|h| *h == (pep_id_2, 8, 15)));
    }

    /// Test that we can use a carefully-crafted regex to permit treating n-term met. excision as an allowed
    /// terminus without having to implement special-case logic.
    #[test]
    fn test_annotate_filtered_nterm_m() {
        let mut tree = PeptideTrie::new();

        let pep_id = tree.insert("APEPTIDEK".as_bytes());
        let pep_id_2 = tree.insert("PEPTIDER".as_bytes());

        let fasta = Fasta::new(">HEADER\nMPEPTIDER\nAPEPT\nIDEKANOTHER".as_bytes().into());

        let res = annotate_fasta_filtered(
            &fasta,
            tree,
            Regex::new("(?:(?<=[KR])(?!P))|(?<=^M)").unwrap(),
            2
        );

        assert!(res.is_some());

        let coll_res: Vec<_> = res.unwrap().collect();

        println!("Filtered peptides: {:?}", coll_res);

        assert_eq!(coll_res.len(), 1);

        let prepped_entry = coll_res.iter().next().unwrap();

        assert_eq!(prepped_entry.peptides().len(), 2);
        assert!(prepped_entry.peptide_indices().any(|h| *h == (pep_id, 9, 17)));
        assert!(prepped_entry.peptide_indices().any(|h| *h == (pep_id_2, 1, 8)));
    }
}
