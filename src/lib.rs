pub mod io;
pub mod fasta;

// Reexports
pub use fasta::read_fasta;

// Imports
use std::ffi::CString;
use blart::TreeMap;
use blart::visitor::{TreeStats, TreeStatsCollector};

// Types
type PeptideId = u64;

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

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn it_works() {
    //     let result = add(2, 2);
    //     assert_eq!(result, 4);
    // }
}
