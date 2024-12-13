use std::collections::BTreeMap;
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyDict, PyIterator, PyList, PyString, PyTuple};
use pyo3::wrap_pyfunction;
use preppers::{PeptideTrie, PeptideId};
use preppers::fasta::{read_fasta, Fasta, FastaEntry, PreppedFastaEntry};

/// Annotates a FASTA file.
#[pyfunction]
fn annotate_fasta<'a>(peptides: &Bound<'a, PyList>, input: &str) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let fasta = read_fasta(input.into());

    _annotate_fasta(peptides, fasta)
}

/// Annotates a FASTA file.
#[pyfunction]
fn annotate_fasta_bytes(peptides: &Bound<PyList>, input: &[u8]) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let fasta = Fasta::new(input.to_vec());

    _annotate_fasta(peptides, fasta)
}

/// Annotates a FASTA file.
fn _annotate_fasta(peptides: &Bound<PyList>, fasta: Fasta) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let mut trie = PeptideTrie::new();

    let mut peptide_ids = Vec::new();

    for peptide in peptides.iter() {
        let sequence = peptide.str()?;

        let id = trie.insert(sequence.to_str()?.as_bytes());

        peptide_ids.push((sequence.to_string(), id));
    }

    let iter = preppers::fasta::annotate_fasta(&fasta, trie);

    Ok(
        (
            peptide_ids,
            iter.map(|e| PreppedFastaEntryCopy::from(e)).collect()
        )
    )
}

#[pyclass]
#[derive(Debug)]
struct PreppedFastaEntryCopy {
    header: String,
    sequence: String,
    peptides: Vec<PeptideId>,
}

impl From<PreppedFastaEntry<'_>> for PreppedFastaEntryCopy {
    fn from(entry: PreppedFastaEntry) -> Self {
        PreppedFastaEntryCopy {
            header: String::from_utf8(entry.header().to_vec()).unwrap(),
            sequence: String::from_utf8(entry.sequence().to_vec()).unwrap(),
            peptides: entry.peptides().iter().map(|p| p.to_owned()).collect(),
        }
    }
}

#[pymethods]
impl PreppedFastaEntryCopy {
    #[getter]
    fn header(&self) -> &str {
        &self.header
    }

    #[getter]
    fn sequence(&self) -> &str {
        &self.sequence
    }

    #[getter]
    fn peptides(&self) -> Vec<PeptideId> {
        self.peptides.clone()
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self)
    }
}

/// This module is implemented in Rust.
#[pymodule]
fn preppyrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(annotate_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(annotate_fasta_bytes, m)?)
}