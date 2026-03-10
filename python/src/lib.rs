use std::path::PathBuf;
use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::wrap_pyfunction;

use ::preppers::fasta::{read_fasta, Fasta, FastaEntry, PreppedFastaEntry};
use ::preppers::{PeptideId, PeptideTrie};

/// Annotates a FASTA file
///
/// Parameters
/// ----------
/// peptides: [str]
///     A list of peptide sequence strings
/// input: str
///     The path to the FASTA file
#[pyfunction]
fn annotate_fasta<'a>(peptides: &Bound<'a, PyList>, input: &Bound<'a, PyAny>) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let fasta = read_fasta(PathBuf::from(input.str()?.to_str()?));

    _annotate_fasta(peptides, fasta)
}

/// Annotates a FASTA file, given as a byte array
///
/// Parameters
/// ----------
/// peptides: [str]
///     A list of peptide sequence strings
/// input: str
///     The path to the FASTA file
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

    let opt_iter = ::preppers::fasta::annotate_fasta(&fasta, trie);

    if opt_iter.is_none() {
        return Err(pyo3::exceptions::PyValueError::new_err("Unable to annotate FASTA file"));
    }

    let iter = opt_iter.unwrap();

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
            peptides: entry.peptides().into_iter().map(|p| p.to_owned()).collect(),
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
    fn peptides(&self) -> &Vec<PeptideId> {
        &self.peptides
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self)
    }
}

/// Python wrapper around the "Preppers" Rust library
#[pymodule(name="preppers")]
fn preppers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(annotate_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(annotate_fasta_bytes, m)?)
}
