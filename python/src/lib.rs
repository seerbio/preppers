use std::iter;
use pyo3::prelude::*;
use pyo3::types::{PyInt, PyList, PyString};
use pyo3::wrap_pyfunction;
use std::path::PathBuf;

use fancy_regex::Regex;

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
/// enzyme_patt: str
///     A RegEx string that will give a zero-width match at enzymatic cleavages (optional).
///
///     Examples:
///     - Strict Trypsin: ``r"(?<=[KR])(?!P)"``"
///     - Trypsin(/P): ``r"(?<=[KR])"``
///     - Strict Trypsin, allowing N-term met. excision: ``r"(?<=[KR])(?!P)|(?<=^M)"``
/// require_termini: int
///     An integer in ``[0, 2]``, defining the required number of enzymatic matches for a peptide-
///     protein connection to be reported. Ignored when ``enzyme_patt`` is ``None``. When ``0``,
///     behavior is identical to ``enzyme_patt = None``, and ``enzyme_patt`` is ignored.
///     Default: 2
#[pyfunction(signature = (peptides, input, enzyme_patt=None, require_termini=None))]
fn annotate_fasta<'a>(peptides: &Bound<'a, PyList>, input: &Bound<'a, PyAny>, enzyme_patt: Option<&Bound<'a, PyString>>, require_termini: Option<&Bound<'a, PyInt>>) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let fasta = read_fasta(PathBuf::from(input.str()?.to_str()?));

    _annotate_fasta(peptides, fasta, enzyme_patt, require_termini)
}

/// Annotates a FASTA file, given as a byte array
///
/// Parameters
/// ----------
/// peptides: [str]
///     A list of peptide sequence strings
/// input: str
///     The path to the FASTA file
/// enzyme_patt: str
///     A RegEx string that will give a zero-width match at enzymatic cleavages (optional).
///
///     Examples:
///     - Strict Trypsin: ``r"(?<=[KR])(?!P)"``"
///     - Trypsin(/P): ``r"(?<=[KR])"``
///     - Strict Trypsin, allowing N-term met. excision: ``r"(?<=[KR])(?!P)|(?<=^M)"``
/// require_termini: int
///     An integer in ``[0, 2]``, defining the required number of enzymatic matches for a peptide-
///     protein connection to be reported. Ignored when ``enzyme_patt`` is ``None``. When ``0``,
///     behavior is identical to ``enzyme_patt = None``, and ``enzyme_patt`` is ignored.
///     Default: 2
#[pyfunction(signature = (peptides, input, enzyme_patt=None, require_termini=None))]
fn annotate_fasta_bytes(peptides: &Bound<PyList>, input: &[u8], enzyme_patt: Option<&Bound<PyString>>, require_termini: Option<&Bound<PyInt>>) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let fasta = Fasta::new(input.to_vec());

    _annotate_fasta(peptides, fasta, enzyme_patt, require_termini)
}

/// Annotates a FASTA file.
fn _annotate_fasta(peptides: &Bound<PyList>, fasta: Fasta, enzyme_patt: Option<&Bound<PyString>>, require_termini: Option<&Bound<PyInt>>) -> PyResult<(Vec<(String, PeptideId)>, Vec<PreppedFastaEntryCopy>)> {
    let n_req_termini = match require_termini {
        Some(req) => req.extract()?,
        None => 2,
    };

    if enzyme_patt.is_some() && (n_req_termini < 0 || n_req_termini > 2) {
        return Err(pyo3::exceptions::PyValueError::new_err("require_termini must be in [0, 2]"));
    }

    let mut trie = PeptideTrie::new();

    let mut peptide_ids = Vec::new();

    for peptide in peptides.iter() {
        let sequence = peptide.str()?;

        let id = trie.insert(sequence.to_str()?.as_bytes());

        peptide_ids.push((sequence.to_string(), id));
    }

    if peptide_ids.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("No peptides!"));
    }

    Ok((
        peptide_ids,
        match (enzyme_patt, n_req_termini) {
            (None, _) | (_, 0) => {
                let iter = match ::preppers::fasta::annotate_fasta(&fasta, trie) {
                    Some(i) => i,
                    None => return Err(pyo3::exceptions::PyValueError::new_err("No peptides!")),
                };

                iter.map(|e| PreppedFastaEntryCopy::from(e)).collect()
            },
            (Some(patt), _) => {
                let regex = match Regex::new(patt.str()?.to_str()?) {
                    Ok(regex) => regex,
                    Err(e) => {
                        return Err(pyo3::exceptions::PyValueError::new_err(format!("Invalid regular expression: {}", e)));
                    }
                };

                let iter = match ::preppers::fasta::annotate_fasta_filtered(&fasta, trie, regex, n_req_termini) {
                    Some(iter) => iter,
                    None => return Err(pyo3::exceptions::PyValueError::new_err("No peptides!")),
                };

                iter.map(|e| PreppedFastaEntryCopy::from(e)).collect()
            },
        }
    ))
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
