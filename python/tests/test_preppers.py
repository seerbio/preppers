import preppers

def test_preppers_basic_bytes(peptides, fasta_bytes):
    result = preppers.annotate_fasta_bytes(
        peptides,
        fasta_bytes,
    )
    assert_annotation_result(result)

def test_annotate_fasta(peptides, fasta_file):
    result = preppers.annotate_fasta(
        peptides,
        fasta_file,
    )
    assert_annotation_result(result)

def assert_annotation_result(result):
    assert isinstance(result, tuple), "Result should be a tuple"
    assert len(result) == 2, "Result tuple should have two elements"

    peptides, proteins = result

    # Check peptides
    assert isinstance(peptides, list), "Peptides should be a list"
    for peptide in peptides:
        assert isinstance(peptide, tuple), "Each peptide should be a tuple"
        assert len(peptide) == 2, "Each peptide tuple should have two elements"
        assert isinstance(peptide[0], str), "First element of peptide tuple should be a string"
        assert isinstance(peptide[1], int), "Second element of peptide tuple should be an integer"

    # Check proteins
    assert isinstance(proteins, list), "Proteins should be a list"
    for protein in proteins:
        assert hasattr(protein, 'header'), "Protein object should have a header attribute"
        assert hasattr(protein, 'sequence'), "Protein object should have a sequence attribute"
        assert hasattr(protein, 'peptides'), "Protein object should have a peptides attribute"
        assert isinstance(protein.header, str), "Protein header should be a string"
        assert isinstance(protein.sequence, str), "Protein sequence should be a string"
        assert isinstance(protein.peptides, list), "Protein peptides should be a list"
        for peptide_id in protein.peptides:
            assert isinstance(peptide_id, int), "Each peptide ID should be an integer"

import pytest

@pytest.fixture()
def peptides():
    return [
        "APEPTIDEK",
        "ANOTHER",
    ]

import pytest

@pytest.fixture()
def fasta_bytes() -> bytes:
    fasta_content = """>sequence1
APEPTIDEKANOTHER
"""
    return fasta_content.encode('utf-8')

@pytest.fixture()
def fasta_file(tmp_path, fasta_bytes):
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "wb") as f:
        f.write(fasta_bytes)
    return fasta_path
