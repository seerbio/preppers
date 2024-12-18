import preppers

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


def test_preppers_basic_bytes(peptides, fasta_bytes):
    preppers.annotate_fasta_bytes(
        peptides,
        fasta_bytes,
    )