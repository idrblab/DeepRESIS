from __future__ import annotations

import pytest

from deepresis.io_utils import read_fasta_records, read_pairs_records, read_smiles_records


def test_read_fasta_records(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">a\nACTG\n>b\nAUGC\n", encoding="utf-8")
    df = read_fasta_records(fasta)
    assert df["ncrna_id"].tolist() == ["a", "b"]
    assert df["sequence"].tolist() == ["ACTG", "AUGC"]


def test_read_fasta_records_duplicate_id(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">a\nACTG\n>a\nACTG\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate ncRNA ID"):
        read_fasta_records(fasta)


def test_read_fasta_records_empty_sequence(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">a\n\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No FASTA records found|empty"):
        read_fasta_records(fasta)


def test_read_smiles_records_space_and_tab(tmp_path):
    smiles = tmp_path / "smiles.txt"
    smiles.write_text("drug1 CCO\ndrug2\tCCN(CC)CC\n", encoding="utf-8")
    df = read_smiles_records(smiles)
    assert df["drug_id"].tolist() == ["drug1", "drug2"]


def test_read_smiles_records_invalid_format(tmp_path):
    smiles = tmp_path / "smiles.txt"
    smiles.write_text("drug1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Invalid SMILES format"):
        read_smiles_records(smiles)


def test_read_pairs_records_skip_header(tmp_path):
    pairs = tmp_path / "pairs.txt"
    pairs.write_text("ncrna_id\tdrug_id\na\td1\n", encoding="utf-8")
    df = read_pairs_records(pairs)
    assert df.to_dict("records") == [{"ncrna_id": "a", "drug_id": "d1"}]
