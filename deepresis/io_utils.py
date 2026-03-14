from __future__ import annotations

import csv
import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO


LOGGER = logging.getLogger(__name__)
VALID_SEQUENCE_CHARS = set("ACGUTN")
HEADER_TOKENS = {
    "ncrna_id",
    "ncrna",
    "drug_id",
    "pubchem_id",
    "smiles",
    "cid",
}


@dataclass(frozen=True)
class InputData:
    fasta_path: Path
    smiles_path: Path
    pairs_path: Path
    ncrna_records: pd.DataFrame
    smiles_records: pd.DataFrame
    valid_pairs: pd.DataFrame
    skipped_pairs: pd.DataFrame


def prepare_input_data(
    fasta_path: str | Path,
    smiles_path: str | Path,
    pairs_path: str | Path,
) -> InputData:
    fasta = Path(fasta_path).expanduser().resolve()
    smiles = Path(smiles_path).expanduser().resolve()
    pairs = Path(pairs_path).expanduser().resolve()

    ncrna_records = read_fasta_records(fasta)
    smiles_records = read_smiles_records(smiles)
    pair_records = read_pairs_records(pairs)

    ncrna_ids = set(ncrna_records["ncrna_id"])
    drug_ids = set(smiles_records["drug_id"])
    valid_rows: list[dict[str, str]] = []
    skipped_rows: list[dict[str, str]] = []

    for row in pair_records.itertuples(index=False):
        missing: list[str] = []
        if row.ncrna_id not in ncrna_ids:
            missing.append(f"ncRNA {row.ncrna_id!r} not found in FASTA")
        if row.drug_id not in drug_ids:
            missing.append(f"drug {row.drug_id!r} not found in SMILES")
        if missing:
            message = "; ".join(missing)
            LOGGER.warning("Skipping pair (%s, %s): %s", row.ncrna_id, row.drug_id, message)
            skipped_rows.append(
                {
                    "ncrna_id": row.ncrna_id,
                    "drug_id": row.drug_id,
                    "message": message,
                }
            )
            continue
        valid_rows.append({"ncrna_id": row.ncrna_id, "drug_id": row.drug_id})

    valid_pairs = pd.DataFrame(valid_rows, columns=["ncrna_id", "drug_id"])
    skipped_pairs = pd.DataFrame(skipped_rows, columns=["ncrna_id", "drug_id", "message"])
    if valid_pairs.empty:
        raise ValueError("No valid ncRNA-drug pairs remain after input validation.")

    return InputData(
        fasta_path=fasta,
        smiles_path=smiles,
        pairs_path=pairs,
        ncrna_records=ncrna_records,
        smiles_records=smiles_records,
        valid_pairs=valid_pairs,
        skipped_pairs=skipped_pairs,
    )


def read_fasta_records(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for record in SeqIO.parse(str(path), "fasta"):
        record_id = str(record.id).strip()
        sequence = str(record.seq).strip().upper()
        if not record_id:
            raise ValueError(f"Found FASTA entry with empty ID in {path}")
        if record_id in seen:
            raise ValueError(f"Duplicate ncRNA ID {record_id!r} in FASTA file {path}")
        if not sequence:
            raise ValueError(f"FASTA sequence for {record_id!r} is empty in {path}")
        illegal = sorted(set(sequence) - VALID_SEQUENCE_CHARS)
        if illegal:
            raise ValueError(
                f"FASTA sequence for {record_id!r} contains invalid characters {illegal} in {path}"
            )
        rows.append({"ncrna_id": record_id, "sequence": sequence})
        seen.add(record_id)

    if not rows:
        raise ValueError(f"No FASTA records found in {path}")
    return pd.DataFrame(rows)


def read_smiles_records(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"SMILES file not found: {path}")

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if _is_header(parts):
            continue
        if len(parts) < 2:
            raise ValueError(f"Invalid SMILES format at {path}:{line_no}: expected 2 columns")
        drug_id, smiles = parts[0], parts[1]
        if drug_id in seen:
            raise ValueError(f"Duplicate drug ID {drug_id!r} in SMILES file {path}")
        rows.append({"drug_id": drug_id, "smiles": smiles})
        seen.add(drug_id)

    if not rows:
        raise ValueError(f"No SMILES records found in {path}")
    return pd.DataFrame(rows)


def read_pairs_records(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Pairs file not found: {path}")

    rows: list[dict[str, str]] = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if _is_header(parts):
            continue
        if len(parts) < 2:
            raise ValueError(f"Invalid pairs format at {path}:{line_no}: expected 2 columns")
        rows.append({"ncrna_id": parts[0], "drug_id": parts[1]})

    if not rows:
        raise ValueError(f"No pairs found in {path}")
    return pd.DataFrame(rows)


def write_normalized_fasta(records: pd.DataFrame, destination: str | Path) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for row in records.itertuples(index=False):
            handle.write(f">{row.ncrna_id}\n{row.sequence}\n")
    return destination


def write_normalized_smiles(records: pd.DataFrame, destination: str | Path) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for row in records.itertuples(index=False):
            handle.write(f"{row.drug_id}\t{row.smiles}\n")
    return destination


def create_work_files(
    ncrna_records: pd.DataFrame,
    smiles_records: pd.DataFrame,
) -> tuple[tempfile.TemporaryDirectory[str], Path, Path]:
    temp_dir = tempfile.TemporaryDirectory(prefix="deepresis_")
    tmp_path = Path(temp_dir.name)
    fasta_path = write_normalized_fasta(ncrna_records, tmp_path / "input.fasta")
    smiles_path = write_normalized_smiles(smiles_records, tmp_path / "input_smiles.tsv")
    return temp_dir, fasta_path, smiles_path


def write_tsv(dataframe: pd.DataFrame, destination: str | Path) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(destination, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    return destination


def _is_header(parts: Iterable[str]) -> bool:
    tokens = [part.strip().lower() for part in parts]
    if len(tokens) < 2:
        return False
    return tokens[0] in HEADER_TOKENS and tokens[1] in HEADER_TOKENS
