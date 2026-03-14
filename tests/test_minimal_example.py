from __future__ import annotations

import pandas as pd

from deepresis.pipeline import run_prediction


def _make_model_dir(root):
    for fold in range(5):
        fold_dir = root / f"fold{fold}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        (fold_dir / "ne_student.ckpt").write_text("stub", encoding="utf-8")
    return root


def _make_gene_dir(root):
    root.mkdir(parents=True, exist_ok=True)
    (root / "drug_gene.csv").write_text(",g1\ndrug1,0.8\n", encoding="utf-8")
    (root / "ncrna_gene.csv").write_text(",g1\nrna1,0.6\n", encoding="utf-8")
    return root


def test_minimal_example_runs(tmp_path, monkeypatch):
    fasta = tmp_path / "test.fasta"
    smiles = tmp_path / "smiles.txt"
    pairs = tmp_path / "pairs.txt"
    fasta.write_text(">rna1\nACTG\n", encoding="utf-8")
    smiles.write_text("drug1 CCO\n", encoding="utf-8")
    pairs.write_text("rna1 drug1\n", encoding="utf-8")

    monkeypatch.setattr(
        "deepresis.pipeline._generate_ncrna_features",
        lambda path: pd.DataFrame([[1.0] * 638], index=["rna1"]),
    )
    monkeypatch.setattr(
        "deepresis.pipeline._generate_drug_descriptors",
        lambda path: pd.DataFrame([[1.0] * 1345], index=["drug1"]),
    )
    monkeypatch.setattr(
        "deepresis.pipeline._generate_drug_fingerprints",
        lambda path: pd.DataFrame([[1] * 16204], index=["drug1"]),
    )
    monkeypatch.setattr(
        "deepresis.pipeline._predict_resistance",
        lambda **kwargs: pd.DataFrame(
            [
                {
                    "ncrna_id": "rna1",
                    "drug_id": "drug1",
                    "score": 0.75,
                    "label": "resistant",
                    "label_id": 1,
                }
            ]
        ),
    )

    outputs = run_prediction(
        fasta_path=fasta,
        smiles_path=smiles,
        pairs_path=pairs,
        output_dir=tmp_path / "outputs",
        model_dir=_make_model_dir(tmp_path / "models"),
        gene_matrix_dir=_make_gene_dir(tmp_path / "genes"),
        device="cpu",
        topk=1,
    )
    assert outputs.predictions_path.exists()
    assert outputs.topk_genes_path.exists()
    assert outputs.gene_warnings_path.exists()
    predictions = pd.read_csv(outputs.predictions_path, sep="\t")
    assert predictions.to_dict("records") == [
        {
            "ncrna_id": "rna1",
            "drug_id": "drug1",
            "score": 0.75,
            "label": "resistant",
            "label_id": 1,
        }
    ]
