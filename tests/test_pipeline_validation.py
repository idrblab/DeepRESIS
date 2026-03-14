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
    (root / "drug_gene.csv").write_text(",g1,g2\nd1,0.8,0.1\nd2,0.6,0.4\n", encoding="utf-8")
    (root / "ncrna_gene.csv").write_text(",g1,g2\nr1,0.7,0.3\nr2,0.1,0.9\n", encoding="utf-8")
    return root


def test_pipeline_outputs_requested_pairs_only(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fasta"
    smiles = tmp_path / "input.txt"
    pairs = tmp_path / "pairs.txt"
    fasta.write_text(">r1\nACTG\n>r2\nACTG\n", encoding="utf-8")
    smiles.write_text("d1 CCO\nd2 CCC\n", encoding="utf-8")
    pairs.write_text("r1 d1\nr1 d1\nr2 d2\n", encoding="utf-8")

    monkeypatch.setattr(
        "deepresis.pipeline._generate_ncrna_features",
        lambda path: pd.DataFrame([[1.0] * 638, [2.0] * 638], index=["r1", "r2"]),
    )
    monkeypatch.setattr(
        "deepresis.pipeline._generate_drug_descriptors",
        lambda path: pd.DataFrame([[1.0] * 1345, [2.0] * 1345], index=["d1", "d2"]),
    )
    monkeypatch.setattr(
        "deepresis.pipeline._generate_drug_fingerprints",
        lambda path: pd.DataFrame([[1] * 16204, [0] * 16204], index=["d1", "d2"]),
    )
    monkeypatch.setattr(
        "deepresis.pipeline._predict_resistance",
        lambda **kwargs: pd.DataFrame(
            [
                {"ncrna_id": "r1", "drug_id": "d1", "score": 0.9, "label": "resistant", "label_id": 1},
                {"ncrna_id": "r1", "drug_id": "d1", "score": 0.9, "label": "resistant", "label_id": 1},
                {"ncrna_id": "r2", "drug_id": "d2", "score": 0.1, "label": "non-resistant", "label_id": 0},
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
    predictions = pd.read_csv(outputs.predictions_path, sep="\t")
    assert predictions.shape[0] == 3
    assert predictions[["ncrna_id", "drug_id"]].values.tolist() == [["r1", "d1"], ["r1", "d1"], ["r2", "d2"]]
    topk = pd.read_csv(outputs.topk_genes_path, sep="\t")
    assert topk["rank"].max() == 1
