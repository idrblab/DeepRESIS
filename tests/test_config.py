from __future__ import annotations

from pathlib import Path

import pytest

from deepresis.config import ConfigError, load_config


@pytest.fixture(autouse=True)
def _isolate_config_resolution(monkeypatch, tmp_path):
    for key in (
        "DEEPRESIS_CONFIG",
        "DEEPRESIS_MODEL_DIR",
        "DEEPRESIS_GENE_MATRIX_DIR",
    ):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.chdir(tmp_path)


def _make_model_dir(root: Path) -> Path:
    for fold in range(5):
        fold_dir = root / f"fold{fold}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        (fold_dir / "ne_student.ckpt").write_text("stub", encoding="utf-8")
    return root


def _make_gene_dir(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "drug_gene.csv").write_text(",g1\n1,0.1\n", encoding="utf-8")
    (root / "ncrna_gene.csv").write_text(",g1\na,0.2\n", encoding="utf-8")
    return root


def test_load_config_from_file(tmp_path):
    model_dir = _make_model_dir(tmp_path / "models")
    gene_dir = _make_gene_dir(tmp_path / "genes")
    config = tmp_path / "deepresis.toml"
    config.write_text(
        f"[paths]\nmodel_dir = \"{model_dir}\"\ngene_matrix_dir = \"{gene_dir}\"\n[run]\nseed = 7\ntopk = 10\ndevice = \"cpu\"\n",
        encoding="utf-8",
    )
    resolved = load_config(config_path=config)
    assert resolved.model_dir == model_dir.resolve()
    assert resolved.gene_matrix_dir == gene_dir.resolve()
    assert resolved.seed == 7
    assert resolved.topk == 10
    assert resolved.device == "cpu"


def test_cli_overrides_config(tmp_path):
    model_dir = _make_model_dir(tmp_path / "models")
    gene_dir = _make_gene_dir(tmp_path / "genes")
    resolved = load_config(model_dir=model_dir, gene_matrix_dir=gene_dir, topk=9, device="cpu", seed=5)
    assert resolved.topk == 9
    assert resolved.device == "cpu"
    assert resolved.seed == 5


def test_missing_model_dir_raises(tmp_path):
    gene_dir = _make_gene_dir(tmp_path / "genes")
    with pytest.raises(ConfigError, match="Model directory is not configured"):
        load_config(gene_matrix_dir=gene_dir)


def test_missing_gene_matrix_dir_raises(tmp_path):
    model_dir = _make_model_dir(tmp_path / "models")
    with pytest.raises(ConfigError, match="Gene matrix directory is not configured"):
        load_config(model_dir=model_dir)


def test_missing_checkpoint_raises(tmp_path):
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    gene_dir = _make_gene_dir(tmp_path / "genes")
    with pytest.raises(ConfigError, match="Missing model checkpoint files"):
        load_config(model_dir=model_dir, gene_matrix_dir=gene_dir)


def test_missing_gene_matrix_file_raises(tmp_path):
    model_dir = _make_model_dir(tmp_path / "models")
    gene_dir = tmp_path / "genes"
    gene_dir.mkdir()
    with pytest.raises(ConfigError, match="Missing gene matrix files"):
        load_config(model_dir=model_dir, gene_matrix_dir=gene_dir)
