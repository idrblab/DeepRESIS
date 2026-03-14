from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ResolvedConfig:
    model_dir: Path
    gene_matrix_dir: Path
    topk: int
    device: str
    seed: int


@dataclass(frozen=True)
class PredictionOutputs:
    output_dir: Path
    predictions_path: Path
    topk_genes_path: Path
    gene_warnings_path: Path
    valid_pair_count: int
    skipped_pair_count: int
