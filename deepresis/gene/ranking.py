from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd


LOGGER = logging.getLogger(__name__)


def rank_genes(
    *,
    predictions: pd.DataFrame,
    gene_matrix_dir: str | Path,
    topk: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_matrix_dir = Path(gene_matrix_dir)
    gene_drug = pd.read_csv(gene_matrix_dir / "drug_gene.csv", index_col=0)
    gene_rna = pd.read_csv(gene_matrix_dir / "ncrna_gene.csv", index_col=0)

    gene_drug.index = gene_drug.index.astype(str)
    gene_drug.columns = gene_drug.columns.astype(str)
    gene_rna.index = gene_rna.index.astype(str)
    gene_rna.columns = gene_rna.columns.astype(str)

    topk_frames: list[pd.DataFrame] = []
    warning_rows: list[dict[str, str]] = []

    unique_pairs = predictions[["ncrna_id", "drug_id"]].drop_duplicates()
    for row in unique_pairs.itertuples(index=False):
        topk_df, warning = rank_genes_for_pair(
            gene_drug=gene_drug,
            gene_rna=gene_rna,
            drug_id=row.drug_id,
            ncrna_id=row.ncrna_id,
            topk=topk,
        )
        if warning is not None:
            LOGGER.warning(warning)
            warning_rows.append(
                {"ncrna_id": row.ncrna_id, "drug_id": row.drug_id, "message": warning}
            )
            continue
        topk_frames.append(topk_df)

    if topk_frames:
        topk_result = pd.concat(topk_frames, axis=0, ignore_index=True)
    else:
        topk_result = pd.DataFrame(
            columns=[
                "ncrna_id",
                "drug_id",
                "rank",
                "gene_id",
                "score",
                "drug_gene_score",
                "ncrna_gene_score",
            ]
        )

    warning_result = pd.DataFrame(
        warning_rows,
        columns=["ncrna_id", "drug_id", "message"],
    )
    return topk_result, warning_result


def rank_genes_for_pair(
    *,
    gene_drug: pd.DataFrame,
    gene_rna: pd.DataFrame,
    drug_id: str,
    ncrna_id: str,
    topk: int,
) -> tuple[pd.DataFrame, str | None]:
    missing: list[str] = []
    if drug_id not in gene_drug.index:
        missing.append(f"drug {drug_id}")
    if ncrna_id not in gene_rna.index:
        missing.append(f"ncRNA {ncrna_id}")
    if missing:
        return pd.DataFrame(), (
            f"Skipped pair ({ncrna_id}, {drug_id}): {' and '.join(missing)} not in network."
        )

    drug_row = gene_drug.loc[drug_id].dropna().astype(float)
    rna_row = gene_rna.loc[ncrna_id].dropna().astype(float)
    common_genes = drug_row.index.intersection(rna_row.index)
    if len(common_genes) == 0:
        return pd.DataFrame(), (
            f"Skipped pair ({ncrna_id}, {drug_id}): "
            "no shared genes between drug-gene and ncrna-gene networks."
        )

    pair_df = pd.DataFrame(
        {
            "gene_id": common_genes,
            "drug_gene_score": drug_row.loc[common_genes].values,
            "ncrna_gene_score": rna_row.loc[common_genes].values,
        }
    )
    pair_df["score"] = (pair_df["drug_gene_score"] + pair_df["ncrna_gene_score"]) / 2.0
    pair_df = pair_df.sort_values("score", ascending=False).head(topk).reset_index(drop=True)
    pair_df["rank"] = pair_df.index + 1
    pair_df["ncrna_id"] = ncrna_id
    pair_df["drug_id"] = drug_id
    return pair_df[
        [
            "ncrna_id",
            "drug_id",
            "rank",
            "gene_id",
            "score",
            "drug_gene_score",
            "ncrna_gene_score",
        ]
    ], None
