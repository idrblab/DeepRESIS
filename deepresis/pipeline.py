from __future__ import annotations

import logging
from pathlib import Path

from .config import load_config
from .io_utils import create_work_files, prepare_input_data, write_tsv
from .logging_utils import configure_logging
from .types import PredictionOutputs


LOGGER = logging.getLogger(__name__)


def run_prediction(
    fasta_path: str | Path,
    smiles_path: str | Path,
    pairs_path: str | Path,
    output_dir: str | Path,
    topk: int | None = None,
    device: str | None = None,
    seed: int | None = None,
    model_dir: str | Path | None = None,
    gene_matrix_dir: str | Path | None = None,
    config_path: str | Path | None = None,
    verbose: bool = False,
) -> PredictionOutputs:
    configure_logging(verbose)
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    config = load_config(
        config_path=config_path,
        model_dir=model_dir,
        gene_matrix_dir=gene_matrix_dir,
        topk=topk,
        device=device,
        seed=seed,
    )
    input_data = prepare_input_data(fasta_path, smiles_path, pairs_path)

    temp_dir, normalized_fasta, normalized_smiles = create_work_files(
        input_data.ncrna_records,
        input_data.smiles_records,
    )
    try:
        ncrna_feature = _generate_ncrna_features(normalized_fasta)
        drug_descriptors = _generate_drug_descriptors(normalized_smiles)
        drug_fingerprint = _generate_drug_fingerprints(normalized_smiles)

        predictions = _predict_resistance(
            ncrna_feature=ncrna_feature,
            drug_descriptors=drug_descriptors,
            drug_fingerprint=drug_fingerprint,
            pairs=input_data.valid_pairs,
            model_dir=config.model_dir,
            device=config.device,
            seed=config.seed,
        )
        topk_result, gene_warnings = _rank_genes(
            predictions=predictions,
            gene_matrix_dir=config.gene_matrix_dir,
            topk=config.topk,
        )
    finally:
        temp_dir.cleanup()

    predictions_path = write_tsv(predictions, output_dir / "resistance_predictions.tsv")
    topk_genes_path = write_tsv(topk_result, output_dir / "topk_genes.tsv")
    gene_warnings_path = write_tsv(gene_warnings, output_dir / "gene_ranking_warnings.tsv")

    LOGGER.info("Prediction results written to %s", predictions_path)
    LOGGER.info("Top-k gene results written to %s", topk_genes_path)
    LOGGER.info("Gene warnings written to %s", gene_warnings_path)

    return PredictionOutputs(
        output_dir=output_dir,
        predictions_path=predictions_path,
        topk_genes_path=topk_genes_path,
        gene_warnings_path=gene_warnings_path,
        valid_pair_count=len(input_data.valid_pairs),
        skipped_pair_count=len(input_data.skipped_pairs),
    )


def _generate_ncrna_features(fasta_path: Path):
    from deepresis.legacy.ncresis_backend.ncrna_feature import get_ncrna_feature

    return get_ncrna_feature(str(fasta_path))


def _generate_drug_descriptors(smiles_path: Path):
    from deepresis.legacy.ncresis_backend.drug_descriptor import get_descriptor

    return get_descriptor(str(smiles_path))


def _generate_drug_fingerprints(smiles_path: Path):
    from deepresis.legacy.ncresis_backend.drug_fingerprint import get_fingerprint

    return get_fingerprint(str(smiles_path))


def _predict_resistance(**kwargs):
    from .predictors.resistance import predict_resistance

    return predict_resistance(**kwargs)


def _rank_genes(**kwargs):
    from .gene.ranking import rank_genes

    return rank_genes(**kwargs)
