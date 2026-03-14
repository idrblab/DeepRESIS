from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from deepresis.models.transformer import Classifier


LOGGER = logging.getLogger(__name__)
FOLDS = 5


def predict_resistance(
    *,
    ncrna_feature: pd.DataFrame,
    drug_descriptors: pd.DataFrame,
    drug_fingerprint: pd.DataFrame,
    pairs: pd.DataFrame,
    model_dir: str | Path,
    device: str,
    seed: int,
) -> pd.DataFrame:
    torch.manual_seed(seed)
    np.random.seed(seed)

    resolved_device = _resolve_device(device)
    model = build_model(resolved_device)
    states = load_checkpoints(model_dir, resolved_device)

    ncrna_feature = _sanitize_features(ncrna_feature)
    drug_descriptors = _sanitize_features(drug_descriptors)
    drug_fingerprint = _sanitize_features(drug_fingerprint)

    result = pairs.copy().reset_index(drop=True)
    result["score"] = np.nan
    result["label"] = ""
    result["label_id"] = -1
    for ncrna_id, group in pairs.groupby("ncrna_id", sort=False):
        unique_drug_ids = list(dict.fromkeys(group["drug_id"].tolist()))
        LOGGER.info("Predicting %d pairs for ncRNA %s", len(unique_drug_ids), ncrna_id)

        ncrna_vec = ncrna_feature.loc[ncrna_id].values.astype(np.float32)
        ncrna_batch = np.repeat(ncrna_vec[None, None, :], len(unique_drug_ids), axis=0)

        desc_values = torch.tensor(
            drug_descriptors.loc[unique_drug_ids].values[:, None, :],
            dtype=torch.float32,
            device=resolved_device,
        )
        fp_values = torch.tensor(
            drug_fingerprint.loc[unique_drug_ids].values[:, None, :],
            dtype=torch.float32,
            device=resolved_device,
        )
        ncrna_tensor = torch.tensor(ncrna_batch, dtype=torch.float32, device=resolved_device)

        fold_preds: list[np.ndarray] = []
        with torch.no_grad():
            for state in states:
                model.load_state_dict(state)
                pred = model(ncrna_tensor, desc_values, fp_values)
                fold_preds.append(pred.squeeze(1).detach().cpu().numpy())

        mean_pred = np.mean(np.stack(fold_preds, axis=0), axis=0)
        score_map = {drug_id: float(score) for drug_id, score in zip(unique_drug_ids, mean_pred)}

        for index in group.index:
            score = score_map[pairs.loc[index, "drug_id"]]
            label_id = int(score >= 0.5)
            result.loc[index, "score"] = score
            result.loc[index, "label"] = "resistant" if label_id else "non-resistant"
            result.loc[index, "label_id"] = label_id

    return result[["ncrna_id", "drug_id", "score", "label", "label_id"]]


def build_model(device: torch.device) -> Classifier:
    model = Classifier(
        dim_1=638,
        dim_2=1345,
        dim_3=16204,
        d_model=256,
        dropout=0.01,
        nhead=8,
        num_encoder_layers=12,
        num_decoder_layers=12,
    ).to(device)
    model.eval()
    return model


def load_checkpoints(model_dir: str | Path, device: torch.device) -> list[dict[str, torch.Tensor]]:
    model_dir = Path(model_dir)
    states = []
    for fold in range(FOLDS):
        ckpt = model_dir / f"fold{fold}" / "ne_student.ckpt"
        LOGGER.debug("Loading checkpoint %s", ckpt)
        states.append(torch.load(ckpt, map_location=device))
    return states


def _resolve_device(device: str) -> torch.device:
    if device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available.")
    return torch.device(device)


def _sanitize_features(dataframe: pd.DataFrame) -> pd.DataFrame:
    result = dataframe.copy()
    result.index = result.index.astype(str)
    return result.replace([np.inf, -np.inf], np.nan).fillna(0)
