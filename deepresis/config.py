from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from .types import ResolvedConfig

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib


DEFAULT_TOPK = 300
DEFAULT_DEVICE = "auto"
DEFAULT_SEED = 123
SUPPORTED_DEVICES = {"auto", "cpu", "cuda"}


class ConfigError(ValueError):
    """Raised when configuration is invalid."""


def load_config(
    *,
    config_path: str | Path | None = None,
    model_dir: str | Path | None = None,
    gene_matrix_dir: str | Path | None = None,
    topk: int | None = None,
    device: str | None = None,
    seed: int | None = None,
) -> ResolvedConfig:
    config_data = _load_config_file(config_path)

    resolved_model_dir = _resolve_path_value(
        cli_value=model_dir,
        env_var="DEEPRESIS_MODEL_DIR",
        config_value=_nested_get(config_data, "paths", "model_dir"),
    )
    resolved_gene_matrix_dir = _resolve_path_value(
        cli_value=gene_matrix_dir,
        env_var="DEEPRESIS_GENE_MATRIX_DIR",
        config_value=_nested_get(config_data, "paths", "gene_matrix_dir"),
    )
    resolved_topk = _resolve_scalar_value(
        cli_value=topk,
        config_value=_nested_get(config_data, "run", "topk"),
        default=DEFAULT_TOPK,
    )
    resolved_device = _resolve_scalar_value(
        cli_value=device,
        config_value=_nested_get(config_data, "run", "device"),
        default=DEFAULT_DEVICE,
    )
    resolved_seed = _resolve_scalar_value(
        cli_value=seed,
        config_value=_nested_get(config_data, "run", "seed"),
        default=DEFAULT_SEED,
    )

    if resolved_model_dir is None:
        raise ConfigError(
            "Model directory is not configured. Set --model-dir, "
            "DEEPRESIS_MODEL_DIR, or paths.model_dir in deepresis.toml."
        )
    if resolved_gene_matrix_dir is None:
        raise ConfigError(
            "Gene matrix directory is not configured. Set --gene-matrix-dir, "
            "DEEPRESIS_GENE_MATRIX_DIR, or paths.gene_matrix_dir in deepresis.toml."
        )
    if int(resolved_topk) <= 0:
        raise ConfigError(f"topk must be > 0, got {resolved_topk!r}.")
    if str(resolved_device) not in SUPPORTED_DEVICES:
        raise ConfigError(
            f"Unsupported device {resolved_device!r}. Choose from: "
            f"{', '.join(sorted(SUPPORTED_DEVICES))}."
        )

    model_dir_path = Path(resolved_model_dir).expanduser().resolve()
    gene_matrix_dir_path = Path(resolved_gene_matrix_dir).expanduser().resolve()

    _validate_model_dir(model_dir_path)
    _validate_gene_matrix_dir(gene_matrix_dir_path)

    return ResolvedConfig(
        model_dir=model_dir_path,
        gene_matrix_dir=gene_matrix_dir_path,
        topk=int(resolved_topk),
        device=str(resolved_device),
        seed=int(resolved_seed),
    )


def _load_config_file(config_path: str | Path | None) -> dict[str, Any]:
    candidate_paths: list[Path] = []
    env_config = os.getenv("DEEPRESIS_CONFIG")
    if config_path is not None:
        candidate_paths.append(Path(config_path).expanduser())
    elif env_config:
        candidate_paths.append(Path(env_config).expanduser())
    else:
        candidate_paths.extend(
            [
                Path.cwd() / "deepresis.toml",
                Path.home() / ".config" / "deepresis" / "config.toml",
            ]
        )

    for candidate in candidate_paths:
        if candidate.exists():
            with candidate.open("rb") as handle:
                data = tomllib.load(handle)
            if not isinstance(data, dict):
                raise ConfigError(f"Configuration file {candidate} must define a TOML table.")
            return data
    return {}


def _resolve_path_value(
    *,
    cli_value: str | Path | None,
    env_var: str,
    config_value: Any,
) -> str | Path | None:
    if cli_value is not None:
        return cli_value
    env_value = os.getenv(env_var)
    if env_value:
        return env_value
    return config_value


def _resolve_scalar_value(*, cli_value: Any, config_value: Any, default: Any) -> Any:
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return default


def _nested_get(data: dict[str, Any], *keys: str) -> Any:
    current: Any = data
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def _validate_model_dir(model_dir: Path) -> None:
    if not model_dir.is_dir():
        raise ConfigError(f"Model directory does not exist: {model_dir}")
    missing = [
        str(model_dir / f"fold{fold}" / "ne_student.ckpt")
        for fold in range(5)
        if not (model_dir / f"fold{fold}" / "ne_student.ckpt").is_file()
    ]
    if missing:
        raise ConfigError(
            "Missing model checkpoint files:\n" + "\n".join(missing)
        )


def _validate_gene_matrix_dir(gene_matrix_dir: Path) -> None:
    if not gene_matrix_dir.is_dir():
        raise ConfigError(f"Gene matrix directory does not exist: {gene_matrix_dir}")
    missing = [
        str(path)
        for path in (
            gene_matrix_dir / "drug_gene.csv",
            gene_matrix_dir / "ncrna_gene.csv",
        )
        if not path.is_file()
    ]
    if missing:
        raise ConfigError("Missing gene matrix files:\n" + "\n".join(missing))
