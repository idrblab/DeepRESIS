#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_ROOT="${1:-${ROOT_DIR}/artifacts}"
CONFIG_PATH="${2:-${ROOT_DIR}/deepresis.toml}"

MODELS_DIR="${TARGET_ROOT}/models"
GENE_MATRIX_DIR="${TARGET_ROOT}/gene_matrix"

DEFAULT_MODELS_URL="https://huggingface.co/swallow-design/DeepRESIS/tree/main/model_pharameter"
DEFAULT_GENE_MATRIX_URL="https://huggingface.co/swallow-design/DeepRESIS/tree/main/gene_matrix"

MODELS_URL="${DEEPRESIS_MODELS_URL:-${DEFAULT_MODELS_URL}}"
GENE_MATRIX_URL="${DEEPRESIS_GENE_MATRIX_URL:-${DEFAULT_GENE_MATRIX_URL}}"

if command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="python"
else
  echo "python3 or python is required to download Hugging Face assets" >&2
  exit 1
fi

validate_models() {
  local model_root="$1"
  local missing=0

  for fold in 0 1 2 3 4; do
    if [[ ! -f "${model_root}/fold${fold}/ne_student.ckpt" ]]; then
      echo "Missing model checkpoint: ${model_root}/fold${fold}/ne_student.ckpt" >&2
      missing=1
    fi
  done

  if [[ "${missing}" -ne 0 ]]; then
    exit 1
  fi
}

validate_gene_matrix() {
  local gene_root="$1"

  if [[ ! -f "${gene_root}/drug_gene.csv" ]]; then
    echo "Missing gene matrix file: ${gene_root}/drug_gene.csv" >&2
    exit 1
  fi

  if [[ ! -f "${gene_root}/ncrna_gene.csv" ]]; then
    echo "Missing gene matrix file: ${gene_root}/ncrna_gene.csv" >&2
    exit 1
  fi
}

write_config() {
  local models_dir="$1"
  local gene_dir="$2"
  local config_path="$3"

  mkdir -p "$(dirname "${config_path}")"

  cat > "${config_path}" <<EOF
[paths]
model_dir = "${models_dir}"
gene_matrix_dir = "${gene_dir}"

[run]
device = "auto"
seed = 123
topk = 300
EOF
}

mkdir -p "${MODELS_DIR}" "${GENE_MATRIX_DIR}"

echo "[1/4] Preparing asset directories under ${TARGET_ROOT}"
echo "[2/4] Downloading required model checkpoints and gene files from Hugging Face"

"${PYTHON_BIN}" - "${MODELS_URL}" "${GENE_MATRIX_URL}" "${MODELS_DIR}" "${GENE_MATRIX_DIR}" <<'PY'
from __future__ import annotations

import shutil
import sys
import urllib.request
from pathlib import Path


MODEL_FILES = [
    "fold0/ne_student.ckpt",
    "fold1/ne_student.ckpt",
    "fold2/ne_student.ckpt",
    "fold3/ne_student.ckpt",
    "fold4/ne_student.ckpt",
]

GENE_FILES = [
    "drug_gene.csv",
    "ncrna_gene.csv",
]


def normalize_base(url: str) -> str:
    normalized = url.rstrip("/")
    normalized = normalized.replace("/tree/main/", "/resolve/main/")
    normalized = normalized.replace("/blob/main/", "/resolve/main/")
    return normalized


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        print(f"Skip existing file: {destination}")
        return

    tmp_destination = destination.with_suffix(destination.suffix + ".part")
    print(f"Downloading {url}")
    try:
        with urllib.request.urlopen(url) as response, tmp_destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    except Exception as exc:
        if tmp_destination.exists():
            tmp_destination.unlink()
        raise SystemExit(f"Failed to download {url} -> {destination}: {exc}") from exc

    tmp_destination.replace(destination)


def main() -> int:
    model_base = normalize_base(sys.argv[1])
    gene_base = normalize_base(sys.argv[2])
    model_dir = Path(sys.argv[3])
    gene_dir = Path(sys.argv[4])

    for relative_path in MODEL_FILES:
        download(f"{model_base}/{relative_path}", model_dir / relative_path)

    for relative_path in GENE_FILES:
        download(f"{gene_base}/{relative_path}", gene_dir / relative_path)

    return 0


raise SystemExit(main())
PY

echo "[3/4] Validating asset layout"
validate_models "${MODELS_DIR}"
validate_gene_matrix "${GENE_MATRIX_DIR}"

echo "[4/4] Writing config file to ${CONFIG_PATH}"
write_config "${MODELS_DIR}" "${GENE_MATRIX_DIR}" "${CONFIG_PATH}"

cat <<EOF

Asset bootstrap finished.

Config file:
  ${CONFIG_PATH}

Resolved asset paths:
  models      ${MODELS_DIR}
  gene_matrix ${GENE_MATRIX_DIR}

Example:
  deepresis predict --config ${CONFIG_PATH} --fasta <fasta> --smiles <smiles> --pairs <pairs> --output-dir <output_dir>
EOF
