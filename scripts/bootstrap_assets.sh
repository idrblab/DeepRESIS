#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_ROOT="${1:-${ROOT_DIR}/artifacts}"
CONFIG_PATH="${2:-${ROOT_DIR}/deepresis.toml}"

MODELS_DIR="${TARGET_ROOT}/models"
GENE_MATRIX_DIR="${TARGET_ROOT}/gene_matrix"
TMP_DIR="${TARGET_ROOT}/.tmp"

MODELS_URL="${DEEPRESIS_MODELS_URL:-}"
GENE_MATRIX_URL="${DEEPRESIS_GENE_MATRIX_URL:-}"
LOCAL_MODEL_DIR="${DEEPRESIS_LOCAL_MODEL_DIR:-}"
LOCAL_GENE_MATRIX_DIR="${DEEPRESIS_LOCAL_GENE_MATRIX_DIR:-}"

download_and_extract() {
  local url="$1"
  local dest="$2"

  python3 - "$url" "$dest" <<'PY'
from __future__ import annotations

import shutil
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path


def main() -> int:
    url = sys.argv[1]
    dest = Path(sys.argv[2])
    dest.mkdir(parents=True, exist_ok=True)

    archive_name = url.rsplit("/", 1)[-1] or "downloaded_asset"
    archive_path = dest / archive_name

    with urllib.request.urlopen(url) as response, archive_path.open("wb") as handle:
        shutil.copyfileobj(response, handle)

    suffixes = archive_path.suffixes
    name = archive_path.name.lower()

    if name.endswith(".tar.gz") or name.endswith(".tgz"):
        with tarfile.open(archive_path, "r:gz") as tar:
            tar.extractall(dest)
    elif name.endswith(".tar"):
        with tarfile.open(archive_path, "r:") as tar:
            tar.extractall(dest)
    elif name.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(dest)
    else:
        raise SystemExit(
            "Unsupported asset archive format. Use .tar.gz, .tgz, .tar, or .zip."
        )

    archive_path.unlink()
    return 0


raise SystemExit(main())
PY
}

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

rm -rf "${TMP_DIR}"
mkdir -p "${MODELS_DIR}" "${GENE_MATRIX_DIR}" "${TMP_DIR}"

echo "[1/4] Preparing asset directories under ${TARGET_ROOT}"

if [[ -n "${LOCAL_MODEL_DIR}" && -n "${LOCAL_GENE_MATRIX_DIR}" ]]; then
  echo "[2/4] Copying local model and gene assets"
  rm -rf "${MODELS_DIR}" "${GENE_MATRIX_DIR}"
  mkdir -p "${MODELS_DIR}" "${GENE_MATRIX_DIR}"
  cp -a "${LOCAL_MODEL_DIR}/." "${MODELS_DIR}/"
  cp -a "${LOCAL_GENE_MATRIX_DIR}/." "${GENE_MATRIX_DIR}/"
elif [[ -n "${MODELS_URL}" && -n "${GENE_MATRIX_URL}" ]]; then
  echo "[2/4] Downloading model and gene assets"
  rm -rf "${MODELS_DIR}" "${GENE_MATRIX_DIR}"
  mkdir -p "${MODELS_DIR}" "${GENE_MATRIX_DIR}"
  download_and_extract "${MODELS_URL}" "${MODELS_DIR}"
  download_and_extract "${GENE_MATRIX_URL}" "${GENE_MATRIX_DIR}"
else
  echo "Asset source is not configured." >&2
  echo "Use either:" >&2
  echo "  DEEPRESIS_LOCAL_MODEL_DIR + DEEPRESIS_LOCAL_GENE_MATRIX_DIR" >&2
  echo "or:" >&2
  echo "  DEEPRESIS_MODELS_URL + DEEPRESIS_GENE_MATRIX_URL" >&2
  exit 1
fi

echo "[3/4] Validating asset layout"
validate_models "${MODELS_DIR}"
validate_gene_matrix "${GENE_MATRIX_DIR}"

echo "[4/4] Writing config file to ${CONFIG_PATH}"
write_config "${MODELS_DIR}" "${GENE_MATRIX_DIR}" "${CONFIG_PATH}"

rm -rf "${TMP_DIR}"

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
