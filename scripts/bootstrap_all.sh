#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${1:-deepresis}"
ASSET_ROOT="${2:-${ROOT_DIR}/artifacts}"
CONFIG_PATH="${3:-${ROOT_DIR}/deepresis.toml}"

echo "[bootstrap-all] Step 1/2: environment"
bash "${ROOT_DIR}/scripts/bootstrap_env.sh" "${ENV_NAME}"

echo "[bootstrap-all] Step 2/2: assets"
bash "${ROOT_DIR}/scripts/bootstrap_assets.sh" "${ASSET_ROOT}" "${CONFIG_PATH}"

cat <<EOF

All bootstrap steps finished.

Next:
  conda activate ${ENV_NAME}
  cd ${ROOT_DIR}
  deepresis --help
  deepresis predict --config ${CONFIG_PATH} --fasta <fasta> --smiles <smiles> --pairs <pairs> --output-dir <output_dir>
EOF
