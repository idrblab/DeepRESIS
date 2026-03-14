#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${1:-deepresis}"
ENV_FILE="${ROOT_DIR}/environment.yml"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found on PATH" >&2
  exit 1
fi

echo "[1/5] Creating or updating conda environment: ${ENV_NAME}"
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune --solver libmamba
else
  conda env create -n "${ENV_NAME}" -f "${ENV_FILE}" --solver libmamba
fi

echo "[2/5] Upgrading pip inside ${ENV_NAME}"
conda run -n "${ENV_NAME}" pip install --upgrade pip

echo "[3/5] Installing R package LncFinder inside ${ENV_NAME}"
conda run -n "${ENV_NAME}" R -q -e 'options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")); if (!requireNamespace("LncFinder", quietly = TRUE)) install.packages("LncFinder", dependencies = TRUE)'

echo "[4/5] Installing DeepRESIS in editable mode"
conda run -n "${ENV_NAME}" bash -lc "cd '${ROOT_DIR}' && pip install -e ."

echo "[5/5] Running runtime checks"
conda run -n "${ENV_NAME}" bash -lc "cd '${ROOT_DIR}' && python scripts/check_runtime.py"

cat <<EOF

Bootstrap finished.

Activate the environment:
  conda activate ${ENV_NAME}

Then run:
  cd ${ROOT_DIR}
  deepresis --help
EOF
