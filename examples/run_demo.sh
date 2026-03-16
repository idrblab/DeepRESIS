#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f ./deepresis.toml ]]; then
  echo "Missing ./deepresis.toml. Run 'bash scripts/bootstrap_all.sh' first." >&2
  exit 1
fi

deepresis predict \
  --config ./deepresis.toml \
  --fasta ./examples/test.fasta \
  --smiles ./examples/test_cid.txt \
  --pairs ./examples/test_pairs.txt \
  --output-dir ./demo_outputs \
  --topk 300
