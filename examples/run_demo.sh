#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f ./deepresis.toml ]]; then
  echo "Missing ./deepresis.toml. Run 'bash scripts/bootstrap_all.sh' first." >&2
  exit 1
fi

deepresis predict \
  --config ./deepresis.toml \
  --fasta /mnt/e/ncdresis/ncRESIS_Files/test/test.fasta \
  --smiles /mnt/e/ncdresis/ncRESIS_Files/test/test_cid.txt \
  --pairs /mnt/e/ncdresis/ncRESIS_Files/test/test_pairs.txt \
  --output-dir ./demo_outputs \
  --topk 300
