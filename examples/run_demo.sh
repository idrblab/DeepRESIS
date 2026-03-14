#!/usr/bin/env bash
set -euo pipefail

if [[ -f ./deepresis.toml ]]; then
  deepresis predict \
    --config ./deepresis.toml \
    --fasta /mnt/e/ncdresis/ncRESIS_Files/test/test.fasta \
    --smiles /mnt/e/ncdresis/ncRESIS_Files/test/test_cid.txt \
    --pairs /mnt/e/ncdresis/ncRESIS_Files/test/test_pairs.txt \
    --output-dir ./demo_outputs \
    --topk 300
else
  deepresis predict \
    --fasta /mnt/e/ncdresis/ncRESIS_Files/test/test.fasta \
    --smiles /mnt/e/ncdresis/ncRESIS_Files/test/test_cid.txt \
    --pairs /mnt/e/ncdresis/ncRESIS_Files/test/test_pairs.txt \
    --output-dir ./demo_outputs \
    --model-dir /mnt/e/ncdresis/ncRESIS_Files/model_pharameter \
    --gene-matrix-dir /mnt/e/ncdresis/ncRESIS_Files/gene_matrix \
    --topk 300
fi
