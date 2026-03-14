# DeepRESIS

DeepRESIS is a lightweight Python package that wraps the existing `ncRESIS_Files` inference logic into an installable package and CLI.

## Repository Design

This GitHub-targeted repository is intentionally lightweight.

- Included: Python package code, vendored feature backend, tests, docs, examples
- Not included: model checkpoints and large gene-matrix CSV files

Those large assets should be downloaded separately and configured with a TOML config file or CLI flags.

## Large Files

- `TODO: Hugging Face model link`
- `TODO: Hugging Face gene matrix link`

Expected local layout after download:

```text
/data/deepresis/
├── models/
│   ├── fold0/ne_student.ckpt
│   ├── fold1/ne_student.ckpt
│   ├── fold2/ne_student.ckpt
│   ├── fold3/ne_student.ckpt
│   └── fold4/ne_student.ckpt
└── gene_matrix/
    ├── drug_gene.csv
    └── ncrna_gene.csv
```

### Asset Bootstrap

There are now two asset bootstrap modes:

1. local copy mode for development and validation
2. URL download mode for future Hugging Face archives

Local copy mode:

```bash
export DEEPRESIS_LOCAL_MODEL_DIR=/path/to/model_pharameter
export DEEPRESIS_LOCAL_GENE_MATRIX_DIR=/path/to/gene_matrix
bash scripts/bootstrap_assets.sh
```

URL download mode:

```bash
export DEEPRESIS_MODELS_URL="https://..."
export DEEPRESIS_GENE_MATRIX_URL="https://..."
bash scripts/bootstrap_assets.sh
```

By default this creates:

- `./artifacts/models`
- `./artifacts/gene_matrix`
- `./deepresis.toml`

### Full Bootstrap

If you want environment setup and asset setup in one chain:

```bash
bash scripts/bootstrap_all.sh
```

## Installation

### One-command Bootstrap

On a fresh server, the intended setup path is:

```bash
git clone <your-github-repo> DeepRESIS
cd DeepRESIS
bash scripts/bootstrap_env.sh
```

This script will:

1. create or update a conda environment from `environment.yml`
2. install the R runtime and conda dependencies
3. install the R package `LncFinder`
4. run `pip install -e .`
5. verify that Python modules, `RNAfold`, and `LncFinder` are available

This bootstrap path has been validated in a fresh conda environment on the current machine.

Default environment name: `deepresis`

You can override it:

```bash
bash scripts/bootstrap_env.sh deepresis_prod
```

### Manual Equivalent

```bash
conda env create -n deepresis -f environment.yml
conda run -n deepresis pip install --upgrade pip
conda run -n deepresis R -q -e 'options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")); if (!requireNamespace("LncFinder", quietly = TRUE)) install.packages("LncFinder", dependencies = TRUE)'
conda run -n deepresis bash -lc 'cd /path/to/DeepRESIS && pip install -e .'
conda run -n deepresis bash -lc 'cd /path/to/DeepRESIS && python scripts/check_runtime.py'
```

If you also want assets and a config file prepared automatically, run:

```bash
bash scripts/bootstrap_assets.sh
```

## Configuration

Example `deepresis.toml`:

```toml
[paths]
model_dir = "/data/deepresis/models"
gene_matrix_dir = "/data/deepresis/gene_matrix"

[run]
device = "auto"
seed = 123
topk = 300
```

Resolution priority:

1. CLI flags
2. Environment variables
3. `--config`
4. `./deepresis.toml`
5. `~/.config/deepresis/config.toml`

Supported environment variables:

- `DEEPRESIS_CONFIG`
- `DEEPRESIS_MODEL_DIR`
- `DEEPRESIS_GENE_MATRIX_DIR`

## Input Formats

### FASTA

```text
>ncRNA_1
AUGCUAGCUAGCUA
>ncRNA_2
GGCUAAGCUU
```

- record ID becomes `ncrna_id`
- sequence must be non-empty
- allowed characters: `A/C/G/U/T/N`

### SMILES

```text
drug1    CCO
drug2    CCN(CC)CC
```

- first column: `drug_id`
- second column: `smiles`
- whitespace-separated input is accepted

### Pairs

```text
ncRNA_1    drug1
ncRNA_2    drug2
```

- first column: `ncrna_id`
- second column: `drug_id`
- pairs missing from FASTA or SMILES are skipped with warnings

## CLI Usage

```bash
deepresis predict \
  --fasta /mnt/e/ncdresis/ncRESIS_Files/test/test.fasta \
  --smiles /mnt/e/ncdresis/ncRESIS_Files/test/test_cid.txt \
  --pairs /mnt/e/ncdresis/ncRESIS_Files/test/test_pairs.txt \
  --output-dir ./outputs \
  --model-dir /mnt/e/ncdresis/ncRESIS_Files/model_pharameter \
  --gene-matrix-dir /mnt/e/ncdresis/ncRESIS_Files/gene_matrix \
  --topk 300
```

## Python Usage

```python
from deepresis.pipeline import run_prediction

run_prediction(
    fasta_path="test.fasta",
    smiles_path="test_cid.txt",
    pairs_path="test_pairs.txt",
    output_dir="outputs",
    model_dir="/data/deepresis/models",
    gene_matrix_dir="/data/deepresis/gene_matrix",
    topk=300,
)
```

## Outputs

### `resistance_predictions.tsv`

- `ncrna_id`
- `drug_id`
- `score`
- `label`
- `label_id`

### `topk_genes.tsv`

- `ncrna_id`
- `drug_id`
- `rank`
- `gene_id`
- `score`
- `drug_gene_score`
- `ncrna_gene_score`

### `gene_ranking_warnings.tsv`

- `ncrna_id`
- `drug_id`
- `message`

## Minimal Example

The repo includes lightweight tests that stub out heavy feature/model work. Use them to confirm the editable install and package wiring:

```bash
cd /mnt/e/ncdresis/DeepRESIS
conda run -n deepresis pytest tests/test_minimal_example.py
```

## Real Example

Use your local large assets:

```bash
cd /mnt/e/ncdresis/DeepRESIS
conda run -n deepresis deepresis predict \
  --fasta /mnt/e/ncdresis/ncRESIS_Files/test/test.fasta \
  --smiles /mnt/e/ncdresis/ncRESIS_Files/test/test_cid.txt \
  --pairs /mnt/e/ncdresis/ncRESIS_Files/test/test_pairs.txt \
  --output-dir ./outputs_real \
  --model-dir /mnt/e/ncdresis/ncRESIS_Files/model_pharameter \
  --gene-matrix-dir /mnt/e/ncdresis/ncRESIS_Files/gene_matrix \
  --topk 300
```

Or, if you already ran `bootstrap_assets.sh`, you can use the generated config file:

```bash
cd /mnt/e/ncdresis/DeepRESIS
conda run -n deepresis deepresis predict \
  --config ./deepresis.toml \
  --fasta /mnt/e/ncdresis/ncRESIS_Files/test/test.fasta \
  --smiles /mnt/e/ncdresis/ncRESIS_Files/test/test_cid.txt \
  --pairs /mnt/e/ncdresis/ncRESIS_Files/test/test_pairs.txt \
  --output-dir ./outputs_real
```

Expected behavior for the current sample:

- `resistance_predictions.tsv` should contain exactly the 3 requested pairs
- `topk_genes.tsv` may be empty except for the header
- `gene_ranking_warnings.tsv` may contain 3 warnings because the sample circRNA IDs are not present in `ncrna_gene.csv`

## Common Errors

- `Model directory is not configured`
- `Gene matrix directory is not configured`
- `Missing model checkpoint files`
- `Missing gene matrix files`
- `CUDA was requested but is not available`
