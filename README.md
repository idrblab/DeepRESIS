# DeepRESIS

DeepRESIS is a command-line tool for predicting ncRNA-drug resistance relationships. You provide an ncRNA FASTA file, a drug SMILES file, and an ncRNA-drug pair file. The tool produces resistance prediction results and top-k gene ranking outputs.

## Quick Start

On a fresh machine:

```bash
git clone https://github.com/idrblab/DeepRESIS.git
cd DeepRESIS
bash scripts/bootstrap_all.sh
conda activate deepresis
deepresis --help
```

This is the main installation path. It installs the environment, downloads the required large files from Hugging Face, generates `deepresis.toml`, and validates that the runtime is ready.

## What the Bootstrap Script Does

`bash scripts/bootstrap_all.sh` will automatically:

1. create a conda environment
2. install Python dependencies
3. install the R runtime and required R packages
4. install the R package `LncFinder`
5. install DeepRESIS in editable mode
6. download 5 model checkpoint files from Hugging Face
7. download `drug_gene.csv` and `ncrna_gene.csv` from Hugging Face
8. generate `deepresis.toml`
9. validate the runtime and downloaded assets

Large model and gene matrix files are not stored in the GitHub repository. They are fetched automatically during bootstrap.

## Files You Need

DeepRESIS needs three input files.

### 1. ncRNA FASTA file

Example:

```text
>ncRNA_1
UGGGGGAGUGAAGAGUAGAUAAAAUAUUGGUACCUGAUGAAUCUGAGGCCAGGUUUCAAU
ACUUUAUCUGCUCUUCAUUUCCCCAUAUCUACUUAC
>ncRNA_2
UAGCACCAUCUGAAAUCGGUUA
```

The FASTA record ID becomes the `ncrna_id`.

### 2. Drug SMILES file

Example:

```text
drug1    CCO
drug2    CCN(CC)CC
```

The first column is `drug_id`, and the second column is the SMILES string.

### 3. ncRNA-drug pair file

Example:

```text
ncRNA_1    drug1
ncRNA_2    drug2
```

The first column is `ncrna_id`, and the second column is `drug_id`.

Important: every `ncrna_id` in the pair file must exist in the FASTA file, and every `drug_id` in the pair file must exist in the SMILES file.

## Run Your First Prediction

After bootstrap finishes, run prediction with the generated config file:

```bash
deepresis predict \
  --config ./deepresis.toml \
  --fasta /path/to/test.fasta \
  --smiles /path/to/test_cid.txt \
  --pairs /path/to/test_pairs.txt \
  --output-dir ./outputs \
  --topk 300
```

You do not need to pass `--model-dir` or `--gene-matrix-dir` in the normal workflow.

## Output Files

DeepRESIS writes three output files in the output directory.

### 1. `resistance_predictions.tsv`

This file contains resistance prediction results for each requested ncRNA-drug pair.

Columns:

- `ncrna_id`
- `drug_id`
- `score`
- `label`
- `label_id`

### 2. `topk_genes.tsv`

This file contains the top-k ranked genes for each pair when gene ranking is available.

Columns:

- `ncrna_id`
- `drug_id`
- `rank`
- `gene_id`
- `score`
- `drug_gene_score`
- `ncrna_gene_score`

### 3. `gene_ranking_warnings.tsv`

This file records why gene ranking could not be generated for specific pairs.

Columns:

- `ncrna_id`
- `drug_id`
- `message`

## Demo Example

You can run the repository demo after installation:

```bash
bash examples/run_demo.sh
```

Or run the same example explicitly:

```bash
deepresis predict \
  --config ./deepresis.toml \
  --fasta ./examples/test.fasta \
  --smiles ./examples/test_cid.txt \
  --pairs ./examples/test_pairs.txt \
  --output-dir ./demo_outputs
```

Expected behavior for the current sample:

- `resistance_predictions.tsv` contains exactly 3 requested pairs
- `topk_genes.tsv` may contain only the header
- `gene_ranking_warnings.tsv` may contain 3 warnings because the sample circRNA IDs are not present in `ncrna_gene.csv`

## Troubleshooting

### `conda: command not found`

Conda is not installed or not on your `PATH`.

Next step: install Miniconda or Anaconda first, then rerun `bash scripts/bootstrap_all.sh`.

### Hugging Face download failed

The machine cannot access Hugging Face, or the connection was interrupted.

Next step: check internet access and rerun the bootstrap command.

### Missing model checkpoint files

One or more downloaded checkpoint files are missing under `artifacts/models`.

Next step: delete the incomplete `artifacts` directory and rerun `bash scripts/bootstrap_all.sh`.

### Missing gene matrix files

`drug_gene.csv` or `ncrna_gene.csv` is missing under `artifacts/gene_matrix`.

Next step: delete the incomplete `artifacts` directory and rerun `bash scripts/bootstrap_all.sh`.

### `RNAfold` not found

The ViennaRNA binary is not available in the installed environment.

Next step: rerun the bootstrap command and confirm that environment creation completed successfully.

### `LncFinder` check failed

The R package `LncFinder` was not installed correctly.

Next step: rerun the bootstrap command. If it still fails, inspect the R installation step printed by the script.

## Advanced Options

### Custom environment name

```bash
bash scripts/bootstrap_all.sh deepresis_prod
```

### Custom environment name, asset directory, and config path

```bash
bash scripts/bootstrap_all.sh deepresis_prod /data/deepresis_assets /data/deepresis.toml
```

### Override Hugging Face source URLs

Normally you do not need this. If needed, you can override the Hugging Face directory URLs:

```bash
export DEEPRESIS_MODELS_URL="https://huggingface.co/swallow-design/DeepRESIS/tree/main/model_pharameter"
export DEEPRESIS_GENE_MATRIX_URL="https://huggingface.co/swallow-design/DeepRESIS/tree/main/gene_matrix"
bash scripts/bootstrap_all.sh
```

The bootstrap script will automatically normalize `tree/main` or `blob/main` URLs to real `resolve/main` download URLs.

### Python API

```python
from deepresis.pipeline import run_prediction

run_prediction(
    fasta_path="test.fasta",
    smiles_path="test_cid.txt",
    pairs_path="test_pairs.txt",
    output_dir="outputs",
    config_path="deepresis.toml",
    topk=300,
)
```
