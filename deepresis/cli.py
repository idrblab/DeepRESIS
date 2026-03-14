from __future__ import annotations

import argparse

from .pipeline import run_prediction


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="deepresis", description="DeepRESIS ncRNA-drug prediction CLI.")
    subparsers = parser.add_subparsers(dest="command")
    predict_parser = subparsers.add_parser("predict", help="Run resistance prediction.")
    _add_shared_arguments(parser)
    _add_shared_arguments(predict_parser)
    return parser


def _add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--fasta", help="Path to the ncRNA FASTA file.")
    parser.add_argument("--smiles", help="Path to the drug SMILES file.")
    parser.add_argument("--pairs", help="Path to the ncRNA-drug pairs file.")
    parser.add_argument("--output-dir", help="Directory for output files.")
    parser.add_argument("--topk", type=int, default=None, help="Top-k genes to output.")
    parser.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda"],
        default=None,
        help="Inference device. Defaults to auto.",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed.")
    parser.add_argument("--config", help="Path to deepresis TOML config.")
    parser.add_argument("--model-dir", help="Directory containing fold checkpoints.")
    parser.add_argument("--gene-matrix-dir", help="Directory containing gene matrix CSV files.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None and not all([args.fasta, args.smiles, args.pairs, args.output_dir]):
        parser.print_help()
        return 0

    if args.command == "predict" or args.command is None:
        if not all([args.fasta, args.smiles, args.pairs, args.output_dir]):
            parser.error("--fasta, --smiles, --pairs, and --output-dir are required.")
        run_prediction(
            fasta_path=args.fasta,
            smiles_path=args.smiles,
            pairs_path=args.pairs,
            output_dir=args.output_dir,
            topk=args.topk,
            device=args.device,
            seed=args.seed,
            model_dir=args.model_dir,
            gene_matrix_dir=args.gene_matrix_dir,
            config_path=args.config,
            verbose=args.verbose,
        )
        return 0

    parser.error(f"Unsupported command: {args.command}")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
