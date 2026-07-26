#!/usr/bin/env python
from argparse import ArgumentParser
from pathlib import Path

from atac_hic.config import AUTOSOMES, DATASETS
from atac_hic.data import extract_encode_dataset


def parse_args():
    parser = ArgumentParser(
        description="Extract matched 100-kb ATAC/Hi-C arrays from ENCODE."
    )
    parser.add_argument(
        "--dataset",
        choices=sorted(DATASETS),
        required=True,
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--source-directory",
        type=Path,
        default=Path("data/source"),
    )
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=AUTOSOMES,
        help="Chromosome names; defaults to chr1–chr22.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output = extract_encode_dataset(
        config=DATASETS[args.dataset],
        chromosomes=args.chromosomes,
        output_path=args.output,
        source_directory=args.source_directory / args.dataset,
    )
    print("Saved:", output)


if __name__ == "__main__":
    main()
