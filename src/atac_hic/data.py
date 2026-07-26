from __future__ import annotations

import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .config import DatasetConfig


def file_md5(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Return an MD5 digest for source-file provenance checks."""

    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: str | Path) -> Path:
    """Download a file only when it is not already present."""

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        urllib.request.urlretrieve(url, destination)
    return destination


def validate_processed_dataset(
    data: np.lib.npyio.NpzFile,
    chromosomes: list[str],
    symmetry_tolerance: float = 1e-5,
) -> list[dict]:
    """Validate paired ATAC vectors and symmetric Hi-C matrices."""

    rows = []
    for chromosome in chromosomes:
        atac = data[f"{chromosome}_ATAC"]
        hic = data[f"{chromosome}_HiC"]

        if hic.ndim != 2 or hic.shape[0] != hic.shape[1]:
            raise ValueError(f"{chromosome}: Hi-C matrix is not square")
        if len(atac) != hic.shape[0]:
            raise ValueError(f"{chromosome}: ATAC/Hi-C dimensions do not match")
        if not np.isfinite(atac).all() or not np.isfinite(hic).all():
            raise ValueError(f"{chromosome}: non-finite values detected")
        if (atac < 0).any() or (hic < 0).any():
            raise ValueError(f"{chromosome}: negative values detected")

        symmetry_error = float(np.max(np.abs(hic - hic.T)))
        if symmetry_error > symmetry_tolerance:
            raise ValueError(
                f"{chromosome}: symmetry error {symmetry_error} exceeds "
                f"{symmetry_tolerance}"
            )

        rows.append({
            "chromosome": chromosome,
            "bins": len(atac),
            "nonzero_hic_percent": float(np.mean(hic > 0) * 100),
            "symmetry_error": symmetry_error,
        })
    return rows


def extract_encode_dataset(
    config: DatasetConfig,
    chromosomes: list[str],
    output_path: str | Path,
    source_directory: str | Path,
) -> Path:
    """
    Extract matched chromosome-wide ATAC and Hi-C arrays from ENCODE.

    `hic-straw` and `pyBigWig` are imported lazily so the core analysis package
    can be tested without network/data-access dependencies.
    """

    import hicstraw
    import pyBigWig

    output_path = Path(output_path)
    source_directory = Path(source_directory)
    source_directory.mkdir(parents=True, exist_ok=True)

    atac_path = download_file(
        config.atac_url,
        source_directory / f"{config.atac_accession}.bigWig",
    )
    observed_atac_md5 = file_md5(atac_path)
    if config.atac_md5 and observed_atac_md5 != config.atac_md5:
        raise ValueError(
            f"{config.atac_accession}: expected MD5 {config.atac_md5}, "
            f"observed {observed_atac_md5}"
        )

    hic = hicstraw.HiCFile(config.hic_url)
    chromosome_lengths = {
        chromosome.name: chromosome.length
        for chromosome in hic.getChromosomes()
    }
    missing = sorted(set(chromosomes) - set(chromosome_lengths))
    if missing:
        raise ValueError(f"Chromosomes absent from Hi-C file: {missing}")

    atac_file = pyBigWig.open(str(atac_path))
    processed: dict[str, np.ndarray] = {}

    try:
        for chromosome in chromosomes:
            length = chromosome_lengths[chromosome]
            hic_data = hic.getMatrixZoomData(
                chromosome,
                chromosome,
                "observed",
                config.hic_normalisation,
                "BP",
                config.resolution_bp,
            )
            hic_matrix = np.asarray(
                hic_data.getRecordsAsMatrix(0, length, 0, length),
                dtype=np.float32,
            )
            hic_matrix = np.nan_to_num(
                hic_matrix,
                nan=0.0,
                posinf=0.0,
                neginf=0.0,
            )

            atac_values = atac_file.stats(
                chromosome,
                0,
                length,
                nBins=hic_matrix.shape[0],
                type="mean",
            )
            atac_values = np.asarray(
                [0.0 if value is None else value for value in atac_values],
                dtype=np.float32,
            )
            atac_values = np.nan_to_num(atac_values)

            processed[f"{chromosome}_ATAC"] = atac_values
            processed[f"{chromosome}_HiC"] = hic_matrix
            print(
                f"{config.cell_line} {chromosome}: "
                f"ATAC {atac_values.shape}, Hi-C {hic_matrix.shape}"
            )
    finally:
        atac_file.close()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **processed)

    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "cell_line": config.cell_line,
        "genome_build": config.genome_build,
        "hic_accession": config.hic_accession,
        "hic_url": config.hic_url,
        "hic_expected_md5": config.hic_md5,
        "atac_accession": config.atac_accession,
        "atac_url": config.atac_url,
        "atac_expected_md5": config.atac_md5,
        "hic_output_type": config.hic_output_type,
        "atac_output_type": config.atac_output_type,
        "hic_normalisation": config.hic_normalisation,
        "resolution_bp": config.resolution_bp,
        "chromosomes": chromosomes,
        "atac_source_md5": observed_atac_md5,
        "processed_npz_md5": file_md5(output_path),
    }
    with open(output_path.with_suffix(".provenance.json"), "w") as handle:
        json.dump(provenance, handle, indent=2)

    return output_path
