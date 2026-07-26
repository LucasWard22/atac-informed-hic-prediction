from __future__ import annotations

import json
import platform
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .data import file_md5

TRACKED_PACKAGES = (
    "hic-straw",
    "joblib",
    "matplotlib",
    "numpy",
    "pandas",
    "pyBigWig",
    "scikit-learn",
    "scipy",
)


def runtime_versions() -> dict[str, str]:
    packages = {}
    for package in TRACKED_PACKAGES:
        try:
            packages[package] = version(package)
        except PackageNotFoundError:
            packages[package] = "not installed"
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        **packages,
    }


def write_run_metadata(
    output_path: str | Path,
    analysis: str,
    input_paths: dict[str, str | Path],
    parameters: dict,
) -> Path:
    """Record input checksums, runtime versions and analysis parameters."""

    output_path = Path(output_path)
    payload = {
        "analysis": analysis,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            name: {
                "path": str(Path(path).resolve()),
                "md5": file_md5(path),
            }
            for name, path in input_paths.items()
        },
        "parameters": parameters,
        "runtime": runtime_versions(),
    }
    with open(output_path, "w") as handle:
        json.dump(payload, handle, indent=2)
    return output_path
