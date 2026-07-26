from dataclasses import dataclass


@dataclass(frozen=True)
class DatasetConfig:
    """Source metadata required to build a matched ATAC/Hi-C dataset."""

    name: str
    cell_line: str
    genome_build: str
    hic_accession: str
    atac_accession: str
    hic_output_type: str = "mapping quality thresholded contact matrix"
    atac_output_type: str = "signal p-value"
    hic_normalisation: str = "VC_SQRT"
    resolution_bp: int = 100_000
    hic_md5: str = ""
    atac_md5: str = ""

    @property
    def hic_url(self) -> str:
        return (
            f"https://www.encodeproject.org/files/{self.hic_accession}/"
            f"@@download/{self.hic_accession}.hic"
        )

    @property
    def atac_url(self) -> str:
        return (
            f"https://www.encodeproject.org/files/{self.atac_accession}/"
            f"@@download/{self.atac_accession}.bigWig"
        )


DATASETS = {
    "gm12878": DatasetConfig(
        name="gm12878",
        cell_line="GM12878",
        genome_build="GRCh38",
        hic_accession="ENCFF916GWS",
        atac_accession="ENCFF667MDI",
        hic_md5="49db6d4b5b6a85bad3f42e935406b3b6",
        atac_md5="94e11903753580c8e39d77b331f90caf",
    ),
    "k562": DatasetConfig(
        name="k562",
        cell_line="K562",
        genome_build="GRCh38",
        hic_accession="ENCFF080DPJ",
        atac_accession="ENCFF600FDO",
        hic_md5="17a4977bb78e267f77199e915ec9243b",
        atac_md5="044f8568fb4c7735cf01e4f8d0d3e5b9",
    ),
}

AUTOSOMES = [f"chr{chromosome}" for chromosome in range(1, 23)]
