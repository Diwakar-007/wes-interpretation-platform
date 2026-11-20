"""Read QC metrics from a pipeline sample output folder (results/<sample>/qc/...)."""

from __future__ import annotations

import gzip
from pathlib import Path


def _mosdepth(qc_dir: Path, sample: str) -> dict:
    out: dict = {}
    summary = qc_dir / "mosdepth" / f"{sample}.mosdepth.summary.txt"
    if summary.exists():
        for line in summary.read_text().splitlines():
            parts = line.split("\t")
            if parts[0] == "total_region":
                out["mean_coverage"] = float(parts[3])
    thresholds = qc_dir / "mosdepth" / f"{sample}.thresholds.bed.gz"
    if thresholds.exists():
        total = 0
        covered: dict[str, int] = {}
        with gzip.open(thresholds, "rt") as fh:
            header = fh.readline().lstrip("#").split()
            cols = header[4:]                          # e.g. 10X 20X 30X 50X
            for line in fh:
                p = line.split()
                total += int(p[2]) - int(p[1])
                for name, value in zip(cols, p[4:]):
                    covered[name] = covered.get(name, 0) + int(value)
        if total:
            out["target_bp"] = total
            out["pct_covered"] = {name: 100 * n / total for name, n in covered.items()}
    return out


def _duplication(qc_dir: Path, sample: str) -> dict:
    f = qc_dir / "markduplicates" / f"{sample}.md_metrics.txt"
    if not f.exists():
        return {}
    lines = f.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY"):
            row = dict(zip(line.split("\t"), lines[i + 1].split("\t")))
            return {"pct_duplication": 100 * float(row["PERCENT_DUPLICATION"]),
                    "read_pairs": int(row["READ_PAIRS_EXAMINED"])}
    return {}


def _bcftools(qc_dir: Path, sample: str) -> dict:
    f = qc_dir / "bcftools" / f"{sample}.bcftools_stats.txt"
    if not f.exists():
        return {}
    out: dict = {}
    for line in f.read_text().splitlines():
        p = line.split("\t")
        if p[0] == "SN" and len(p) >= 4:
            if p[2] == "number of SNPs:":
                out["snvs"] = int(p[3])
            elif p[2] == "number of indels:":
                out["indels"] = int(p[3])
        elif p[0] == "TSTV" and len(p) >= 5:
            out["ts_tv"] = float(p[4])
    return out


def read_qc(sample_dir: str | Path | None, sample: str) -> dict:
    if not sample_dir:
        return {}
    qc_dir = Path(sample_dir) / "qc"
    return {**_mosdepth(qc_dir, sample), **_duplication(qc_dir, sample), **_bcftools(qc_dir, sample)}
