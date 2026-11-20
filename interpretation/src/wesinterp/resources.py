"""Loaders for gene- and protein-level resources used by the ACMG rules."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


def load_loeuf(path: str | Path) -> dict[str, float]:
    """Gene symbol -> LOEUF from the gnomAD v4.1 constraint table.

    Uses the MANE Select transcript row for each gene, falling back to the canonical one.
    """
    mane: dict[str, float] = {}
    canonical: dict[str, float] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "lof.oe_ci.upper" not in (reader.fieldnames or []):
            raise ValueError(f"{path}: no 'lof.oe_ci.upper' (LOEUF) column")
        for row in reader:
            value = row["lof.oe_ci.upper"]
            if value in ("", "NA", "nan"):
                continue
            if row.get("mane_select") == "true":
                mane[row["gene"]] = float(value)
            elif row.get("canonical") == "true":
                canonical.setdefault(row["gene"], float(value))
    return {**canonical, **mane}


@dataclass(frozen=True)
class ClinVarProteinRecord:
    """A pathogenic/likely pathogenic missense variant in ClinVar, keyed by its protein change."""

    gene: str
    ref_aa: str
    aa_pos: int
    alt_aa: str
    chrom: str
    pos: int
    ref: str
    alt: str
    significance: str
    stars: int
    variation_id: str

    @property
    def protein_change(self) -> str:
        return f"p.{self.ref_aa}{self.aa_pos}{self.alt_aa}"

    @property
    def genomic_key(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


# gene -> amino-acid position -> records at that residue
ClinVarProteinIndex = dict[str, dict[int, list[ClinVarProteinRecord]]]


def load_clinvar_protein_index(path: str | Path) -> ClinVarProteinIndex:
    """Load the index written by `wesinterp build-clinvar-index`."""
    index: ClinVarProteinIndex = defaultdict(lambda: defaultdict(list))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rec = ClinVarProteinRecord(
                gene=row["gene"], ref_aa=row["ref_aa"], aa_pos=int(row["aa_pos"]), alt_aa=row["alt_aa"],
                chrom=row["chrom"], pos=int(row["pos"]), ref=row["ref"], alt=row["alt"],
                significance=row["significance"], stars=int(row["stars"]), variation_id=row["variation_id"],
            )
            index[rec.gene][rec.aa_pos].append(rec)
    return index
