"""Build a protein-change index of pathogenic ClinVar missense variants (for PS1 and PM5).

Input: ClinVar's variant_summary.txt.gz. Each row's Name looks like
    NM_005228.5(EGFR):c.2573T>G (p.Leu858Arg)
from which we take the gene and protein change.
"""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

NAME_RE = re.compile(r"^\S+?\((?P<gene>[^)]+)\):c\.\S+ \((?P<protein>p\.[^)]+)\)")
MISSENSE_RE = re.compile(r"^p\.(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2})$")
PATHOGENIC_RE = re.compile(r"^(Pathogenic|Likely pathogenic|Pathogenic/Likely pathogenic)(;|$)")

OUTPUT_COLUMNS = ["gene", "ref_aa", "aa_pos", "alt_aa", "chrom", "pos", "ref", "alt",
                  "significance", "stars", "variation_id"]


def review_stars(status: str) -> int:
    """ClinVar review status -> gold stars (0-4). Works for both the VCF (underscore) and TSV forms."""
    s = status.replace("_", " ").replace("&", " ").lower()
    if "no assertion" in s or "no classification" in s:   # contains "criteria provided" too, so check first
        return 0
    if "practice guideline" in s:
        return 4
    if "expert panel" in s:
        return 3
    if "multiple submitters" in s and "no conflicts" in s:
        return 2
    if "criteria provided" in s:      # single submitter, or conflicting classifications
        return 1
    return 0


def parse_missense(protein: str) -> tuple[str, int, str] | None:
    """'p.Leu858Arg' -> ('Leu', 858, 'Arg'); None for anything that is not a simple missense change."""
    m = MISSENSE_RE.match(protein)
    if not m or m["ref"] == m["alt"] or m["alt"] == "Ter":
        return None
    return m["ref"], int(m["pos"]), m["alt"]


def build_index(variant_summary: str | Path, output: str | Path, min_stars: int = 1) -> int:
    """Write the index as a gzipped TSV and return the number of records."""
    seen: set[tuple] = set()
    n = 0
    with gzip.open(variant_summary, "rt", newline="") as fin, gzip.open(output, "wt", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        cols = reader.fieldnames or []
        # ClinVar renamed this column in 2024; support both names
        sig_col = next(c for c in ("ClinicalSignificance", "GermlineClassification") if c in cols)
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(OUTPUT_COLUMNS)
        for row in reader:
            if row["Assembly"] != "GRCh38" or row["Type"] != "single nucleotide variant":
                continue
            if not PATHOGENIC_RE.match(row[sig_col]):
                continue
            stars = review_stars(row["ReviewStatus"])
            if stars < min_stars:
                continue
            name = NAME_RE.match(row["Name"])
            if not name:
                continue
            aa = parse_missense(name["protein"])
            if not aa:
                continue
            chrom = row["Chromosome"]
            chrom = "chrM" if chrom == "MT" else f"chr{chrom}"
            key = (name["gene"], *aa, chrom, row["PositionVCF"], row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"])
            if key in seen:
                continue
            seen.add(key)
            writer.writerow([name["gene"], aa[0], aa[1], aa[2], chrom, row["PositionVCF"],
                             row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"],
                             row[sig_col], stars, row["VariationID"]])
            n += 1
    return n
