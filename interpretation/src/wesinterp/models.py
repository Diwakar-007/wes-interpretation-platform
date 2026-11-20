"""Core data structures: an annotated variant and a piece of ACMG evidence."""

from __future__ import annotations

from dataclasses import dataclass, field

# Points per evidence strength (Tavtigian et al. 2020, Hum Mutat 41:1734).
# Pathogenic evidence adds points, benign evidence subtracts them.
STRENGTH_POINTS = {"VeryStrong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1, "StandAlone": 0}

# The strength each criterion has by default in the original ACMG 2015 guideline.
# If a rule is applied at a different strength, its label gets a suffix, e.g. "PM2_Supporting".
DEFAULT_STRENGTH = {
    "PVS": "VeryStrong", "PS": "Strong", "PM": "Moderate", "PP": "Supporting",
    "BA": "StandAlone", "BS": "Strong", "BP": "Supporting",
}


@dataclass
class Variant:
    """One variant allele with the VEP consequence picked for one gene."""

    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str = "."
    gene: str = ""
    transcript: str = ""                    # MANE Select RefSeq ID when available
    biotype: str = "protein_coding"
    consequence: list[str] = field(default_factory=list)
    hgvsc: str = ""
    hgvsp: str = ""
    exon: str = ""                          # "5/12" = exon 5 of 12
    popmax_af: float | None = None          # highest continental gnomAD AF; None = absent from gnomAD
    gnomad_af: float | None = None          # overall gnomAD exomes AF
    revel: float | None = None
    am_class: str = ""
    am_score: float | None = None
    clinvar_sig: str = ""
    clinvar_review: str = ""
    clinvar_disease: str = ""
    rsid: str = ""

    @property
    def key(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


@dataclass(frozen=True)
class Evidence:
    """One ACMG criterion that was met, with its strength and a human-readable reason."""

    code: str           # e.g. "PVS1", "PM2", "BA1"
    strength: str       # VeryStrong / Strong / Moderate / Supporting / StandAlone
    reason: str

    @property
    def is_benign(self) -> bool:
        return self.code.startswith("B")

    @property
    def points(self) -> int:
        p = STRENGTH_POINTS[self.strength]
        return -p if self.is_benign else p

    @property
    def label(self) -> str:
        prefix = self.code.rstrip("0123456789")
        if DEFAULT_STRENGTH.get(prefix) == self.strength:
            return self.code
        return f"{self.code}_{self.strength}"
