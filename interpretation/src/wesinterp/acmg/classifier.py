"""Combine ACMG evidence into a 5-tier classification using the Bayesian points system
(Tavtigian et al. 2020): Pathogenic >= 10, Likely pathogenic 6..9, VUS 0..5,
Likely benign -1..-6, Benign <= -7. BA1 is stand-alone benign.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..clinvar_index import review_stars
from ..config import BENIGN, CLASS_POINTS
from ..models import Evidence, Variant
from .criteria import ALL_RULES, Context


@dataclass
class Classification:
    variant: Variant
    evidence: list[Evidence] = field(default_factory=list)
    points: int = 0
    classification: str = "Uncertain significance"
    clinvar_concordance: str = ""     # "", "concordant", "discordant"

    @property
    def criteria(self) -> str:
        return ";".join(e.label for e in self.evidence)


def points_to_class(points: int, evidence: list[Evidence]) -> str:
    if any(e.code == "BA1" for e in evidence):
        return BENIGN
    for minimum, label in CLASS_POINTS:
        if points >= minimum:
            return label
    return BENIGN


def _tier(label: str) -> str:
    """Collapse a classification (ours or ClinVar's) to P / VUS / B."""
    s = label.lower().replace("_", " ")
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s:
        return "P"
    if "benign" in s:
        return "B"
    if "uncertain" in s:
        return "VUS"
    return ""


def clinvar_concordance(variant: Variant, our_class: str) -> str:
    """Compare with ClinVar when ClinVar's entry is well reviewed (>= 2 stars).
    A 'discordant' result means the variant needs manual review."""
    if not variant.clinvar_sig or review_stars(variant.clinvar_review) < 2:
        return ""
    theirs = _tier(variant.clinvar_sig)
    if theirs in ("", "conflicting"):
        return ""
    return "concordant" if theirs == _tier(our_class) else "discordant"


def classify_variant(variant: Variant, ctx: Context) -> Classification:
    evidence = [e for rule in ALL_RULES if (e := rule(variant, ctx)) is not None]
    points = sum(e.points for e in evidence)
    label = points_to_class(points, evidence)
    return Classification(variant, evidence, points, label, clinvar_concordance(variant, label))
