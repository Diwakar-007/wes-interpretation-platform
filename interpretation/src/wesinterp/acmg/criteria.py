"""ACMG/AMP criteria (Richards et al. 2015) with ClinGen SVI refinements.

Each rule is a function (variant, context) -> Evidence | None, so it can be unit-tested on its own.
Only criteria that can be assessed from annotation alone are automated; criteria that need
case data (PS2/PM6 de novo, PS4 case-control, PP1 segregation, PS3/BS3 functional) are
left for the analyst.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from ..clinvar_index import parse_missense
from ..config import DEFAULT_THRESHOLDS, Thresholds
from ..models import Evidence, Variant
from ..resources import ClinVarProteinIndex

LOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"}
INFRAME_CONSEQUENCES = {"inframe_insertion", "inframe_deletion", "stop_lost"}


@dataclass
class Context:
    """Resources shared by all rules."""

    loeuf: dict[str, float] = field(default_factory=dict)
    clinvar_protein: ClinVarProteinIndex = field(default_factory=dict)
    gene_modes: dict[str, set[str]] = field(default_factory=dict)   # e.g. {"CFTR": {"AR"}}
    thresholds: Thresholds = DEFAULT_THRESHOLDS

    def bs1_threshold(self, gene: str) -> tuple[float, str]:
        """BS1 allele-frequency cut-off for a gene, based on its known inheritance modes."""
        modes = self.gene_modes.get(gene, set())
        t = self.thresholds
        if modes and modes <= {"AR", "XLR"}:
            return t.bs1_af_recessive, "recessive"
        if modes and modes <= {"AD", "XLD"}:
            return t.bs1_af_dominant, "dominant"
        return t.bs1_af, "mixed/unknown inheritance" if modes else "inheritance unknown"


def _has(variant: Variant, consequences: set[str]) -> bool:
    return bool(set(variant.consequence) & consequences)


def _is_missense(variant: Variant) -> bool:
    return "missense_variant" in variant.consequence


def _is_last_exon(variant: Variant) -> bool:
    try:
        number, total = variant.exon.split("/")
        return int(number) == int(total)
    except ValueError:
        return False


# ---------------------------------------------------------------- pathogenic

def pvs1(v: Variant, ctx: Context) -> Evidence | None:
    """Null variant in a gene where loss of function causes disease.

    Simplified ClinGen PVS1 decision tree (Abou Tayoun 2018):
      - LoF mechanism: approximated by gnomAD LOEUF < 0.6
      - stop/frameshift in the LAST exon escapes nonsense-mediated decay -> Moderate
      - start_lost -> Moderate
    """
    is_lof = _has(v, LOF_CONSEQUENCES)
    if not (is_lof or "start_lost" in v.consequence):
        return None
    loeuf = ctx.loeuf.get(v.gene)
    if loeuf is None or loeuf >= ctx.thresholds.loeuf_max:
        return None
    note = f"LOEUF {loeuf:.2f}"
    if not is_lof:
        return Evidence("PVS1", "Moderate", f"start_lost in LoF-intolerant gene ({note})")
    if _has(v, {"stop_gained", "frameshift_variant"}) and _is_last_exon(v):
        return Evidence("PVS1", "Moderate", f"truncating variant in last exon, likely escapes NMD ({note})")
    return Evidence("PVS1", "VeryStrong", f"{'/'.join(sorted(set(v.consequence) & LOF_CONSEQUENCES))} "
                                          f"in LoF-intolerant gene ({note})")


def ps1_pm5(v: Variant, ctx: Context) -> Evidence | None:
    """PS1: same amino-acid change as a known pathogenic variant (different nucleotide change).
    PM5: different amino-acid change at a residue where another missense change is pathogenic.
    PS1 and PM5 are not applied together."""
    if not _is_missense(v):
        return None
    aa = parse_missense(v.hgvsp)
    if not aa:
        return None
    ref_aa, aa_pos, alt_aa = aa
    known = ctx.clinvar_protein.get(v.gene, {}).get(aa_pos, [])
    known = [r for r in known if r.ref_aa == ref_aa and r.genomic_key != v.key]
    same = [r for r in known if r.alt_aa == alt_aa]
    if same:
        r = max(same, key=lambda x: x.stars)
        return Evidence("PS1", "Strong", f"same change {r.protein_change} is {r.significance} in ClinVar "
                                         f"(VCV{r.variation_id}, {r.stars}*) via a different nucleotide change")
    other = [r for r in known if r.alt_aa != alt_aa]
    if other:
        r = max(other, key=lambda x: x.stars)
        return Evidence("PM5", "Moderate", f"different change at same residue, {r.protein_change}, "
                                           f"is {r.significance} in ClinVar (VCV{r.variation_id}, {r.stars}*)")
    return None


def pm2(v: Variant, ctx: Context) -> Evidence | None:
    """Absent or extremely rare in gnomAD. Supporting strength per ClinGen SVI (2020)."""
    if v.popmax_af is None:
        return Evidence("PM2", "Supporting", "absent from gnomAD")
    if v.popmax_af < ctx.thresholds.pm2_af:
        return Evidence("PM2", "Supporting", f"gnomAD popmax AF {v.popmax_af:.2e} < {ctx.thresholds.pm2_af}")
    return None


def pm4(v: Variant, ctx: Context) -> Evidence | None:
    """Protein length change from an in-frame indel or stop-loss (repeat regions not yet excluded)."""
    if _has(v, INFRAME_CONSEQUENCES) and not _has(v, LOF_CONSEQUENCES):
        kind = "/".join(sorted(set(v.consequence) & INFRAME_CONSEQUENCES))
        return Evidence("PM4", "Moderate", f"protein length change ({kind})")
    return None


def pp3_bp4(v: Variant, ctx: Context) -> Evidence | None:
    """Computational evidence from REVEL, calibrated by Pejaver et al. 2022. Missense only."""
    if not _is_missense(v) or v.revel is None:
        return None
    t, s = ctx.thresholds, v.revel
    for cut, strength in ((t.revel_pp3_strong, "Strong"), (t.revel_pp3_moderate, "Moderate"),
                          (t.revel_pp3_supporting, "Supporting")):
        if s >= cut:
            return Evidence("PP3", strength, f"REVEL {s:.3f} >= {cut}")
    for cut, strength in ((t.revel_bp4_strong, "Strong"), (t.revel_bp4_moderate, "Moderate"),
                          (t.revel_bp4_supporting, "Supporting")):
        if s <= cut:
            return Evidence("BP4", strength, f"REVEL {s:.3f} <= {cut}")
    return None


# ---------------------------------------------------------------- benign

def ba1(v: Variant, ctx: Context) -> Evidence | None:
    """Allele frequency > 5% in a continental population: stand-alone benign."""
    if v.popmax_af is not None and v.popmax_af > ctx.thresholds.ba1_af:
        return Evidence("BA1", "StandAlone", f"gnomAD popmax AF {v.popmax_af:.3f} > {ctx.thresholds.ba1_af}")
    return None


def bs1(v: Variant, ctx: Context) -> Evidence | None:
    """Allele frequency greater than expected for the disorder (inheritance-aware; not applied with BA1)."""
    cutoff, basis = ctx.bs1_threshold(v.gene)
    if v.popmax_af is not None and cutoff < v.popmax_af <= ctx.thresholds.ba1_af:
        return Evidence("BS1", "Strong", f"gnomAD popmax AF {v.popmax_af:.4f} > {cutoff} ({basis})")
    return None


def bp7(v: Variant, ctx: Context) -> Evidence | None:
    """Synonymous change outside the splice region (no splice predictor yet, so splice-region excluded)."""
    if "synonymous_variant" in v.consequence and "splice_region_variant" not in v.consequence:
        return Evidence("BP7", "Supporting", "synonymous, outside splice region")
    return None


ALL_RULES: list[Callable[[Variant, Context], Evidence | None]] = [
    pvs1, ps1_pm5, pm2, pm4, pp3_bp4, ba1, bs1, bp7,
]
