"""Rank candidate genes by combining ACMG class, phenotype match and inheritance fit.

    gene score = phenotype score (0..1) x best variant class weight x inheritance fit weight

Kept separate from ACMG on purpose: phenotype match is used to *prioritise* which variants an
analyst looks at first, not as hidden extra evidence inside the classification.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from .acmg.classifier import Classification
from .hpo import GeneAnnotations, Ontology, phenotype_score
from .hpo.annotations import DOMINANT_MODES, RECESSIVE_MODES

CLASS_WEIGHT = {"Pathogenic": 1.0, "Likely pathogenic": 0.9, "Uncertain significance": 0.5,
                "Likely benign": 0.1, "Benign": 0.0}

# Consequences that can plausibly cause Mendelian disease (others are not ranked)
RELEVANT_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
    "start_lost", "stop_lost", "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "splice_region_variant", "splice_donor_5th_base_variant",
}


@dataclass
class Candidate:
    gene: str
    modes: list[str]
    phenotype_score: float
    matched_terms: list[tuple[str, str, float]]
    variants: list[Classification] = field(default_factory=list)
    inheritance_fit: str = ""
    fit_weight: float = 1.0
    score: float = 0.0


def zygosity(genotype: str) -> str:
    alleles = genotype.replace("|", "/").split("/")
    alt = [a for a in alleles if a not in ("0", ".")]
    if len(alleles) == 1:
        return "hemizygous" if alt else "ref"
    if len(alt) == 2:
        return "homozygous"
    return "heterozygous" if alt else "ref"


def inheritance_fit(modes: set[str], variants: list[Classification]) -> tuple[str, float]:
    """Does the number/zygosity of qualifying variants fit the gene's inheritance mode?"""
    zyg = [zygosity(r.variant.genotype) for r in variants]
    hom = any(z in ("homozygous", "hemizygous") for z in zyg)
    n_het = zyg.count("heterozygous")
    if not modes:
        return "inheritance unknown", 0.7
    if modes & (DOMINANT_MODES | {"XL"}):
        return f"fits {'/'.join(sorted(modes & (DOMINANT_MODES | {'XL'})))}", 1.0
    if modes & RECESSIVE_MODES:
        if hom:
            return "fits AR (homozygous)", 1.0
        if n_het >= 2:
            return "possible AR (compound heterozygous, phase unknown)", 0.8
        return "AR gene but single heterozygous variant (carrier?)", 0.3
    return f"modes: {'/'.join(sorted(modes))}", 0.7


def is_relevant(r: Classification) -> bool:
    v = r.variant
    return (v.biotype == "protein_coding"
            and bool(set(v.consequence) & RELEVANT_CONSEQUENCES)
            and r.classification not in ("Benign", "Likely benign"))


def prioritize(results: list[Classification], patient_terms: list[str],
               annotations: GeneAnnotations, ontology: Ontology) -> list[Candidate]:
    by_gene: dict[str, list[Classification]] = defaultdict(list)
    for r in results:
        if is_relevant(r):
            by_gene[r.variant.gene].append(r)

    candidates = []
    for gene, variants in by_gene.items():
        modes = annotations.gene_modes.get(gene, set())
        if patient_terms:
            pheno, matches = phenotype_score(patient_terms, gene, annotations, ontology)
        else:
            pheno, matches = 1.0, []
        fit, fit_w = inheritance_fit(modes, variants)
        best = max(CLASS_WEIGHT[r.classification] for r in variants)
        variants.sort(key=lambda r: -CLASS_WEIGHT[r.classification])
        candidates.append(Candidate(gene, sorted(modes), pheno, matches, variants, fit, fit_w,
                                    round(pheno * best * fit_w, 4)))
    candidates.sort(key=lambda c: -c.score)
    return candidates
