"""Phenotype similarity between a patient's HPO terms and a gene's known phenotypes.

Method: Resnik similarity with a one-sided best-match average (patient -> gene).
  sim(q, d)  = IC of the most informative common ancestor of terms q and d
  score      = mean over patient terms q of  max over gene terms d of sim(q, d)
  normalised = score / mean IC(q)          -> 0..1 (1 = every patient term is explained exactly)

Trick: max_d IC(MICA(q, d)) = max IC over (ancestors(q) ∩ union of ancestors of all d),
so one set intersection per patient term is enough.
"""

from __future__ import annotations

from .annotations import GeneAnnotations
from .ontology import Ontology


def phenotype_score(patient_terms: list[str], gene: str, annotations: GeneAnnotations,
                    ontology: Ontology) -> tuple[float, list[tuple[str, str, float]]]:
    """Return (normalised score 0..1, per-term matches [(patient_term, best_shared_term, IC)])."""
    gene_set = annotations.gene_ancestor_set(gene, ontology)
    ic = annotations.information_content
    if not patient_terms or not gene_set:
        return 0.0, []
    matches = []
    total = maximum = 0.0
    for q in patient_terms:
        q_ic = ic.get(q, 0.0)
        shared = ontology.ancestors(q) & gene_set
        best = max(shared, key=lambda t: ic.get(t, 0.0), default=None)
        best_ic = ic.get(best, 0.0) if best else 0.0
        matches.append((q, best or "", best_ic))
        total += best_ic
        maximum += q_ic
    return (total / maximum if maximum else 0.0), matches
