"""Gene -> phenotype terms, gene -> inheritance modes, and term information content."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .ontology import Ontology

# HPO inheritance terms -> short codes
INHERITANCE_CODES = {
    "HP:0000006": "AD",   # Autosomal dominant
    "HP:0000007": "AR",   # Autosomal recessive
    "HP:0001417": "XL",   # X-linked
    "HP:0001419": "XLR",  # X-linked recessive
    "HP:0001423": "XLD",  # X-linked dominant
    "HP:0001427": "MT",   # Mitochondrial
    "HP:0001450": "YL",   # Y-linked
}
DOMINANT_MODES = {"AD", "XLD"}
RECESSIVE_MODES = {"AR", "XLR"}


@dataclass
class GeneAnnotations:
    gene_terms: dict[str, set[str]] = field(default_factory=dict)        # phenotype terms per gene
    gene_modes: dict[str, set[str]] = field(default_factory=dict)        # inheritance codes per gene
    information_content: dict[str, float] = field(default_factory=dict)  # IC per term

    def gene_ancestor_set(self, gene: str, ontology: Ontology) -> frozenset[str]:
        """All terms implied by the gene's annotations (annotations + their ancestors)."""
        result: set[str] = set()
        for t in self.gene_terms.get(gene, ()):
            result |= ontology.ancestors(t)
        return frozenset(result)


def _association_types(genes_to_disease: str | Path) -> dict[tuple[str, str], str]:
    """(gene, disease) -> MENDELIAN / POLYGENIC / UNKNOWN"""
    with open(genes_to_disease, newline="") as fh:
        return {(r["gene_symbol"], r["disease_id"]): r["association_type"]
                for r in csv.DictReader(fh, delimiter="\t")}


def load_gene_annotations(genes_to_phenotype: str | Path, ontology: Ontology,
                          genes_to_disease: str | Path | None = None) -> GeneAnnotations:
    """If genes_to_disease is given, POLYGENIC (risk-factor) associations are ignored, and
    inheritance modes are taken only from MENDELIAN associations. Without this, CFTR would be
    labelled 'autosomal dominant' because of hereditary-pancreatitis susceptibility."""
    assoc = _association_types(genes_to_disease) if genes_to_disease else None
    gene_terms: dict[str, set[str]] = defaultdict(set)
    gene_modes: dict[str, set[str]] = defaultdict(set)
    with open(genes_to_phenotype, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene, raw_term = row["gene_symbol"], row["hpo_id"]
            kind = assoc.get((gene, row["disease_id"]), "UNKNOWN") if assoc else "MENDELIAN"
            if kind == "POLYGENIC":
                continue
            if raw_term in INHERITANCE_CODES:
                if kind == "MENDELIAN":
                    gene_modes[gene].add(INHERITANCE_CODES[raw_term])
                continue
            term = ontology.resolve(raw_term)
            if term and ontology.is_phenotype(term):
                gene_terms[gene].add(term)

    # Information content: IC(t) = -log(fraction of genes annotated with t or a descendant of t).
    # Rare, specific terms (e.g. "Exocrine pancreatic insufficiency") get high IC;
    # broad terms (e.g. "Abnormality of the nervous system") get low IC.
    counts: dict[str, int] = defaultdict(int)
    for terms in gene_terms.values():
        implied: set[str] = set()
        for t in terms:
            implied |= ontology.ancestors(t)
        for t in implied:
            counts[t] += 1
    n_genes = max(len(gene_terms), 1)
    ic = {t: -math.log(c / n_genes) for t, c in counts.items()}
    return GeneAnnotations(dict(gene_terms), dict(gene_modes), ic)
