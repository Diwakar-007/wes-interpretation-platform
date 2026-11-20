"""Minimal parser and graph utilities for the Human Phenotype Ontology (hp.obo)."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

PHENOTYPIC_ABNORMALITY = "HP:0000118"
MODE_OF_INHERITANCE = "HP:0000005"


@dataclass
class Term:
    id: str
    name: str = ""
    synonyms: list[str] = field(default_factory=list)
    parents: list[str] = field(default_factory=list)
    alt_ids: list[str] = field(default_factory=list)
    obsolete: bool = False
    replaced_by: str = ""


class Ontology:
    def __init__(self, terms: dict[str, Term]):
        self.terms = terms
        # map alternative / obsolete IDs to their current term
        self._alias: dict[str, str] = {}
        for t in terms.values():
            for alt in t.alt_ids:
                self._alias[alt] = t.id
            if t.obsolete and t.replaced_by:
                self._alias[t.id] = t.replaced_by

    def resolve(self, term_id: str) -> str | None:
        """Current ID for a term (follows alt_id / replaced_by); None if unknown or obsolete."""
        term_id = self._alias.get(term_id, term_id)
        t = self.terms.get(term_id)
        return None if t is None or t.obsolete else term_id

    def name(self, term_id: str) -> str:
        t = self.terms.get(term_id)
        return t.name if t else term_id

    @lru_cache(maxsize=None)
    def ancestors(self, term_id: str) -> frozenset[str]:
        """The term itself plus all terms above it (following is_a)."""
        result = {term_id}
        for parent in self.terms[term_id].parents if term_id in self.terms else ():
            result |= self.ancestors(parent)
        return frozenset(result)

    def is_phenotype(self, term_id: str) -> bool:
        return PHENOTYPIC_ABNORMALITY in self.ancestors(term_id)

    def search(self, query: str, limit: int = 15) -> list[tuple[str, str]]:
        """Simple text search over names and synonyms (for autocomplete). Phenotype terms only."""
        q = query.lower().strip()
        if not q:
            return []
        if q.startswith("hp:"):
            tid = self.resolve(query.upper())
            return [(tid, self.name(tid))] if tid else []
        scored = []
        for t in self.terms.values():
            if t.obsolete or not self.is_phenotype(t.id):
                continue
            name = t.name.lower()
            if name == q:
                rank = 0
            elif name.startswith(q):
                rank = 1
            elif q in name:
                rank = 2
            elif any(q in s.lower() for s in t.synonyms):
                rank = 3
            else:
                continue
            scored.append((rank, len(t.name), t.id, t.name))
        scored.sort()
        return [(tid, name) for _, _, tid, name in scored[:limit]]


def load_ontology(path: str | Path) -> Ontology:
    terms: dict[str, Term] = {}
    current: Term | None = None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "[Term]":
                current = Term(id="")
                continue
            if line.startswith("[") and line.endswith("]"):   # [Typedef] etc.
                current = None
                continue
            if current is None or ": " not in line:
                if not line and current is not None and current.id:
                    terms[current.id] = current
                    current = None
                continue
            key, value = line.split(": ", 1)
            if key == "id":
                current.id = value
            elif key == "name":
                current.name = value
            elif key == "is_a":
                current.parents.append(value.split(" ! ")[0].strip())
            elif key == "synonym":
                current.synonyms.append(value.split('"')[1] if '"' in value else value)
            elif key == "alt_id":
                current.alt_ids.append(value)
            elif key == "is_obsolete":
                current.obsolete = value == "true"
            elif key == "replaced_by":
                current.replaced_by = value
    if current is not None and current.id:
        terms[current.id] = current
    return Ontology(terms)
