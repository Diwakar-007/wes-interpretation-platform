import pytest

from wesinterp.acmg import Context
from wesinterp.models import Variant
from wesinterp.resources import ClinVarProteinRecord


def make_variant(**kwargs) -> Variant:
    """A rare missense variant by default; override any field."""
    base = dict(chrom="chr1", pos=1000, ref="C", alt="T", genotype="0/1", gene="GENE1",
                consequence=["missense_variant"], hgvsp="p.Arg100Trp", exon="3/10", popmax_af=None)
    base.update(kwargs)
    return Variant(**base)


@pytest.fixture
def ctx() -> Context:
    known = [
        # pathogenic p.Arg100Trp caused by a DIFFERENT nucleotide change at chr1:999
        ClinVarProteinRecord("GENE1", "Arg", 100, "Trp", "chr1", 999, "G", "A", "Pathogenic", 2, "111"),
        # pathogenic p.Gly200Asp
        ClinVarProteinRecord("GENE1", "Gly", 200, "Asp", "chr1", 2000, "G", "A", "Likely pathogenic", 1, "222"),
    ]
    index: dict = {}
    for r in known:
        index.setdefault(r.gene, {}).setdefault(r.aa_pos, []).append(r)
    return Context(loeuf={"GENE1": 0.2, "TOLERANT": 1.4}, clinvar_protein=index)
