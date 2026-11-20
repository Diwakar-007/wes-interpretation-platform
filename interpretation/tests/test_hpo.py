"""HPO parsing, information content, phenotype similarity, inheritance-aware BS1 and ranking."""

import pytest
from conftest import make_variant

from wesinterp.acmg import Context, classify_variant
from wesinterp.acmg import criteria as c
from wesinterp.hpo import load_gene_annotations, load_ontology, phenotype_score
from wesinterp.prioritize import inheritance_fit, prioritize, zygosity

# Tiny ontology:
#   HP:0000118 Phenotypic abnormality
#     ├── HP:1 Seizure
#     │     └── HP:2 Focal seizure
#     └── HP:3 Lung disease
#           └── HP:4 Bronchiectasis
OBO = """format-version: 1.2

[Term]
id: HP:0000118
name: Phenotypic abnormality

[Term]
id: HP:1
name: Seizure
synonym: "Epileptic seizure" EXACT []
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:2
name: Focal seizure
is_a: HP:1 ! Seizure

[Term]
id: HP:3
name: Lung disease
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:4
name: Bronchiectasis
alt_id: HP:44
is_a: HP:3 ! Lung disease

[Term]
id: HP:5
name: Old term
is_obsolete: true
replaced_by: HP:4
"""

G2P = """ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id
1\tEPI1\tHP:2\tFocal seizure\t-\tOMIM:1
1\tEPI1\tHP:0000006\tAutosomal dominant inheritance\t-\tOMIM:1
2\tLUNG1\tHP:4\tBronchiectasis\t-\tOMIM:2
2\tLUNG1\tHP:0000007\tAutosomal recessive inheritance\t-\tOMIM:2
3\tOTHER\tHP:3\tLung disease\t-\tOMIM:3
"""


@pytest.fixture
def hpo(tmp_path):
    (tmp_path / "hp.obo").write_text(OBO)
    (tmp_path / "g2p.txt").write_text(G2P)
    onto = load_ontology(tmp_path / "hp.obo")
    return onto, load_gene_annotations(tmp_path / "g2p.txt", onto)


def test_ontology_parsing_and_aliases(hpo):
    onto, _ = hpo
    assert onto.ancestors("HP:2") == {"HP:2", "HP:1", "HP:0000118"}
    assert onto.resolve("HP:44") == "HP:4"          # alt_id
    assert onto.resolve("HP:5") == "HP:4"           # obsolete -> replaced_by
    assert onto.search("seiz")[0][1] == "Seizure"
    assert onto.search("epileptic")[0][0] == "HP:1"  # synonym match


def test_inheritance_modes_and_information_content(hpo):
    _, ann = hpo
    assert ann.gene_modes == {"EPI1": {"AD"}, "LUNG1": {"AR"}}
    # Bronchiectasis (1 of 3 genes) is more specific than Lung disease (2 of 3 genes)
    assert ann.information_content["HP:4"] > ann.information_content["HP:3"] > 0
    assert ann.information_content["HP:0000118"] == 0      # every gene has it


def test_polygenic_associations_ignored(tmp_path):
    (tmp_path / "hp.obo").write_text(OBO)
    g2p = G2P + "2\tLUNG1\tHP:0000006\tAutosomal dominant inheritance\t-\tOMIM:9\n"
    (tmp_path / "g2p.txt").write_text(g2p)
    (tmp_path / "g2d.txt").write_text(
        "ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n"
        "2\tLUNG1\tMENDELIAN\tOMIM:2\tx\n"
        "2\tLUNG1\tPOLYGENIC\tOMIM:9\tx\n")
    onto = load_ontology(tmp_path / "hp.obo")
    ann = load_gene_annotations(tmp_path / "g2p.txt", onto, tmp_path / "g2d.txt")
    assert ann.gene_modes["LUNG1"] == {"AR"}     # the AD label came from a risk-factor association


def test_phenotype_score(hpo):
    onto, ann = hpo
    exact, _ = phenotype_score(["HP:4"], "LUNG1", ann, onto)
    partial, _ = phenotype_score(["HP:4"], "OTHER", ann, onto)   # only shares "Lung disease"
    none, _ = phenotype_score(["HP:4"], "EPI1", ann, onto)
    assert exact == pytest.approx(1.0)
    assert 0 < partial < exact
    assert none == 0


def test_bs1_is_inheritance_aware():
    ctx = Context(gene_modes={"RECGENE": {"AR"}, "DOMGENE": {"AD"}})
    # 1.5% popmax: like CFTR p.Phe508del
    assert c.bs1(make_variant(gene="RECGENE", popmax_af=0.015), ctx) is None        # carriers are common
    assert c.bs1(make_variant(gene="DOMGENE", popmax_af=0.015), ctx).points == -4   # far too common for AD
    assert c.bs1(make_variant(gene="UNKNOWN", popmax_af=0.015), ctx) is not None    # generic 1%


def test_zygosity():
    assert zygosity("0/1") == "heterozygous"
    assert zygosity("1|1") == "homozygous"
    assert zygosity("1") == "hemizygous"


def test_inheritance_fit(ctx):
    het = classify_variant(make_variant(genotype="0/1"), ctx)
    het2 = classify_variant(make_variant(pos=2000, genotype="0/1"), ctx)
    hom = classify_variant(make_variant(genotype="1/1"), ctx)
    assert inheritance_fit({"AR"}, [hom])[1] == 1.0
    assert inheritance_fit({"AR"}, [het, het2])[1] == 0.8     # possible compound het
    assert inheritance_fit({"AR"}, [het])[1] == 0.3           # likely just a carrier
    assert inheritance_fit({"AD"}, [het])[1] == 1.0


def test_prioritize_ranks_matching_gene_first(hpo):
    onto, ann = hpo
    ctx = Context(gene_modes=ann.gene_modes)
    results = [
        classify_variant(make_variant(gene="EPI1", genotype="0/1"), ctx),
        classify_variant(make_variant(gene="LUNG1", pos=5, genotype="1/1"), ctx),
    ]
    ranked = prioritize(results, ["HP:4"], ann, onto)   # patient has bronchiectasis
    assert ranked[0].gene == "LUNG1"
    assert ranked[0].inheritance_fit == "fits AR (homozygous)"
