from conftest import make_variant

from wesinterp.acmg import classify_variant
from wesinterp.acmg.classifier import clinvar_concordance, points_to_class
from wesinterp.clinvar_index import parse_missense, review_stars
from wesinterp.models import Evidence


def test_point_boundaries():
    assert points_to_class(10, []) == "Pathogenic"
    assert points_to_class(6, []) == "Likely pathogenic"
    assert points_to_class(5, []) == "Uncertain significance"
    assert points_to_class(0, []) == "Uncertain significance"
    assert points_to_class(-1, []) == "Likely benign"
    assert points_to_class(-6, []) == "Likely benign"
    assert points_to_class(-7, []) == "Benign"


def test_ba1_overrides_pathogenic_evidence():
    ev = [Evidence("PVS1", "VeryStrong", ""), Evidence("BA1", "StandAlone", "")]
    assert points_to_class(8, ev) == "Benign"


def test_rare_nonsense_in_intolerant_gene_is_likely_pathogenic(ctx):
    # PVS1 (8) + PM2_Supporting (1) = 9 -> Likely pathogenic
    r = classify_variant(make_variant(consequence=["stop_gained"], hgvsp="p.Arg100Ter"), ctx)
    assert r.criteria == "PVS1;PM2_Supporting"
    assert r.points == 9 and r.classification == "Likely pathogenic"


def test_ps1_plus_strong_revel_plus_absent_is_pathogenic(ctx):
    # PS1 (4) + PM2_Supporting (1) + PP3_Strong (4) = 9 -> Likely pathogenic
    r = classify_variant(make_variant(hgvsp="p.Arg100Trp", revel=0.95), ctx)
    assert r.points == 9 and r.classification == "Likely pathogenic"


def test_common_damaging_missense_is_benign(ctx):
    # like MTHFR p.Ala222Val: REVEL says damaging, but 30% allele frequency -> BA1 wins
    r = classify_variant(make_variant(hgvsp="p.Ala222Val", revel=0.84, popmax_af=0.32), ctx)
    assert "BA1" in r.criteria and "PP3_Moderate" in r.criteria
    assert r.classification == "Benign"


def test_clinvar_discordance_flag():
    v = make_variant(clinvar_sig="Pathogenic", clinvar_review="criteria_provided&_multiple_submitters&_no_conflicts")
    assert clinvar_concordance(v, "Likely pathogenic") == "concordant"
    assert clinvar_concordance(v, "Likely benign") == "discordant"
    one_star = make_variant(clinvar_sig="Pathogenic", clinvar_review="criteria_provided&_single_submitter")
    assert clinvar_concordance(one_star, "Benign") == ""


def test_review_stars():
    assert review_stars("practice guideline") == 4
    assert review_stars("reviewed_by_expert_panel") == 3
    assert review_stars("criteria provided, multiple submitters, no conflicts") == 2
    assert review_stars("criteria provided, single submitter") == 1
    assert review_stars("no assertion criteria provided") == 0


def test_parse_missense():
    assert parse_missense("p.Leu858Arg") == ("Leu", 858, "Arg")
    assert parse_missense("p.Arg100Ter") is None     # nonsense, not missense
    assert parse_missense("p.Phe508del") is None
    assert parse_missense("p.Val276=") is None
