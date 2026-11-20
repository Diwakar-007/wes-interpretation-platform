from conftest import make_variant

from wesinterp.acmg import criteria as c


# ---------------- PVS1
def test_pvs1_stop_gained_in_intolerant_gene_is_very_strong(ctx):
    e = c.pvs1(make_variant(consequence=["stop_gained"]), ctx)
    assert e.code == "PVS1" and e.strength == "VeryStrong" and e.points == 8


def test_pvs1_not_applied_in_lof_tolerant_gene(ctx):
    assert c.pvs1(make_variant(gene="TOLERANT", consequence=["stop_gained"]), ctx) is None


def test_pvs1_not_applied_when_gene_constraint_unknown(ctx):
    assert c.pvs1(make_variant(gene="UNKNOWN", consequence=["frameshift_variant"]), ctx) is None


def test_pvs1_last_exon_truncation_downgraded_to_moderate(ctx):
    e = c.pvs1(make_variant(consequence=["frameshift_variant"], exon="10/10"), ctx)
    assert e.strength == "Moderate" and e.label == "PVS1_Moderate"


def test_pvs1_splice_donor_in_last_exon_stays_very_strong(ctx):
    e = c.pvs1(make_variant(consequence=["splice_donor_variant"], exon="10/10"), ctx)
    assert e.strength == "VeryStrong"


def test_pvs1_start_lost_is_moderate(ctx):
    assert c.pvs1(make_variant(consequence=["start_lost"]), ctx).strength == "Moderate"


def test_pvs1_missense_not_applied(ctx):
    assert c.pvs1(make_variant(), ctx) is None


# ---------------- PS1 / PM5
def test_ps1_same_amino_acid_change_different_nucleotide(ctx):
    e = c.ps1_pm5(make_variant(hgvsp="p.Arg100Trp"), ctx)
    assert e.code == "PS1" and e.points == 4


def test_ps1_not_applied_to_the_same_genomic_variant(ctx):
    # the query IS the ClinVar variant -> not PS1 (that would be circular)
    assert c.ps1_pm5(make_variant(pos=999, ref="G", alt="A", hgvsp="p.Arg100Trp"), ctx) is None


def test_pm5_different_change_same_residue(ctx):
    e = c.ps1_pm5(make_variant(hgvsp="p.Gly200Val"), ctx)
    assert e.code == "PM5" and e.points == 2


def test_ps1_pm5_need_missense(ctx):
    assert c.ps1_pm5(make_variant(consequence=["synonymous_variant"], hgvsp="p.Arg100="), ctx) is None


# ---------------- frequency: PM2 / BS1 / BA1
def test_pm2_supporting_when_absent(ctx):
    e = c.pm2(make_variant(popmax_af=None), ctx)
    assert e.label == "PM2_Supporting" and e.points == 1


def test_pm2_not_applied_when_common(ctx):
    assert c.pm2(make_variant(popmax_af=0.001), ctx) is None


def test_ba1_above_5_percent(ctx):
    assert c.ba1(make_variant(popmax_af=0.30), ctx).code == "BA1"
    assert c.ba1(make_variant(popmax_af=0.04), ctx) is None


def test_bs1_between_1_and_5_percent(ctx):
    assert c.bs1(make_variant(popmax_af=0.02), ctx).points == -4
    assert c.bs1(make_variant(popmax_af=0.30), ctx) is None   # BA1 territory, not double counted
    assert c.bs1(make_variant(popmax_af=0.005), ctx) is None


# ---------------- PP3 / BP4 (REVEL, Pejaver 2022)
def test_pp3_strengths(ctx):
    assert c.pp3_bp4(make_variant(revel=0.95), ctx).label == "PP3_Strong"
    assert c.pp3_bp4(make_variant(revel=0.80), ctx).label == "PP3_Moderate"
    assert c.pp3_bp4(make_variant(revel=0.70), ctx).label == "PP3"


def test_bp4_strengths(ctx):
    assert c.pp3_bp4(make_variant(revel=0.25), ctx).label == "BP4"
    assert c.pp3_bp4(make_variant(revel=0.10), ctx).label == "BP4_Moderate"
    assert c.pp3_bp4(make_variant(revel=0.01), ctx).label == "BP4_Strong"


def test_revel_indeterminate_range_gives_nothing(ctx):
    assert c.pp3_bp4(make_variant(revel=0.5), ctx) is None


# ---------------- PM4 / BP7
def test_pm4_inframe_deletion(ctx):
    assert c.pm4(make_variant(consequence=["inframe_deletion"]), ctx).points == 2


def test_bp7_synonymous_outside_splice_region(ctx):
    assert c.bp7(make_variant(consequence=["synonymous_variant"]), ctx).points == -1
    assert c.bp7(make_variant(consequence=["synonymous_variant", "splice_region_variant"]), ctx) is None
