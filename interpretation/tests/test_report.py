import json

from wesinterp.report import build_report, select_findings


def _variant(pos, gt, cls, gene="G"):
    return {"chrom": "chr1", "pos": pos, "ref": "A", "alt": "T", "genotype": gt, "gene": gene,
            "classification": cls, "transcript": "NM_1.1", "hgvsc": "c.1A>T", "hgvsp": "p.Met1Leu",
            "consequence": ["missense_variant"], "popmax_af": None, "revel": 0.9, "am_class": "",
            "clinvar_sig": "", "clinvar_review": "", "clinvar_concordance": "", "rsid": "",
            "points": 6, "evidence": [{"code": "PM2_Supporting", "points": 1, "reason": "absent"}]}


def _data():
    return {
        "patient_terms": [{"id": "HP:0001250", "name": "Seizure"}],
        "variants": [_variant(1, "0/1", "Likely pathogenic", "DOM"),
                     _variant(2, "0/1", "Pathogenic", "REC"),
                     _variant(3, "1/1", "Uncertain significance", "MATCH"),
                     _variant(4, "0/1", "Uncertain significance", "WEAK")],
        "candidates": [
            {"gene": "DOM", "variants": ["chr1-1-A-T"], "phenotype_score": 0.9, "fit_weight": 1.0,
             "inheritance": ["AD"], "inheritance_fit": "fits AD", "score": 0.8, "matched_phenotypes": []},
            {"gene": "REC", "variants": ["chr1-2-A-T"], "phenotype_score": 0.2, "fit_weight": 0.3,
             "inheritance": ["AR"], "inheritance_fit": "AR gene but single heterozygous variant (carrier?)",
             "score": 0.1, "matched_phenotypes": []},
            {"gene": "MATCH", "variants": ["chr1-3-A-T"], "phenotype_score": 0.8, "fit_weight": 1.0,
             "inheritance": ["AR"], "inheritance_fit": "fits AR (homozygous)", "score": 0.4, "matched_phenotypes": []},
            {"gene": "WEAK", "variants": ["chr1-4-A-T"], "phenotype_score": 0.1, "fit_weight": 1.0,
             "inheritance": ["AD"], "inheritance_fit": "fits AD", "score": 0.05, "matched_phenotypes": []},
        ],
    }


def test_select_findings_tiers():
    f = select_findings(_data())
    assert [(x["candidate"]["gene"], x["tier"]) for x in f["primary"]] == [("DOM", "Reportable"), ("MATCH", "VUS of interest")]
    assert [x["candidate"]["gene"] for x in f["carrier"]] == ["REC"]
    assert [x["candidate"]["gene"] for x in f["other"]] == ["WEAK"]


def test_same_position_different_genes_not_mixed_up():
    data = _data()
    # an antisense-gene record at the same position as DOM's variant
    data["variants"].append({**_variant(1, "0/1", "Benign", "DOM-AS1"), "transcript": "WRONG"})
    f = select_findings(data)
    assert f["primary"][0]["variant"]["transcript"] == "NM_1.1"


def test_build_report_html(tmp_path):
    path = tmp_path / "c.json"
    path.write_text(json.dumps(_data()))
    ctx = build_report(path, tmp_path / "r.html", sample_id="S1", indication="Epilepsy")
    html = (tmp_path / "r.html").read_text()
    assert "POSITIVE" in html and "DOM" in html and "Carrier findings" in html and "Seizure" in html
    assert ctx["sample_id"] == "S1"
