"""Build an HTML + PDF report from the classification JSON (output of `wesinterp classify --json`).

Which variants go where:
  Primary findings   P/LP variants in a candidate gene whose inheritance fits (not a lone het in an AR gene),
                     plus VUS in genes with a strong phenotype match (>= 0.6) and a fitting inheritance
                     ("VUS of interest": these need expert review before any clinical use).
  Carrier findings   P/LP heterozygous variants in autosomal recessive genes.
  Other candidates   the next highest-ranked genes, for the analyst.
"""

from __future__ import annotations

import json
import uuid
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, PackageLoader, select_autoescape

from .. import __version__
from ..config import DEFAULT_THRESHOLDS
from .qc import read_qc

PLP = ("Pathogenic", "Likely pathogenic")
CLASS_CSS = {"Pathogenic": "p", "Likely pathogenic": "lp", "Uncertain significance": "vus",
             "Likely benign": "lb", "Benign": "b"}


def zygosity(gt: str) -> str:
    alleles = gt.replace("|", "/").split("/")
    alt = [a for a in alleles if a not in ("0", ".")]
    if len(alleles) == 1:
        return "Hemizygous" if alt else "Reference"
    return "Homozygous" if len(alt) == 2 else "Heterozygous"


def select_findings(data: dict, max_other: int = 10, min_vus_phenotype: float = 0.6) -> dict[str, list[dict]]:
    # one variant position can have records for several genes (e.g. CFTR and CFTR-AS1),
    # so look records up by position AND gene
    variants = {(f"{v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}", v["gene"]): v for v in data["variants"]}
    primary, carrier, other = [], [], []
    for rank, cand in enumerate(data.get("candidates", []), 1):
        for key in cand["variants"]:
            v = variants[(key, cand["gene"])]
            row = {"rank": rank, "candidate": cand, "variant": v, "zygosity": zygosity(v["genotype"])}
            is_carrier = "carrier" in cand["inheritance_fit"]
            if v["classification"] in PLP and not is_carrier:
                primary.append({**row, "tier": "Reportable"})
            elif v["classification"] in PLP and is_carrier:
                carrier.append(row)
            elif cand["phenotype_score"] >= min_vus_phenotype and cand.get("fit_weight", 0) >= 0.8:
                primary.append({**row, "tier": "VUS of interest"})
            elif len(other) < max_other:
                other.append(row)
    return {"primary": primary, "carrier": carrier, "other": other}


def build_report(classification_json: str | Path, out_html: str | Path, out_pdf: str | Path | None = None,
                 sample_id: str = "SAMPLE", sex: str = "Not provided", indication: str = "",
                 sample_dir: str | Path | None = None) -> dict:
    data = json.loads(Path(classification_json).read_text())
    findings = select_findings(data)
    env = Environment(loader=PackageLoader("wesinterp.report", "templates"),
                      autoescape=select_autoescape(["html", "j2"]))
    env.filters["af"] = lambda x: "absent" if x is None else (f"{x:.2e}" if x < 0.001 else f"{x:.4f}")
    env.filters["num"] = lambda x, d=1: f"{x:,.{d}f}" if isinstance(x, (int, float)) else "-"

    context = {
        "report_id": uuid.uuid4().hex[:10].upper(),
        "generated": datetime.now().strftime("%d %b %Y, %H:%M"),
        "sample_id": sample_id, "sex": sex, "indication": indication,
        "patient_terms": data.get("patient_terms", []),
        "resources": {**data.get("resources", {}), "wesinterp": __version__},
        "findings": findings,
        "n_variants": len(data["variants"]),
        "n_candidates": len(data.get("candidates", [])),
        "qc": read_qc(sample_dir, sample_id),
        "thresholds": DEFAULT_THRESHOLDS,
        "class_css": CLASS_CSS,
    }
    html = env.get_template("report.html.j2").render(**context)
    Path(out_html).write_text(html)
    if out_pdf:
        from weasyprint import HTML
        HTML(string=html, base_url=str(Path(out_html).parent)).write_pdf(out_pdf)
    return context
