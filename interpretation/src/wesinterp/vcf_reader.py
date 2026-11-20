"""Read a VEP-annotated VCF into Variant objects (one per variant allele per PICKed gene consequence)."""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path
from urllib.parse import unquote

from .config import POPMAX_POPULATIONS
from .models import Variant


def _max_float(value: str | None) -> float | None:
    """VEP joins multiple values with '&' (e.g. several co-located known variants). Take the max."""
    if not value:
        return None
    numbers = []
    for part in value.split("&"):
        try:
            numbers.append(float(part))
        except ValueError:
            continue
    return max(numbers) if numbers else None


def _strip_hgvs(value: str) -> str:
    """'ENST00000318948.7:c.2987A>C' -> 'c.2987A>C'; also decode %3D -> '='."""
    return unquote(value.split(":", 1)[1]) if ":" in value else unquote(value)


def csq_fields_from_header(description: str) -> list[str]:
    m = re.search(r"Format: ([^\"]+)", description)
    if not m:
        raise ValueError("CSQ header has no 'Format:' field list; is this a VEP VCF?")
    return m.group(1).strip().split("|")


def variant_from_csq(chrom: str, pos: int, ref: str, alt: str, genotype: str, csq: dict[str, str]) -> Variant:
    popmax_values = [
        _max_float(csq.get(f"gnomAD{src}_{pop}_AF"))
        for src in ("e", "g") for pop in POPMAX_POPULATIONS
    ]
    popmax_values = [v for v in popmax_values if v is not None]
    overall = [v for v in (_max_float(csq.get("gnomADe_AF")), _max_float(csq.get("gnomADg_AF"))) if v is not None]
    return Variant(
        chrom=chrom, pos=pos, ref=ref, alt=alt, genotype=genotype,
        gene=csq.get("SYMBOL", ""),
        transcript=csq.get("MANE_SELECT") or csq.get("Feature", ""),
        consequence=[c for c in csq.get("Consequence", "").split("&") if c],
        hgvsc=_strip_hgvs(csq.get("HGVSc", "")),
        hgvsp=_strip_hgvs(csq.get("HGVSp", "")),
        exon=csq.get("EXON", ""),
        popmax_af=max(popmax_values) if popmax_values else (max(overall) if overall else None),
        gnomad_af=max(overall) if overall else None,
        revel=_max_float(csq.get("REVEL")),
        am_class=csq.get("am_class", ""),
        am_score=_max_float(csq.get("am_pathogenicity")),
        clinvar_sig=csq.get("ClinVar_CLNSIG", ""),
        clinvar_review=csq.get("ClinVar_CLNREVSTAT", ""),
        clinvar_disease=csq.get("ClinVar_CLNDN", "").replace("_", " "),
        rsid="&".join(x for x in csq.get("Existing_variation", "").split("&") if x.startswith("rs")),
    )


def read_vep_vcf(path: str | Path, pass_only: bool = True, picked_only: bool = True) -> Iterator[Variant]:
    from cyvcf2 import VCF   # imported here so the rest of the package works without cyvcf2

    vcf = VCF(str(path))
    fields = csq_fields_from_header(vcf.get_header_type("CSQ")["Description"])
    for rec in vcf:
        if pass_only and rec.FILTER is not None:      # cyvcf2: None means PASS
            continue
        gt = rec.genotypes[0] if rec.genotypes else [-1, -1, False]
        alleles = ["." if a < 0 else str(a) for a in gt[:-1]]
        genotype = ("|" if gt[-1] else "/").join(alleles)
        csq_raw = rec.INFO.get("CSQ")
        if not csq_raw:
            continue
        for entry in csq_raw.split(","):
            csq = dict(zip(fields, entry.split("|")))
            if picked_only and csq.get("PICK") != "1":
                continue
            yield variant_from_csq(rec.CHROM, rec.POS, rec.REF, rec.ALT[0], genotype, csq)
