"""Command line interface.

    wesinterp build-clinvar-index --variant-summary variant_summary.txt.gz --out clinvar_protein_index.tsv.gz
    wesinterp classify --vcf sample.vep.vcf.gz --loeuf constraint.tsv --clinvar-index idx.tsv.gz --out sample.acmg.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter

from . import __version__
from .acmg import Context, classify_variant
from .clinvar_index import build_index
from .resources import load_clinvar_protein_index, load_loeuf
from .vcf_reader import read_vep_vcf

TSV_COLUMNS = ["chrom", "pos", "ref", "alt", "genotype", "gene", "transcript", "consequence", "hgvsc", "hgvsp",
               "rsid", "popmax_af", "revel", "am_class", "clinvar_sig", "clinvar_review", "clinvar_disease",
               "criteria", "points", "classification", "clinvar_concordance", "evidence_detail"]


def _fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


def cmd_classify(args: argparse.Namespace) -> int:
    ctx = Context(
        loeuf=load_loeuf(args.loeuf) if args.loeuf else {},
        clinvar_protein=load_clinvar_protein_index(args.clinvar_index) if args.clinvar_index else {},
    )
    results = [classify_variant(v, ctx) for v in read_vep_vcf(args.vcf, pass_only=not args.all_filters)]

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(TSV_COLUMNS)
        for r in results:
            v = r.variant
            w.writerow([_fmt(x) for x in (
                v.chrom, v.pos, v.ref, v.alt, v.genotype, v.gene, v.transcript, "&".join(v.consequence),
                v.hgvsc, v.hgvsp, v.rsid, v.popmax_af, v.revel, v.am_class, v.clinvar_sig, v.clinvar_review,
                v.clinvar_disease, r.criteria, r.points, r.classification, r.clinvar_concordance,
                " | ".join(f"{e.label}: {e.reason}" for e in r.evidence))])

    if args.json:
        with open(args.json, "w") as fh:
            json.dump([{
                **{k: getattr(r.variant, k) for k in r.variant.__dataclass_fields__},
                "evidence": [{"code": e.label, "points": e.points, "reason": e.reason} for e in r.evidence],
                "points": r.points, "classification": r.classification,
                "clinvar_concordance": r.clinvar_concordance,
            } for r in results], fh, indent=1)

    counts = Counter(r.classification for r in results)
    print(f"Classified {len(results)} variant-gene records -> {args.out}", file=sys.stderr)
    for label in ("Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"):
        print(f"  {label:24s} {counts.get(label, 0)}", file=sys.stderr)
    discordant = sum(r.clinvar_concordance == "discordant" for r in results)
    if discordant:
        print(f"  ! {discordant} discordant with ClinVar (>=2 stars): review manually", file=sys.stderr)
    return 0


def cmd_build_index(args: argparse.Namespace) -> int:
    n = build_index(args.variant_summary, args.out, min_stars=args.min_stars)
    print(f"Wrote {n} pathogenic missense records -> {args.out}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="wesinterp", description="ACMG variant classification")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("classify", help="classify variants in a VEP-annotated VCF")
    c.add_argument("--vcf", required=True)
    c.add_argument("--out", required=True, help="output TSV")
    c.add_argument("--json", help="optional JSON output (used by the web app)")
    c.add_argument("--loeuf", help="gnomAD constraint TSV (enables PVS1)")
    c.add_argument("--clinvar-index", help="ClinVar protein index (enables PS1/PM5)")
    c.add_argument("--all-filters", action="store_true", help="include non-PASS variants")
    c.set_defaults(func=cmd_classify)

    b = sub.add_parser("build-clinvar-index", help="build the PS1/PM5 index from ClinVar variant_summary")
    b.add_argument("--variant-summary", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--min-stars", type=int, default=1)
    b.set_defaults(func=cmd_build_index)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
