# WES Interpretation Platform

An end-to-end whole-exome sequencing (WES) platform. A user uploads paired FASTQ files, the pipeline
produces a VCF, and the user adds HPO phenotypes. The platform then ranks candidate variants,
applies ACMG/AMP classification and generates a clinical-style PDF report.

> Status: work in progress (portfolio / learning project, not for clinical use).

## Architecture

```
FASTQ ─► fastp ─► bwa mem ─► MarkDuplicates ─► BQSR ─► HaplotypeCaller ─► filter/normalise ─► VEP
                                                                                              │
          HPO phenotypes ─► gene ranking ─┐                                                    ▼
                                          ├──► ACMG engine (Bayesian points) ─► review UI ─► PDF report
                     annotated variants ──┘
```

| Folder | Contents |
|---|---|
| `pipeline/` | Nextflow DSL2 pipeline: FASTQ → annotated VCF |
| `interpretation/` | Python package: VCF parsing, HPO ranking, ACMG criteria |
| `reporting/` | Report templates and PDF generation |
| `web/` | FastAPI web app (upload, job status, review, report) |
| `validation/` | Benchmarking against GIAB truth sets with hap.py |
| `envs/`, `containers/` | Reproducible software environments |
| `scripts/setup/` | Installation and reference download scripts |

Data lives outside the repo in `~/ngs-data/` (`refs/`, `test/`, `runs/`).

## Setup (Ubuntu / WSL2)

```bash
bash scripts/setup/01_install_tools.sh       # Miniforge, conda env 'ngs', Nextflow
conda activate ngs
bash scripts/setup/02_download_test_data.sh  # tiny chr22 test set (MB)
bash scripts/setup/03_download_grch38.sh     # GRCh38 + BWA index + known sites (~11 GB)
```

## Hardware target

Designed to run on a 16 GB RAM / 4-core laptop (WSL2 limited to 12 GB): classic `bwa mem`,
HaplotypeCaller restricted to exome target intervals, and the offline VEP cache.
