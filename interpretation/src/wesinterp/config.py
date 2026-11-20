"""All classification thresholds in one place, with their sources.

These are generic (gene-agnostic) defaults. Expert panels (ClinGen VCEPs) publish
gene-specific values, e.g. different allele-frequency cut-offs for CFTR or MYH7.
"""

from dataclasses import dataclass

# gnomAD continental populations used for the "popmax" allele frequency.
# Founder / small populations (ASJ, FIN, MID, AMI, "remaining") are excluded, as ClinGen recommends.
POPMAX_POPULATIONS = ("AFR", "AMR", "EAS", "NFE", "SAS")


@dataclass(frozen=True)
class Thresholds:
    # --- population frequency ---
    ba1_af: float = 0.05        # BA1: > 5% in any continental population (ClinGen SVI, Ghosh 2018)
    # BS1: allele frequency higher than expected for the disorder. Depends on inheritance:
    #  - recessive: healthy carriers are common (CFTR p.Phe508del is ~1.5% in Europeans), so the
    #    cut-off must be high. 0.02 ~ the allele frequency behind a 1/2500 recessive disease.
    #  - dominant: one allele causes disease, so even 0.1% is too common for a rare disorder.
    #  - unknown / mixed inheritance: generic 1%.
    bs1_af: float = 0.01
    bs1_af_recessive: float = 0.02
    bs1_af_dominant: float = 0.001
    pm2_af: float = 0.0001      # PM2_Supporting: absent, or extremely rare (ClinGen SVI 2020 downgrade)

    # --- PP3 / BP4 from REVEL (Pejaver et al. 2022, AJHG 109:2163) ---
    revel_pp3_strong: float = 0.932
    revel_pp3_moderate: float = 0.773
    revel_pp3_supporting: float = 0.644
    revel_bp4_supporting: float = 0.290
    revel_bp4_moderate: float = 0.183
    revel_bp4_strong: float = 0.016

    # --- PVS1: gene must be loss-of-function intolerant ---
    # LOEUF < 0.6 is gnomAD's recommended cut-off for LoF-intolerant genes; used here as a
    # proxy for "loss of function is a known disease mechanism".
    loeuf_max: float = 0.6


DEFAULT_THRESHOLDS = Thresholds()

# Tavtigian 2020 point ranges for the final 5-tier classification
CLASS_POINTS = (
    (10, "Pathogenic"),
    (6, "Likely pathogenic"),
    (0, "Uncertain significance"),
    (-6, "Likely benign"),
)
BENIGN = "Benign"
