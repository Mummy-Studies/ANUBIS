"""
ANUBIS — Ancient DNA damage profiler for MetaPhlAn SGBs.

Integrates MetaPhlAn relative abundance profiles with PyDamage damage
analysis, producing per-SGB damage metrics, BH-corrected q-values,
damage plots, and optionally rescaled/masked BAM files for downstream
StrainPhlAn genotyping.
"""

__version__ = "0.1.0"
__author__ = "Mummy Studies"

from .core import main

__all__ = ["main", "__version__"]
