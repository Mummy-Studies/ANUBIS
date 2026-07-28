# Changelog

All notable changes to ANUBIS are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versions follow [Semantic Versioning](https://semver.org/).

## [0.1.0] — 2025-01-01

### Added
- Per-SGB damage assessment (C→T at 5′ end, G→A at 3′ end) via PyDamage group mode
- Geometric decay damage model with confidence intervals
- Benjamini–Hochberg FDR correction across all analysed SGBs
- Predicted accuracy score (logistic regression GLM from PyDamage)
- Classical aDNA damage profile plots (`--plot`)
- Model-based base-quality rescaling for downstream StrainPhlAn genotyping (`--rescale`)
- Flexible terminal base masking at 5′ and/or 3′ ends (`--mask-5prime`, `--mask-3prime`)
- Clade-specific mode: target one or more SGBs directly (`-c`/`--clade`)
- SAM header index cache for fast repeated runs (`.anubis_idx.json.gz`)
- Multi-threaded BAM sorting and PyDamage analysis (`--threads`)
- Pip-installable Python package (`pip install anubis-adna`)
