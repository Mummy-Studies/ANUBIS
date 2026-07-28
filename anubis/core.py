"""
ANUBIS core — ancient DNA damage profiler for MetaPhlAn SGBs.

Entry point: anubis.core:main  (registered as the 'anubis' console script)
"""

import argparse
import sys
import os
import re
import tempfile
import shutil
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import gzip

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# During development the pydamage clone lives one directory above this file
# (i.e. next to pyproject.toml).  When installed via pip that directory does
# not exist, so the path-insert is silently skipped and the installed pydamage
# package is used instead.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PYDAMAGE_DIR = os.path.join(_REPO_ROOT, "pydamage")
if os.path.isdir(_PYDAMAGE_DIR):
    sys.path.insert(0, _PYDAMAGE_DIR)

from pydamage.damage import test_damage
from pydamage.accuracy_model import glm_predict
from pydamage.models import glm_model_params

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ── MetaPhlAn parsing ────────────────────────────────────────────────────────

def parse_metaphlan(tsv_path: str) -> pd.DataFrame:
    """
    Parse MetaPhlAn output TSV and return a DataFrame with one row per
    SGB (t__ level), containing:
        sgb_id          e.g. 'SGB6653'
        species         human-readable species name
        relative_abundance  (%)
        lineage         full lineage string
    """
    rows = []
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            lineage, ncbi_ids, abundance = parts[0], parts[1], parts[2]
            additional = parts[3] if len(parts) > 3 else ""

            # Only SGB-level rows; strip _group suffix — SAM marker-gene
            # references use the bare numeric SGB ID only.
            m = re.search(r"\|t__SGB(\d+)", lineage)
            if not m:
                continue
            sgb_id = "SGB" + m.group(1)

            # Species name from the s__ level
            s_match = re.search(r"\|s__([^|]+)\|", lineage)
            species = s_match.group(1).replace("_", " ") if s_match else sgb_id

            rows.append(
                {
                    "sgb_id": sgb_id,
                    "species": species,
                    "relative_abundance": float(abundance),
                    "lineage": lineage,
                }
            )

    df = pd.DataFrame(rows).sort_values("relative_abundance", ascending=False)
    return df.reset_index(drop=True)


def apply_filters(
    df: pd.DataFrame,
    min_abundance: float | None,
    top_n: int | None,
) -> pd.DataFrame:
    if min_abundance is not None:
        df = df[df["relative_abundance"] >= min_abundance]
    if top_n is not None:
        df = df.head(top_n)
    return df.reset_index(drop=True)


def _normalize_clade(token: str) -> str | None:
    """
    Accept clade identifiers in multiple formats and return 'SGB<digits>'.

    Supported formats:
        t__SGB12546   (full MetaPhlAn clade suffix)
        SGB12546      (bare SGB ID)
        12546         (numeric ID only)

    Returns None if the token cannot be parsed as an SGB identifier.
    """
    token = token.strip()
    m = re.search(r"SGB(\d+)", token, re.IGNORECASE)
    if m:
        return "SGB" + m.group(1)
    if re.fullmatch(r"\d+", token):
        return "SGB" + token
    return None


# ── SAM header index cache ───────────────────────────────────────────────────
# Scanning 11M+ @SQ entries in Python takes ~12 s on every run.  After the
# first scan we write a small gzip-JSON sidecar (<sam>.anubis_idx.json.gz)
# that maps each SGB to its list of reference dicts.  Subsequent runs load
# this cache and skip the Python scan entirely (pysam still parses the header
# in C, but that is much faster).
#
# The cache stores only the SGBs that have been requested so far; it grows
# incrementally.  It is invalidated automatically when the SAM file changes
# (mtime or size mismatch).

def _sam_cache_path(sam_bz2: str) -> str:
    return sam_bz2 + ".anubis_idx.json.gz"


def _load_sam_cache(sam_bz2: str) -> tuple[dict, dict] | None:
    """
    Return (base_header, sgb_index) if a valid cache exists, else None.
    base_header  — non-SQ header dict (@HD, @RG, @PG …)
    sgb_index    — {sgb_id: [{"SN": ..., "LN": ...}, …]}
    """
    path = _sam_cache_path(sam_bz2)
    if not os.path.exists(path):
        return None
    try:
        stat = os.stat(sam_bz2)
        with gzip.open(path, "rt") as f:
            data = json.load(f)
        if data.get("mtime") != stat.st_mtime or data.get("size") != stat.st_size:
            log.info("SAM index cache is stale — will rebuild")
            return None
        return data["base_header"], data["index"]
    except Exception as e:
        log.warning("Could not read SAM index cache (%s): %s", path, e)
        return None


def _save_sam_cache(
    sam_bz2: str,
    base_header: dict,
    new_index: dict,
    existing_index: dict,
) -> None:
    """Merge new_index into existing_index and write the cache sidecar."""
    path = _sam_cache_path(sam_bz2)
    merged = {**existing_index, **new_index}
    stat = os.stat(sam_bz2)
    try:
        with gzip.open(path, "wt") as f:
            json.dump(
                {
                    "version": 1,
                    "mtime": stat.st_mtime,
                    "size": stat.st_size,
                    "base_header": base_header,
                    "index": merged,
                },
                f,
            )
        log.info("SAM index cache written: %s", os.path.basename(path))
    except Exception as e:
        log.warning("Could not write SAM index cache: %s", e)


# ── SAM streaming → per-SGB BAM files ───────────────────────────────────────

def _extract_sgb_from_ref(ref_name: str) -> str | None:
    """Return 'SGB12345' if the reference name encodes an SGB, else None."""
    m = re.search(r"\|SGB(\d+)$", ref_name)
    return ("SGB" + m.group(1)) if m else None


def split_sam_by_sgb(
    sam_bz2: str,
    target_sgbs: set,
    tmpdir: str,
    threads: int = 1,
) -> dict:
    """
    Stream the bz2-compressed MetaPhlAn SAM file once.  For each target SGB,
    collect all reads that mapped to its marker genes and write a sorted,
    indexed BAM file.

    Returns
    -------
    dict  sgb_id → path_to_sorted_indexed_bam
    """
    log.info("Opening bz2 SAM stream …")
    proc = subprocess.Popen(
        ["bzcat", sam_bz2],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )

    in_sam = pysam.AlignmentFile(proc.stdout, "r")

    # ── Build per-SGB sub-headers (cache-backed) ─────────────────────────
    sgb_sq: dict[str, list] = {s: [] for s in target_sgbs}

    cached = _load_sam_cache(sam_bz2)
    cached_base_header, existing_index = cached if cached else (None, {})
    missing_from_cache = target_sgbs - set(existing_index.keys())

    if not missing_from_cache:
        log.info("SAM index cache hit — skipping header scan")
        base_header = cached_base_header
        for sgb in target_sgbs:
            sgb_sq[sgb] = list(existing_index[sgb])
    else:
        log.info("Parsing SAM header for target SGBs …")
        full_header = in_sam.header.to_dict()
        base_header = {k: v for k, v in full_header.items() if k != "SQ"}
        new_index: dict[str, list] = {}
        for sq in tqdm(full_header.get("SQ", []), desc="Header SQ records", unit=" refs"):
            sn = sq["SN"]
            sgb = _extract_sgb_from_ref(sn)
            if sgb and sgb in sgb_sq:
                sgb_sq[sgb].append(sq)
                new_index.setdefault(sgb, []).append(dict(sq))
        _save_sam_cache(sam_bz2, base_header, new_index, existing_index)

    # Only keep SGBs that actually have marker genes in this SAM
    present = {s for s, sqs in sgb_sq.items() if sqs}
    missing = target_sgbs - present
    if missing:
        log.warning(
            "No marker genes found in SAM for: %s", ", ".join(sorted(missing))
        )

    if not present:
        in_sam.close()
        proc.wait()
        return {}

    # ── Open temporary unsorted BAM writers ─────────────────────────────
    unsorted_paths: dict[str, str] = {}
    writers: dict[str, pysam.AlignmentFile] = {}

    for sgb in present:
        header = {**base_header, "SQ": sgb_sq[sgb]}
        unsorted_path = os.path.join(tmpdir, f"{sgb}.unsorted.bam")
        writers[sgb] = pysam.AlignmentFile(unsorted_path, "wb", header=header)
        unsorted_paths[sgb] = unsorted_path

    # ── Build reference-ID remap table ──────────────────────────────────
    # BAM stores reference as an integer index into the header's SQ list.
    # Each per-SGB BAM has a small sub-header, so we must translate the
    # original reference_id to the new index before writing each read.
    #
    # old_tid (full header) → (sgb_id, new_tid in per-SGB header)
    log.info("Building reference-ID remap table …")
    ref_id_to_sgb: dict[int, tuple] = {}
    for sgb, sq_list in sgb_sq.items():
        if sgb not in present:
            continue
        for new_idx, sq in enumerate(sq_list):
            old_tid = in_sam.get_tid(sq["SN"])
            if old_tid >= 0:
                ref_id_to_sgb[old_tid] = (sgb, new_idx)

    # ── Stream reads → route by SGB ─────────────────────────────────────
    log.info("Streaming reads → writing per-SGB BAMs …")
    routed = {s: 0 for s in present}

    for read in tqdm(in_sam.fetch(until_eof=True), desc="Reads", unit=" reads"):
        if read.is_unmapped:
            continue
        mapping = ref_id_to_sgb.get(read.reference_id)
        if mapping is None:
            continue
        sgb, new_tid = mapping
        read.reference_id = new_tid
        writers[sgb].write(read)
        routed[sgb] += 1

    for w in writers.values():
        w.close()
    in_sam.close()
    proc.wait()

    for s, n in routed.items():
        log.info("  %s → %d reads", s, n)

    # ── Sort and index each BAM (parallel across SGBs) ──────────────────
    log.info("Sorting and indexing per-SGB BAMs …")
    sorted_paths: dict[str, str] = {}

    def _sort_and_index(item):
        sgb, unsorted = item
        if routed.get(sgb, 0) == 0:
            return sgb, None
        sorted_bam = os.path.join(tmpdir, f"{sgb}.bam")
        subprocess.run(
            ["samtools", "sort", "-o", sorted_bam, unsorted],
            check=True, stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            ["samtools", "index", sorted_bam],
            check=True, stderr=subprocess.DEVNULL,
        )
        os.remove(unsorted)
        return sgb, sorted_bam

    with ThreadPoolExecutor(max_workers=threads) as pool:
        for sgb, sorted_bam in pool.map(_sort_and_index, unsorted_paths.items()):
            if sorted_bam is None:
                log.warning("  Skipping %s: no reads routed", sgb)
            else:
                sorted_paths[sgb] = sorted_bam
                log.info("  %s → %s", sgb, sorted_bam)

    return sorted_paths


# ── PyDamage per-SGB analysis ────────────────────────────────────────────────

def run_pydamage_grouped(
    bam_path: str,
    wlen: int = 30,
    processes: int = 1,
) -> dict | None:
    """
    Run pydamage in *grouped* mode: all reads in the BAM are treated as
    a single reference (aggregating across all marker genes of the SGB).

    Returns (result_dict, read_dict), or (None, {}) if fitting fails.
    read_dict maps read names to their damaged-base positions and is needed
    for base-quality rescaling.
    """
    try:
        filt_res, read_dict = test_damage(
            ref=None,       # group mode
            bam=bam_path,
            mode="rb",
            wlen=wlen,
            g2a=True,
            subsample=None,
            show_al=False,
            process=processes,
            verbose=False,
        )
    except Exception as e:
        log.warning("  pydamage failed: %s", e)
        return None, {}
    return filt_res, read_dict


def analyze_all_sgbs(
    sgb_bams: dict,
    wlen: int = 30,
    processes: int = 1,
) -> tuple:
    """
    Run pydamage on each per-SGB BAM.

    Returns
    -------
    pd.DataFrame
        One row per SGB with all damage metrics.
    dict
        sgb_id → read_dict (read name → damaged-base positions array),
        needed for base-quality rescaling.
    """
    records = []
    read_dicts: dict = {}

    # Parallelise across SGBs — pydamage group mode runs over a single
    # reference per BAM, so its internal `process` parameter gives no
    # benefit; better to run multiple SGBs concurrently instead.
    n_workers = min(processes, len(sgb_bams))
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(run_pydamage_grouped, bam, wlen, 1): sgb
            for sgb, bam in sgb_bams.items()
        }
        for fut in tqdm(
            as_completed(futures), total=len(futures), desc="PyDamage", unit=" SGBs"
        ):
            sgb = futures[fut]
            try:
                res, read_dict = fut.result()
            except Exception as e:
                log.warning("  %s: pydamage error: %s", sgb, e)
                continue

            if not res:
                log.warning("  %s: model fitting failed, skipping", sgb)
                continue

            read_dicts[sgb] = read_dict
            record = {"sgb_id": sgb}
            for key, val in res.items():
                if key in ("model_params", "residuals"):
                    continue
                record[key] = val
            records.append(record)

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)

    # Rename pydamage internal keys to cleaner column names
    rename_map = {
        "p0":    "null_model_p0",
        "p0_stdev": "null_model_p0_stdev",
        "p":     "damage_model_p",
        "p_stdev": "damage_model_p_stdev",
        "pmin":  "damage_model_pmin",
        "pmin_stdev": "damage_model_pmin_stdev",
        "pmax":  "damage_model_pmax",
        "pmax_stdev": "damage_model_pmax_stdev",
    }
    df.rename(columns=rename_map, inplace=True)

    # FDR correction across all SGBs
    from statsmodels.stats.multitest import multipletests
    pvals = df["pvalue"].dropna()
    if len(pvals) > 0:
        _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
        df.loc[pvals.index, "qvalue"] = qvals

    # Predicted accuracy (GLM)
    prep = df[["coverage", "reflen", "damage_model_pmax"]].rename(
        columns={"damage_model_pmax": "damage", "reflen": "contiglength", "coverage": "actual_cov"}
    ).astype({"actual_cov": float, "contiglength": float, "damage": float})
    df["predicted_accuracy"] = glm_predict(prep, glm_model_params)["predicted_accuracy"]

    return df, read_dicts


# ── Per-SGB BAM processing: rescaling and/or terminal masking ────────────────

def process_sgb_bams(
    sgb_bams: dict,
    result_df: pd.DataFrame | None,
    read_dicts: dict | None,
    wlen: int,
    rescale_threshold: float,
    rescale_alpha: float,
    mask_5p: int,
    mask_3p: int,
    outdir: str,
    threads: int = 1,
) -> dict:
    """
    Apply rescaling and/or terminal masking directly to the per-SGB BAMs that
    were already produced by split_sam_by_sgb().  This avoids a second full
    scan of the original (huge) bz2 SAM — the per-SGB BAMs are small and fast
    to process.

    For each SGB a sorted, indexed BAM is written to <outdir>/<SGB>.bam.
    These are directly consumable by StrainPhlAn's sample2markers.py.

    rescaling  — activated when result_df and read_dicts are provided; only
                 SGBs whose predicted_accuracy ≥ rescale_threshold AND
                 q-value ≤ rescale_alpha are rescaled.
    masking    — activated when mask_5p > 0 or mask_3p > 0; applied to every
                 mapped read regardless of SGB damage status.

    Returns {sgb_id: path_to_processed_bam}.
    """
    from pydamage.rescale import rescale_qual
    from pydamage.models import damage_model as DamageModel
    from array import array as carray

    do_rescale = result_df is not None and read_dicts is not None
    do_mask    = mask_5p > 0 or mask_3p > 0

    os.makedirs(outdir, exist_ok=True)

    # ── Pre-compute per-SGB damage PMFs ──────────────────────────────────
    sgb_pmf:   dict = {}
    sgb_reads: dict = {}

    if do_rescale:
        x = np.arange(wlen)
        for _, row in result_df.iterrows():
            sgb  = row["sgb_id"]
            acc  = row.get("predicted_accuracy")
            qval = row.get("qvalue")
            if acc is None or pd.isna(acc) or qval is None or pd.isna(qval):
                continue
            if float(acc) >= rescale_threshold and float(qval) <= rescale_alpha:
                p    = row.get("damage_model_p")
                pmin = row.get("damage_model_pmin")
                pmax = row.get("damage_model_pmax")
                if any(v is None or pd.isna(v) for v in (p, pmin, pmax)):
                    continue
                sgb_pmf[sgb]   = DamageModel().fit(
                    x, float(p), float(pmin), float(pmax), wlen=wlen
                )
                sgb_reads[sgb] = read_dicts.get(sgb, {}).get("reference", {})

        log.info(
            "Rescaling: %d/%d SGBs pass thresholds "
            "(predicted_accuracy ≥ %.2f, q ≤ %.2f)",
            len(sgb_pmf), len(result_df), rescale_threshold, rescale_alpha,
        )

    if do_mask:
        log.info(
            "Terminal masking: %d bases at 5' end, %d bases at 3' end",
            mask_5p, mask_3p,
        )

    # ── Process each per-SGB BAM ──────────────────────────────────────────
    def _process_one(sgb: str, bam_in: str) -> tuple:
        bam_out = os.path.join(outdir, f"{sgb}.bam")
        n_rescaled = 0
        n_masked   = 0

        with pysam.AlignmentFile(bam_in, "rb") as f_in, \
             pysam.AlignmentFile(bam_out, "wb", header=f_in.header) as f_out:

            for read in f_in.fetch(until_eof=True):
                if not read.is_unmapped:
                    orig_qual = read.query_qualities
                    if orig_qual is not None:
                        qual    = np.array(orig_qual, dtype=np.uint8)
                        changed = False

                        # 1. Model-based rescaling
                        if do_rescale and sgb in sgb_pmf:
                            positions = sgb_reads[sgb].get(read.query_name)
                            if positions is not None:
                                qual = np.array(
                                    rescale_qual(qual, sgb_pmf[sgb], positions,
                                                 reverse=read.is_reverse),
                                    dtype=np.uint8,
                                )
                                n_rescaled += 1
                                changed = True

                        # 2. Terminal masking
                        if do_mask:
                            rlen = len(qual)
                            m5   = min(mask_5p, rlen)
                            m3   = min(mask_3p, rlen - m5)
                            if m5 > 0:
                                qual[:m5] = 0
                            if m3 > 0:
                                qual[-m3:] = 0
                            n_masked += 1
                            changed = True

                        if changed:
                            read.query_qualities = carray("B", qual)

                f_out.write(read)

        # Input BAM is already sorted; just index the output
        subprocess.run(
            ["samtools", "index", bam_out],
            check=True, stderr=subprocess.DEVNULL,
        )
        return sgb, bam_out, n_rescaled, n_masked

    processed: dict = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(_process_one, sgb, bam): sgb
            for sgb, bam in sgb_bams.items()
        }
        for fut in tqdm(
            as_completed(futures), total=len(futures),
            desc="Processing BAMs", unit=" SGBs",
        ):
            sgb = futures[fut]
            try:
                sgb, bam_out, n_rescaled, n_masked = fut.result()
                processed[sgb] = bam_out
                log.info(
                    "  %s: %d model-rescaled, %d terminal-masked → %s",
                    sgb, n_rescaled, n_masked, os.path.basename(bam_out),
                )
            except Exception as e:
                log.warning("  %s: processing failed: %s", sgb, e)

    log.info("Processed BAMs written to: %s", outdir)
    return processed


# ── Damage profile plots ─────────────────────────────────────────────────────

def plot_damage_profile(row: pd.Series, wlen: int, outdir: str) -> None:
    """
    Classic aDNA damage profile plot for one SGB.

    Left panel  – C→T substitution frequency from the 5' end.
    Right panel – G→A substitution frequency from the 3' end (x-axis mirrored).

    Both panels overlay:
      • Observed data as a filled curve (coloured area + line)
      • Fitted PyDamage damage model
      • Null model (flat baseline)
    """
    from pydamage.models import damage_model as DamageModel, null_model as NullModel

    sgb     = row["sgb_id"]
    species = row.get("species", sgb)

    x  = np.arange(wlen)
    ct = np.array([row.get(f"CtoT-{i}", np.nan) for i in x], dtype=float)
    ga = np.array([row.get(f"GtoA-{i}", np.nan) for i in x], dtype=float)

    # Fitted model parameters
    p    = row.get("damage_model_p",    np.nan)
    pmin = row.get("damage_model_pmin", np.nan)
    pmax_val = row.get("damage_model_pmax", np.nan)
    p0   = row.get("null_model_p0",     np.nan)

    # Confidence-band bounds (±2 σ, clamped to [0, 1])
    pmin_sd  = row.get("damage_model_pmin_stdev", 0) or 0
    pmax_sd  = row.get("damage_model_pmax_stdev", 0) or 0

    def _model_curve(p, pmin, pmax):
        if any(np.isnan(v) for v in (p, pmin, pmax)):
            return None
        return DamageModel().fit(x, p, pmin, pmax, wlen=wlen)

    y_fit  = _model_curve(p, pmin, pmax_val)
    y_low  = _model_curve(p, max(0, pmin - 2 * pmin_sd),
                          max(0, pmax_val - 2 * pmax_sd))
    y_high = _model_curve(p, min(1, pmin + 2 * pmin_sd),
                          min(1, pmax_val + 2 * pmax_sd))
    y_null = NullModel().fit(x, p0) if not np.isnan(p0 or np.nan) else None

    pmax_disp = row.get("damage_model_pmax", np.nan)
    qval      = row.get("qvalue", np.nan)
    acc       = row.get("predicted_accuracy", np.nan)
    n_reads   = row.get("nb_reads_aligned", np.nan)
    cov       = row.get("coverage", np.nan)

    def _fmt(v, decimals=3):
        return f"{v:.{decimals}f}" if (v is not None and not np.isnan(v)) else "n/a"

    title = (
        f"{species}  [{sgb}]\n"
        f"reads={int(n_reads) if (n_reads is not None and not np.isnan(n_reads)) else 'n/a'}   "
        f"coverage={_fmt(cov, 2)}×   "
        f"pmax={_fmt(pmax_disp)}   q={_fmt(qval)}   acc={_fmt(acc)}"
    )

    # ── shared y-axis upper limit ────────────────────────────────────────
    y_upper = max(
        np.nanmax(ct) if not np.all(np.isnan(ct)) else 0,
        np.nanmax(ga) if not np.all(np.isnan(ga)) else 0,
        (pmax_disp or 0) * 1.1,
        0.05,
    ) * 1.15

    fig, (ax_ct, ax_ga) = plt.subplots(
        1, 2, figsize=(12, 4.5), sharey=True,
        gridspec_kw={"wspace": 0.05}
    )
    fig.suptitle(title, fontsize=11, y=1.01)

    # ── helper: draw one panel ───────────────────────────────────────────
    def _draw_panel(ax, y_obs, colour, end_label):
        # Observed: filled area + solid line
        ax.fill_between(x, 0, y_obs, color=colour, alpha=0.20)
        ax.plot(x, y_obs, color=colour, linewidth=1.8,
                marker="o", markersize=3, label="Observed")

        # Null model
        if y_null is not None:
            ax.plot(x, y_null, color="grey", linewidth=1.2,
                    linestyle="--", alpha=0.8, label="Null model")

        # Damage model + CI band
        if y_fit is not None:
            if y_low is not None and y_high is not None:
                ax.fill_between(x, y_low, y_high,
                                color="#D7880F", alpha=0.15)
            ax.plot(x, y_fit, color="#D7880F", linewidth=2.0,
                    linestyle="-", alpha=0.9, label="Damage model")

        ax.set_ylim(0, y_upper)
        ax.set_xlim(-0.5, wlen - 0.5)
        ax.set_xlabel(f"Position from {end_label} end", fontsize=10)
        ax.tick_params(axis="both", labelsize=8)
        ax.grid(axis="y", linestyle=":", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)

    # Left panel: C→T from 5' end (x increases left → right)
    _draw_panel(ax_ct, ct, "#bd0d45", "5'")
    ax_ct.set_ylabel("Substitution frequency", fontsize=10)
    ax_ct.set_title("C→T (5' end)", fontsize=10)
    ax_ct.legend(fontsize=8, loc="upper right")
    ax_ct.set_xticks(np.arange(0, wlen, 5))

    # Right panel: G→A from 3' end (x increases right → left)
    _draw_panel(ax_ga, ga, "#2166ac", "3'")
    ax_ga.set_title("G→A (3' end)", fontsize=10)
    ax_ga.invert_xaxis()
    ax_ga.legend(fontsize=8, loc="upper left")
    ax_ga.set_xticks(np.arange(0, wlen, 5))
    ax_ga.yaxis.set_tick_params(labelleft=False)

    safe_name = re.sub(r"[^\w\-]", "_", sgb)
    fig.savefig(os.path.join(outdir, f"{safe_name}.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)


# ── CLI ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="anubis",
        description="Anubis: integrate MetaPhlAn relative abundances with ancient-DNA "
                    "damage analysis (PyDamage) per SGB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-m", "--metaphlan", required=True,
                   help="MetaPhlAn output TSV (with # header lines)")
    p.add_argument("-s", "--sam", required=True,
                   help="MetaPhlAn SAM file (bz2-compressed)")
    p.add_argument("-o", "--outdir", default="metaphlan2damage_out",
                   help="Output directory")

    # Filters
    filt = p.add_argument_group("Filters")
    filt.add_argument("-c", "--clade", "--clades", dest="clades",
                      nargs="+", metavar="CLADE",
                      help="Restrict analysis to specific clade(s). Accepts one or more "
                           "identifiers separated by spaces or commas. Recognised formats: "
                           "SGB12546, t__SGB12546, or plain 12546. When given, "
                           "--min-abundance and --top-n are ignored.")
    filt.add_argument("--min-reads", type=int, default=None,
                      help="Minimum number of reads per SGB (applied after damage analysis)")
    filt.add_argument("--min-abundance", type=float, default=None,
                      help="Minimum relative abundance %% (e.g. 0.01)")
    filt.add_argument("--top-n", type=int, default=None,
                      help="Only analyse the top N most abundant SGBs")

    # Pydamage options
    dam = p.add_argument_group("Damage analysis")
    dam.add_argument("--wlen", type=int, default=30,
                     help="Window length for damage modelling (bp)")
    dam.add_argument("--threads", type=int, default=4,
                     help="CPU threads for sorting / pydamage")

    # Output options
    out = p.add_argument_group("Output")
    out.add_argument("--plot", action="store_true",
                     help="Save per-SGB damage profile plots to <outdir>/plots/")
    out.add_argument("--rescale", action="store_true",
                     help="Rescale base-quality scores using the pydamage damage "
                          "model; writes sorted+indexed BAMs to <outdir>/rescaled/")
    out.add_argument("--rescale-threshold", type=float, default=0.5, metavar="FLOAT",
                     help="Minimum predicted_accuracy for a BAM to be rescaled "
                          "(SGBs below this are written unchanged)")
    out.add_argument("--rescale-alpha", type=float, default=0.05, metavar="FLOAT",
                     help="Maximum q-value for a BAM to be rescaled")
    out.add_argument("--mask-5prime", type=int, default=0, metavar="N",
                     help="Set base quality to 0 for the first N bases of every "
                          "mapped read (5'-terminal masking, independent of --rescale)")
    out.add_argument("--mask-3prime", type=int, default=0, metavar="N",
                     help="Set base quality to 0 for the last N bases of every "
                          "mapped read (3'-terminal masking, independent of --rescale)")
    out.add_argument("--keep-bams", action="store_true",
                     help="Keep intermediate per-SGB BAM files")

    return p


# ── Main ─────────────────────────────────────────────────────────────────────

def main(argv=None):
    args = build_parser().parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    # ── 1. Parse MetaPhlAn table ──────────────────────────────────────────
    log.info("Parsing MetaPhlAn table: %s", args.metaphlan)
    mpa_df = parse_metaphlan(args.metaphlan)
    log.info("  Found %d SGB-level entries", len(mpa_df))

    if args.clades:
        # Flatten comma-separated tokens so both of these work:
        #   -c SGB12546 t__SGB12547
        #   -c SGB12546,t__SGB12547
        raw_tokens: list[str] = []
        for tok in args.clades:
            raw_tokens.extend(tok.split(","))

        norm: dict[str, str] = {}  # canonical SGB id → original input token
        for tok in raw_tokens:
            sgb = _normalize_clade(tok)
            if sgb:
                norm[sgb] = tok
            else:
                log.warning("Cannot parse clade specifier '%s' — skipping", tok)

        if not norm:
            log.error("No valid clade identifiers could be parsed. Exiting.")
            sys.exit(1)

        mpa_df = mpa_df[mpa_df["sgb_id"].isin(norm)].reset_index(drop=True)

        missing_clades = set(norm) - set(mpa_df["sgb_id"])
        if missing_clades:
            log.warning(
                "Requested clade(s) not found in MetaPhlAn table: %s",
                ", ".join(sorted(missing_clades)),
            )

        log.info(
            "  Clade filter: %d/%d SGBs found in MetaPhlAn table",
            len(mpa_df), len(norm),
        )
    else:
        mpa_df = apply_filters(mpa_df, args.min_abundance, args.top_n)
        log.info("  After filters: %d SGBs selected", len(mpa_df))

    if mpa_df.empty:
        log.error("No SGBs passed the filters. Exiting.")
        sys.exit(1)

    target_sgbs = set(mpa_df["sgb_id"])
    log.info("  Target SGBs: %s", ", ".join(sorted(target_sgbs)))

    # ── 2. Split SAM → per-SGB BAMs ──────────────────────────────────────
    tmpdir = tempfile.mkdtemp(prefix="anubis_", dir=args.outdir)
    try:
        sgb_bams = split_sam_by_sgb(
            args.sam,
            target_sgbs,
            tmpdir,
            threads=args.threads,
        )

        if not sgb_bams:
            log.error("No reads could be routed to any target SGB. Exiting.")
            sys.exit(1)

        # ── 3. Run pydamage ───────────────────────────────────────────────
        log.info("Running PyDamage on %d SGB BAMs …", len(sgb_bams))
        damage_df, read_dicts = analyze_all_sgbs(
            sgb_bams,
            wlen=args.wlen,
            processes=args.threads,
        )

        if damage_df.empty:
            log.error("PyDamage produced no results. Exiting.")
            sys.exit(1)

        # ── 4. Merge and filter early so rescaling uses the final result ──
        merged_early = mpa_df.merge(damage_df, on="sgb_id", how="inner")
        if args.min_reads is not None:
            merged_early = merged_early[
                merged_early["nb_reads_aligned"] >= args.min_reads
            ]

        # ── 5. Rescale/mask per-SGB BAMs (no second SAM scan) ────────────
        _mask_5p = args.mask_5prime or 0
        _mask_3p = args.mask_3prime or 0
        _do_process = args.rescale or _mask_5p > 0 or _mask_3p > 0
        if _do_process and not merged_early.empty:
            rescale_dir = os.path.join(args.outdir, "rescaled")
            process_sgb_bams(
                sgb_bams=sgb_bams,
                result_df=merged_early if args.rescale else None,
                read_dicts=read_dicts   if args.rescale else None,
                wlen=args.wlen,
                rescale_threshold=args.rescale_threshold,
                rescale_alpha=args.rescale_alpha,
                mask_5p=_mask_5p,
                mask_3p=_mask_3p,
                outdir=rescale_dir,
                threads=args.threads,
            )

    finally:
        if not args.keep_bams:
            shutil.rmtree(tmpdir, ignore_errors=True)
        else:
            log.info("Per-SGB BAMs kept in: %s", tmpdir)

    merged = merged_early
    if merged.empty:
        log.error("No SGBs remain after filters. Exiting.")
        sys.exit(1)

    if args.min_reads is not None:
        log.info(
            "min-reads filter (%d): %d SGBs kept",
            args.min_reads, len(merged),
        )

    # ── 6. Select and order output columns ───────────────────────────────
    ct_cols = [f"CtoT-{i}" for i in range(args.wlen) if f"CtoT-{i}" in merged.columns]
    ga_cols = [f"GtoA-{i}" for i in range(args.wlen) if f"GtoA-{i}" in merged.columns]

    core_cols = [
        "sgb_id", "species", "relative_abundance",
        "nb_reads_aligned", "coverage", "reflen",
        "damage_model_pmax", "damage_model_pmax_stdev",
        "damage_model_p",
        "damage_model_pmin", "damage_model_pmin_stdev",
        "null_model_p0",
        "pvalue", "qvalue",
        "predicted_accuracy",
        "RMSE",
    ]
    core_cols = [c for c in core_cols if c in merged.columns]
    out_cols = core_cols + ct_cols + ga_cols

    result = merged[out_cols].copy()
    result = result.round(
        {c: 4 for c in result.select_dtypes("float").columns}
    )
    result.sort_values("relative_abundance", ascending=False, inplace=True)

    # ── 7. Write output table ─────────────────────────────────────────────
    out_tsv = os.path.join(args.outdir, "anubis_results.tsv")
    result.to_csv(out_tsv, sep="\t", index=False)
    log.info("Results written to: %s", out_tsv)

    # ── 8. Print summary ──────────────────────────────────────────────────
    summary_cols = [
        "sgb_id", "species", "relative_abundance",
        "nb_reads_aligned", "coverage",
        "damage_model_pmax", "qvalue", "predicted_accuracy",
    ]
    summary_cols = [c for c in summary_cols if c in result.columns]
    print("\n" + result[summary_cols].to_string(index=False) + "\n")

    # ── 9. Damage profile plots ──────────────────────────────────────────
    if args.plot:
        plotdir = os.path.join(args.outdir, "plots")
        os.makedirs(plotdir, exist_ok=True)
        log.info("Writing damage plots to: %s", plotdir)
        for _, row in result.iterrows():
            plot_damage_profile(row, wlen=args.wlen, outdir=plotdir)
        log.info("  %d plots written", len(result))

    return result


if __name__ == "__main__":
    main()
