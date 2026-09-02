#!/usr/bin/env python3

import os
import glob
import subprocess
import argparse
import shutil


def run(cmd):
    print("\nRunning:")
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--sample",
    required=True,
    help="Sample name"
)

parser.add_argument(
    "--base-dir",
    required=True,
    help="SGB directory"
)

parser.add_argument(
    "--threads",
    type=int,
    default=16
)

parser.add_argument(
    "--reference_genome",
    required=True,
    help="reference genome file path"
)

parser.add_argument(
    "--output-dir",
    required=True,
    help="Directory where PyDamage results will be written"
)

args = parser.parse_args()
sample = args.sample
base = args.base_dir
threads = args.threads
reference = args.reference_genome # "/path/to/reference_genome.fna"
output_dir = args.output_dir


# -----------------------------
# PATHS
# -----------------------------

sample_dir = os.path.abspath(os.path.join(
    base,
    sample
))

out_dir = os.path.join(
    output_dir,
    "alignment"
)

os.makedirs(out_dir, exist_ok=True)


# -----------------------------
# INDEXING REFERENCE GENOME
# -----------------------------

index_base = os.path.splitext(reference)[0]

if not os.path.exists(index_base + ".1.bt2") and not os.path.exists(index_base + ".1.bt2l"):
    run(
        f"bowtie2-build {reference} {index_base}"
    )

if not os.path.exists(reference + ".fai"):
    run(
        f"samtools faidx {reference}"
    )


# -----------------------------
# LOOKING FOR TRIMMED READS
# -----------------------------

reads = []
for pattern in [
    os.path.join(sample_dir, "*.trimmed.fastq.gz"),
    os.path.join(sample_dir, "*.trimmed.fastq"),
    os.path.join(sample_dir, "**", "*.trimmed.fastq.gz"),
    os.path.join(sample_dir, "**", "*.trimmed.fastq"),
]:
    for path in glob.glob(pattern, recursive=True):
        if os.path.isfile(path):
            reads.append(path)

if not reads:
    for root, dirs, files in os.walk(sample_dir):
        for f in files:
            if f.endswith(".trimmed.fastq.gz") or f.endswith(".trimmed.fastq"):
                reads.append(os.path.join(root, f))

reads = sorted(set(reads))

print(
    f"Sample directory: {sample_dir}\n"
    f"Trimmed read paths: {reads}"
)

bam_files = []


# -----------------------------
# ALIGN EACH RUN
# CREATE .SAM
# CREATE .BAM
# INDEX THE .BAM
# -----------------------------

for fq in reads:

    name = os.path.basename(fq)
    run_name = name.replace(
        ".trimmed.fastq.gz",
        ""
    )

    sam = os.path.join(
        out_dir,
        run_name + ".sam"
    )

    bam = os.path.join(
        out_dir,
        run_name + ".sorted.bam"
    )

    run(
        f"""
        bowtie2 \
        -x {index_base} \
        -U {fq} \
        -p {threads} \
        --very-sensitive \
        --no-unal \
        -S {sam}
        """
    )
    # -x The basename of the index for the reference genome
    # -U Comma-separated list of files containing unpaired reads to be aligned
    # --no-unal Suppress SAM records for reads that failed to align

    run(
        f"""
        samtools view \
        -bSq 30 {sam} | \
        samtools sort \
        -@ {threads} \
        -o {bam}
        """
    )

    run(
        f"""
        samtools index {bam}
        """
    )


    bam_files.append(bam)



# -----------------------------
# MERGE ALL RUNS IN A SINGLE .BAM
# SORT IT
# INDEX IT
# -----------------------------

merged = os.path.join(
    out_dir,
    f"{sample}.merged.bam"
)

run(
    "samtools merge -f "
    f"-@ {threads} "
    f"{merged} "
    + " ".join(bam_files)
)

merged_sorted = os.path.join(
    out_dir,
    f"{sample}.merged.sorted.bam"
)

run(
    f"""
    samtools sort \
    -@ {threads} \
    {merged} \
    -o {merged_sorted}
    """
)

run(
    f"""
    samtools index \
    {merged_sorted}
    """
)


# -----------------------------
# RUN PYDAMAGE
# -----------------------------

pydamage_out = os.path.join(
    output_dir,
    "pydamage"
)

os.makedirs(pydamage_out, exist_ok=True)

# ensure pydamage output directory is fresh (pydamage may abort on existing dir)
if os.path.exists(pydamage_out):
    shutil.rmtree(pydamage_out)
os.makedirs(pydamage_out, exist_ok=True)

run(
    f"""
    pydamage --outdir {pydamage_out} analyze --wlen 30 --process {threads} --group --force --plot {merged_sorted}
    """
)

print(f"Finished {sample}")

