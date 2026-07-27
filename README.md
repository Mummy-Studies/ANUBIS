# ANUBIS
Tool for automated species-level authentication of ancient microbial DNA from MetaPhlAn marker alignments.

## Overview
ANUBIS is an automated pipeline that joins MetaPhlAn profiling with PyDamage authentication. The resulting report integrates microbial abundance, authentication metrics, and damage plots in a single workflow. The tool works at  species-level genome bins and rerturns separate information for each SGB in a single sample. 
The parameters allow the filtering of species based on abundance, marker coverage and species name whereas additional parameters allow to modify the .bam files to rescale and mask damaged regions.

## Features
## Installation
## Usage
### Asking Help
```bash
python anubis.py -h
```
### Essential set-up
```bash
python anubis.py \
            -m input.metaphlan.table \
            -s input.sam.file
```
### Optional parameters
```bash
python anubis.py \
            -m input.metaphlan.table \
            -s input.sam.file \
            --min-reads MIN_READS \
            --min-abundance MIN_ABUNDANCE \
            --top-n TOP_N \
            --wlen WLEN \
            --threads THREADS \
            --plot \
            --rescale \
            --rescale-threshold FLOAT \
            --rescale-alpha FLOAT \
            --mask-5prime N \
            --mask-3prime N \
            --keep-bams
```            
### Example and supplementary material
Examples are provided in the `examples/` directory, including supplementary scripts for merging tables across multiple samples and a workflow demonstrating how to plot damage patterns from multiple files.

### Clone repository
```bash
https://github.com/Mummy-Studies/ANUBIS.git
```
