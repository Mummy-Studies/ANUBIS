# ANUBIS
Tool for automated species-level authentication of ancient microbial DNA from MetaPhlAn marker alignments

## Overview
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

### Clone repository
```bash
https://github.com/Mummy-Studies/ANUBIS.git
```
