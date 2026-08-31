# EchoPipe

EchoPipe is an iterative, reproducible pipeline for creating, curating, evaluating and reformatting marker reference databases for environmental DNA (eDNA) metabarcoding.

This repository contains the EchoPipe CLI and examples. For a detailed walkthrough see the full user guide in docs/TUTORIAL.md.

Status: Stable — command-line tool tested with Python 3.10. See the full tutorial for options and troubleshooting.

## Quickstart (recommended minimal steps)

Prerequisites

- Python 3.9+ (3.10 recommended)
- Conda or a Python virtual environment
- NCBI API key (recommended) and a contact email for Entrez

Install

```bash
# create and activate a virtual environment (conda example)
conda create -n echopipe python=3.10 -y
conda activate echopipe

# install EchoPipe from GitHub
pip install git+https://github.com/EivindStensrud/EchoPipe.git

# check CLI
echopipe --help
```

Prepare inputs

- Provide a species list (CSV or TXT) with scientific names.
- Provide primer sequences (forward and reverse, 5'→3').

Quick workflow (replace placeholders with your files / values)

```bash
# 1) Build an initial template (primer sequences required)
echopipe template species_list.csv \
  -f FORWARD_PRIMER \
  -r REVERSE_PRIMER \
  -e your.email@example.com \
  -a YOUR_NCBI_API_KEY

# 2) Inspect and curate the alignment generated under Reference_template_creation/
# (open aligned_sequences_to_curate.fasta in an MSA viewer, remove bad sequences)

# 3) Finalize the template
echopipe template species_list.csv -e your.email@example.com -a YOUR_NCBI_API_KEY -C

# 4) Create the reference database (mines NCBI and extracts regions)
echopipe create species_list.csv reference_template_database.fasta -e your.email@example.com -a YOUR_NCBI_API_KEY

# 5) Curate extracted sequences
echopipe curate BLAST_results/<date>_to_curate.fasta --min_length 200 --max_length 600 -N 0

# 6) Complete and export the official database
echopipe complete -b BLAST_results/<date>_to_curate.fasta -c Database_curation/<date>/<date>_aligned.fasta -u MyDatabase.fasta

# 7) Reformat for downstream classifiers (example: QIIME)
echopipe reformat MyDatabase.fasta qiime
```

Where to find more

- Full user guide and examples: docs/TUTORIAL.md
- If you already have species and primer lists, substitute them into the quickstart commands above.

Contributing

Contributions, bug reports and PRs are welcome. If you want me to add example species lists or helper scripts, tell me which format you prefer and I can add them.
