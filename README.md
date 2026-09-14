# EchoPipe

EchoPipe is an iterative, reproducible pipeline for creating, curating, evaluating and reformatting marker reference databases for environmental DNA (eDNA) metabarcoding.

This repository contains the EchoPipe CLI and examples. For a detailed walkthrough see the full user guide in docs/TUTORIAL.md.

Status: Stable - command-line tool tested with Python 3.10+. See the full tutorial for options and troubleshooting.

## Quickstart (recommended minimal steps)

Prerequisites

- Python 3.10+  
- Conda or a Python virtual environment  
- NCBI API key (recommended) and a contact email for Entrez  

Install

```bash
# create and activate a virtual environment (conda example)
conda create -n echopipe python=3.10 -y
conda activate echopipe

# safely install the older Biopython version using conda-forge 
# (This prevents C-compiler errors on Macs and Windows)
conda install -c conda-forge "biopython<1.81" -y

# install EchoPipe from GitHub
pip install git+https://github.com/EivindStensrud/EchoPipe.git

# check CLI
echopipe --help
```
If you have errors regarding installation, see the Installation Troubleshooting Guide at the end of the README file.  

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

# 3) Finalize the template (Note: Use the "unique_" file generated in Step 1)
echopipe template -C unique_species_list.csv

# 4) Create the reference database (mines NCBI and extracts regions)
echopipe create unique_species_list.csv reference_template_database.fasta \
  -e your.email@example.com \
  -a YOUR_NCBI_API_KEY

# 5) Curate extracted sequences
echopipe curate BLAST_results/<date>_to_curate.fasta --min_length 200 --max_length 600 -N 0

# 6) Complete and export the official database
echopipe complete -b BLAST_results/<date>_to_curate.fasta -c Database_curation/<date>/<date>_aligned.fasta -u MyDatabase.fasta

# 7) Evaluate database quality and completeness
echopipe evaluate MyDatabase.fasta \
  Database_curation/<date>/Curated_content/<date>_post_curation_monophyletic_group.txt \
  -f FORWARD_PRIMER \
  -r REVERSE_PRIMER

# 8) Reformat for downstream classifiers (optional)
echopipe reformat MyDatabase.fasta qiime  
```

---

## Key workflow stages

### Stage 6: Evaluate database (quality assurance)

Before using your reference database, run the **evaluate** step to:

- Test taxonomic **monophyly** (are species forming coherent monophyletic groups?)
- Check **primer binding sites** (do all sequences align with your forward/reverse primers?)
- Assess **GC content** distribution across the database
- Detect **identical sequences** across species (potential misannotations)
- Resolve **NCBI taxonomy IDs** and verify accession annotations
- Generate **diagnostic plots** and summary statistics

Example:

```bash
echopipe evaluate Amphibian_2026-03-17_1.fasta \
  Database_curation/2026-03-17_1/Curated_content/2026-03-17_1_post_curation_monophyletic_group.txt \
  -f ACACCGCCCGTCACCCT \
  -r GTAYACTTACCATGTTACGACTT
```

**Output files** (stored under `Evaluation/<database_name>/`):
- `{database}_GC_content_histogram.png` — GC-content distribution
- `{database}_sequence_summary.csv` — Sequence statistics
- `{database}_evaluation_species_summary.csv` — Species-level metrics
- `{database}_evaluation_family_summary.csv` — Higher-order taxonomic metrics
- Primer mismatch visualizations for forward and reverse primers

This step is essential for **quality control and troubleshooting** before using your database downstream.

---

### Stage 7: Reformat for classifiers (optional)

If you need to export your database for external taxonomic classifiers (QIIME 2, DADA2, SINTAX, RDP, etc.), use the **reformat** step:

```bash
echopipe reformat MyDatabase.fasta qiime    # QIIME 2
echopipe reformat MyDatabase.fasta dadt     # DADA2 assignTaxonomy
echopipe reformat MyDatabase.fasta sintax   # SINTAX
echopipe reformat MyDatabase.fasta rdp      # RDP Classifier
echopipe reformat MyDatabase.fasta dads     # DADA2 assignSpecies
echopipe reformat MyDatabase.fasta idt      # IDTAXA
```

See `docs/TUTORIAL.md` for more details on supported formats.

---

## Where to find more

- **Full user guide with examples**: docs/TUTORIAL.md
- **Installation & requirements**: see Prerequisites above
- **Troubleshooting & performance tips**: docs/TUTORIAL.md § Troubleshooting & performance tips
- If you already have species and primer lists, substitute them into the quickstart commands above.

## Installation Troubleshooting Guide:
### Error 1: "requires a different Python: 3.10.x not in '>=3.11'"  
The Problem:
You will see this error if you are trying to install EchoPipe in an environment running Python 3.10 or older. EchoPipe requires Python 3.11 or higher to ensure compatibility with the underlying bioinformatics libraries (like Biopython).

The Solution:
You need to delete your current environment and create a new one explicitly asking for Python 3.11.

Run these commands one by one:

1. Deactivate your current environment:
```bash
conda deactivate
```
2. Delete the old environment (replace 'echopipe' with your environment name if different):

```bash
conda env remove -n echopipe -y
```
3. Create a brand new environment with Python 3.11:
```bash
conda create -n echopipe python=3.11 -y
```
4. Activate the new environment and install EchoPipe:

```bash
conda activate echopipe
pip install git+https://github.com/EivindStensrud/EchoPipe.git
```

### Error 2: SSL Errors or "Read timed out" during Conda creation
The Problem:
If you are on a strict university/corporate network, a VPN, or a poor internet connection (like a mobile hotspot), Conda might fail to download the required packages and throw an SSL record layer failure or a ReadTimeoutError.

The Solution:
Tell Conda to be more lenient with network timeouts and SSL checks temporarily:

```bash
# Increase the timeout limit
conda config --set remote_read_timeout_secs 120

# (Optional) Disable strict SSL verification if you are behind a firewall
conda config --set ssl_verify false

# Clean the corrupted cache and try again
conda clean -i -y
conda create -n echopipe python=3.11 -y
```

### Error 3: "Command not found: echopipe" after successful installation
The Problem:
Pip successfully installed the tool, but your terminal doesn't recognize the echopipe command.

The Solution:
This usually means you forgot to activate your Conda environment before trying to run the tool. EchoPipe is only accessible inside the environment where it was installed.

Always run this before starting your work:

```bash
conda activate echopipe
echopipe --help
```

