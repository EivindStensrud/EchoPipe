# Tutorial: How to use EchoPipe

> Test run on Scandinavian amphibians and batra primers

## Background

EchoPipe is an iterative, reproducible pipeline for creating, curating, evaluating, and reformatting reference databases for environmental DNA (eDNA) metabarcoding studies. The tool allows for creating, refining, and exporting curated marker reference databases suitable for downstream taxonomic assignment.

By utilizing an iterative workflow, computational workload, memory, and runtime are minimized while maximizing target taxonomic coverage. The entire pipeline is distributed as a Python package installable via pip.

---

## Installation & requirements

### Mandatory inputs & dependencies

- **NCBI API Key**: https://support.nlm.nih.gov/kbArticle/?pn=KA-05317 — increases Entrez rate limits (recommended).
- **Species list** with scientific names (CSV or TXT).
- **Email address** for NCBI Entrez contact.
- **Forward and reverse primer sequences** spanning the target region (5'→3').
- **Python 3.9+** (the tutorial and package were tested with Python 3.10).

### Step-by-step installation

1. Create and activate a virtual environment (example using conda):

```bash
conda create -n echopipe python=3.10 -y
conda activate echopipe
```

2. Install EchoPipe from GitHub via pip:

```bash
pip install git+https://github.com/EivindStensrud/EchoPipe.git
```

3. Verify the installation:

```bash
echopipe --help
```


![Figure 1: EchoPipe workflow](https://github.com/user-attachments/assets/d2088ad1-da3c-4aa4-b14c-b1fa9c9e3785)

_Figure 1. EchoPipe modular workflow for iterative database creation and curation._

---

## Overview of stages

EchoPipe is organized into modular stages (subcommands):

- template — generate an initial (uncurated) reference template using primer matching.
- create — mine NCBI for candidate sequences and extract marker regions via local BLAST.
- curate — align extracted marker regions, apply length/quality filters, and flag problematic sequences.
- complete — build the official database, merging curated content and removing flagged accessions.
- evaluate — check completeness, monophyly, primer binding, and other diagnostics.
- reformat — export reference databases in formats compatible with common taxonomic classifiers.

---

## Stage 1 — Generate reference template (echopipe template)

The `template` subcommand generates an initial, uncurated marker reference database using primer alignment matching.

Common options (template):

```text
input_file            # TXT/CSV with species names, or a FASTA with provided sequences
-f, --forward         # Forward primer (5'-3')
-r, --reverse         # Reverse primer (5'-3')
-e, --email           # Contact email for NCBI Entrez
-a, --api_key         # NCBI API key
-q, --query           # Custom NCBI query additions
-t, --threshold       # Minimum sequence length (default: 150)
-l, --length          # Maximum sequence length (default: 22000)
-m, --max             # Number of sequences downloaded per species (default: 1)
-p, --provided_sequences  # Use a FASTA file as reference template
-z, --longest_amplicon_size # Multiplier for median length (default: 2)
-n, --random_subset   # Number of random species to use
-sf, --subset_file    # Path to file with a specific subset of species
-T, --threads         # Parallel threads (default: auto-detected, max 7)

-C, --Complete        # Complete the reference template database (used during finalization)
```

Example: build an initial unaligned template

```bash
echopipe template Norwegian_swedish_amphibians_250214.csv \
  -f ACACCGCCCGTCACCCT \
  -r GTAYACTTACCATGTTACGACTT \
  -l 1000 \
  -e your.email@example.com \
  -a YOUR_NCBI_API_KEY \
  -t 40 \
  -m 5 \
  -q "AND 12S"
```

Inspect and clean the template alignment

1. Open `Reference_template_creation/aligned_sequences_to_curate.fasta` in an MSA viewer (e.g., Jalview).
2. Color by nucleotide to view conserved regions.
3. Verify sequences overlap forward and reverse primer binding sites.
4. Remove obviously dubious or non-target sequences.
5. Save the file without renaming or moving it (the pipeline expects the same path).

Complete the reference template (finalize)

```bash
echopipe template unique_Norwegian_swedish_amphibians_250214.csv \
  -e your.email@example.com \
  -a YOUR_NCBI_API_KEY \
  -C
```

---

## Stage 2 — Create reference database (echopipe create)

The `create` subcommand mines NCBI for candidate sequences, extracts local marker regions via local BLAST against the template database, and includes automatic retry logic for network disruptions.

Key options (create):

```text
input_file      # TXT/CSV species list
input_database  # Path to reference template FASTA
-e, --email     # Email address
-a, --api_key   # NCBI API key
-s, --sort      # Sort by length (not recommended)
-c, --maxcount  # Max accession numbers per species (default: 10000)
-l, --maxlength # Max sequence length allowed (default: 22000)
-z, --ampliconsize # Minimum amplicon size (default: 50)
-m, --mitochondria  # Target mitochondrial sequences
-r, --ribosomal  # Target mitochondrial 12S rDNA
-q, --query     # Custom NCBI query
-b, --batch_size# Batch size for downloads (default: 5000)
-t, --taxid     # Use saved taxid list
-E, --evalue    # E-value for BLAST (default: 20 => 5e-20)
-R, --repeat    # Repeat curation on previously downloaded sequences
-T, --threads   # Parallel threads (default: auto-detected, max 7)
```

Tip (E-value selection): keep the E-value lower for shorter markers (e.g., 15) and increase for longer markers if needed.

Example: standard database creation

```bash
echopipe create unique_Norwegian_swedish_amphibians_250214.csv reference_template_database.fasta \
  -e your.email@example.com \
  -a YOUR_NCBI_API_KEY \
  -E 20 \
  -b 5000 \
  -T 4
```

Output database header format

Dereplicated sequences track abundance using an appended sequence counter. Example header and sequence:

```text
>gb|Accession_number|NCBI_taxonomy|counter
SEQUENCE_DATA
```

Example:

```text
>gb|KJ858774.1|Eukaryota;Chordata;Amphibia;Anura;Alytidae;Alytes;Alytes_obstetricans|2
ACACCGCCCGTCACCCTCCTCAACTAACTCAACCCCCTAACTAAAAGCTAACTGGTTAACAAGAAGAGGCAAGTCGTAACATGGTAAGTATA
```

---

## Stage 3 — Curate database (echopipe curate)

The `curate` subcommand aligns extracted marker regions, applies strict length/ambiguity filters, builds FastTree gene trees, and flags potential misannotations or paralogs.

Options (curate):

```text
input_file          # Database to revise (FASTA)
-o, --old_database  # Previous database version
-N, --number_ns     # Allowed number of N/ambiguous nucleotides (default: 0)
-M, --mafft_online  # Path to MAFFT online alignment file
--min_length        # Minimum sequence length to keep (default: 150)
--max_length        # Maximum sequence length to keep
```

Example execution with target bounds

```bash
echopipe curate BLAST_results/2026-03-17_1_to_curate.fasta \
  --min_length 200 \
  --max_length 600 \
  -N 0
```

Diagnostic curation output (examples):

- `Database_curation/{date}/{date}_{run}_aligned.fasta` — aligned dataset (* marks new sequences)
- `Database_curation/{date}/{date}_{run}_tree_string.newick` — gene tree for visualization
- `Database_curation/{date}/{date}_{run}_duplicate_sequences.txt` — identical sequences across species
- `Database_curation/{date}/{date}_{run}_sequences_to_delete.txt` — user file to list misannotated accessions to exclude

Options used by the `complete` subcommand

```text
-b, --blast_file     # New database BLAST file
-c, --curated_file   # Curated aligned FASTA
-o, --old_database   # Existing master database
-u, --updated_database # Output database (default: database.fasta)
```

---

## Stage 4 — Complete official database (echopipe complete)

Once misannotated accession numbers are recorded in `sequences_to_delete.txt`, finalize your official reference database.

Initial finalization example:

```bash
echopipe complete \
  -b BLAST_results/2026-03-17_1_to_curate.fasta \
  -c Database_curation/2026-03-17_1/2026-03-17_1_aligned.fasta \
  -u Amphibian_2026-03-17_1.fasta
```

Updating an existing database (merge new data while preserving counters):

```bash
echopipe complete \
  -b BLAST_results/2026-03-17_2_to_curate.fasta \
  -u Amphibian_2026-03-17_2.fasta \
  -o Amphibian_2026-03-17_1.fasta
```

---

## Stage 5 — Evaluate database (echopipe evaluate)

The `evaluate` subcommand tests database completeness, taxid resolution, GC-content distribution, monophyly, and primer binding compatibility.

Options (evaluate):

```text
reference_database      # Path to reference FASTA
monophyletic_group     # Path to monophyletic groups txt file
-f, --forward_primer   # Forward primer (5'-3')
-r, --reverse_primer   # Reverse primer (5'-3')
```

Example:

```bash
echopipe evaluate Amphibian_2026-03-17_1.fasta \
  Database_curation/2026-03-17_1/Curated_content/2026-03-17_1_post_curation_monophyletic_group.txt \
  -f ACACCGCCCGTCACCCT \
  -r GTAYACTTACCATGTTACGACTT
```

After evaluation, the tool provides primer-mismatch statistics and soft suggestions for reformatting the database for downstream classifiers.

---

## Stage 6 — Reformat for taxonomic classifiers (echopipe reformat)

Convert official EchoPipe references into header formats compatible with common taxonomic classifiers.

Supported formats example (code, target, header):

| Format code | Target classifier         | Example header output (brief) |
|-------------|--------------------------|-------------------------------|
| sintax      | SINTAX                   | >ACC;tax=d:Eukaryota,p:Chordata,c:Amphibia... |
| rdp         | RDP Classifier           | >ACC\troot;Eukaryota;Chordata;Amphibia...      |
| dadt        | DADA2 assignTaxonomy     | >Eukaryota;Chordata;Amphibia;Anura;Alytidae;Alytes |
| dads        | DADA2 assignSpecies      | >ACC Alytes obstetricans                        |
| idt         | IDTAXA                   | >Root;Eukaryota;Chordata;Amphibia;Anura;Alytidae... |
| qiime       | QIIME 2                  | >ACC (+ companion taxonomy .txt file)           |

Example reformat commands:

```bash
# Reformat for QIIME 2 (generates a FASTA and companion 2-column taxonomy .txt)
echopipe reformat Amphibian_2026-03-17_1.fasta qiime

# Reformat for DADA2 assignTaxonomy
echopipe reformat Amphibian_2026-03-17_1.fasta dadt

# Reformat for SINTAX
echopipe reformat Amphibian_2026-03-17_1.fasta sintax
```

---

## Troubleshooting & performance tips

### Handling NCBI connection issues (IncompleteRead / 0 bytes read)

High-volume species or large searches can cause NCBI Entrez HTTP streams to drop. EchoPipe includes adaptive download logic:

- Initial run: uses `-b/--batch_size` (default 5000).
- First retry: re-attempt with same batch size after a 10s pause.
- Subsequent retries: dynamically halves the batch size (2500 → 1250 → 250) to avoid timeouts.

If network failures persist, manually reduce threads and batch size:

```bash
# Lower worker threads (-T) and batch size (-b) for unstable connections
echopipe create species.txt ref_template.fasta -e user@mail.com -a API_KEY -T 2 -b 1000
```

---

If you'd like, I can also:

- Split this tutorial into a README-style quickstart and a longer user guide,
- Add code-syntax highlighting for each command and option table, or
- Create small example species lists and dataset fixtures to reproduce the demo run.
