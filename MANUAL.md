# NAME
**echopipe** - automated eDNA reference database generation, curation, evaluation, and reformatting pipeline

# SYNOPSIS
**echopipe** *subcommand* [*options*] [*arguments*]

**echopipe template** *input_file* [**-f** *forward*] [**-r** *reverse*] [**-e** *email*] [**-a** *api_key*] [*options*]

**echopipe create** *input_file* *input_database* [**-e** *email*] [**-a** *api_key*] [*options*]

**echopipe curate** *input_file* [*options*]

**echopipe complete** [**-b** *blast_file*] [**-c** *curated_file*] [**-u** *updated_database*] [*options*]

**echopipe evaluate** *reference_database* *monophyletic_group* [*options*]

**echopipe reformat** *reference_database* *format*

# DESCRIPTION
**EchoPipe** is a comprehensive bioinformatics pipeline designed to create, clean, evaluate, and export custom environmental DNA (eDNA) reference databases. It mines sequence data from NCBI Entrez/BLAST, extracts target regions using local alignment matching, identifies monophyletic and polyphyletic species boundaries, handles network interruptions dynamically, and formats headers for popular taxonomic classifiers.

The pipeline is entirely modular and iterative, ensuring previously validated sequences are preserved (via an abundance counter) rather than re-downloaded.

# SUBCOMMANDS AND OPTIONS

## 1. template
Generates an initial, uncurated marker reference database using primer alignment matching.

**Required arguments for initial generation:**
* **`input_file`**
  Txt or CSV file containing species names, or a FASTA file (requires `-p`).
* **`-f, --forward`** *STRING*
  The forward primer sequence spanning the region of interest (5'-3').
* **`-r, --reverse`** *STRING*
  The reverse primer sequence spanning the region of interest (5'-3').
* **`-e, --email`** *STRING*
  Email address for NCBI Entrez connection.
* **`-a, --api_key`** *STRING*
  NCBI API key to allow faster downloads (10 requests/sec).

**Optional arguments:**
* **`-q, --query`** *STRING*
  Custom NCBI search query additions (e.g., `"AND 12S"` or `"NOT 18S"`).
* **`-t, --threshold`** *INT*
  The minimum sequence length including primers. Shorter sequences are discarded (default: 150).
* **`-l, --length`** *INT*
  The maximum allowed sequence length downloaded from NCBI (default: 22000).
* **`-m, --max`** *INT*
  Maximum number of sequences downloaded per species (default: 1).
* **`-p, --provided_sequences`**
  Flag indicating the input file is an existing FASTA file to be used as a template.
* **`-z, --longest_amplicon_size`** *FLOAT*
  Multiplier for median sequence length to discard abnormally long amplicons (default: 2).
* **`-n, --random_subset`** *INT*
  Number of random species to subset for template creation.
* **`-sf, --subset_file`** *FILE*
  Path to a text file containing a specific subset of species to use.
* **`-T, --threads`** *INT*
  Number of parallel worker threads. Capped at 7 to respect NCBI limits.

**Completion flag:**
* **`-C, --Complete`**
  Finalizes the reference template generation after manual alignment inspection. (Used in a subsequent run).

---

## 2. create
Mines NCBI for candidate sequences, extracts local marker regions via BLAST against the template database, and automatically manages download retries.

**Required arguments:**
* **`input_file`**
  A txt file or CSV with a list of species names.
* **`input_database`**
  Path to the generated reference template FASTA file.

**NCBI Credentials (Required unless repeating):**
* **`-e, --email`** *STRING*
* **`-a, --api_key`** *STRING*

**Optional arguments:**
* **`-c, --maxcount`** *INT*
  Maximum accession numbers requested per species (default: 10000).
* **`-b, --batch_size`** *INT*
  Sequence download chunk size (default: 5000). Automatically halves if NCBI drops the connection.
* **`-l, --maxlength`** *INT*
  Maximum allowed sequence length (default: 22000).
* **`-z, --ampliconsize`** *INT*
  Minimum size an extracted amplicon may be (default: 50).
* **`-E, --evalue`** *INT*
  E-value exponent for BLAST filtering (default: 20, corresponding to 5e-20).
* **`-m, --mitochondria`**
  Appends `AND mitochondrion[filter]` to the NCBI query.
* **`-r, --ribosomal`**
  Appends `AND 12S` to the NCBI query.
* **`-q, --query`** *STRING*
  Custom NCBI search term addition.
* **`-s, --sort`**
  Sort downloads by sequence length (not recommended).
* **`-t, --taxid`**
  Use the previously saved TaxID list from the working directory.
* **`-R, --repeat`**
  Repeat curation on previously downloaded sequences without contacting NCBI.
* **`-T, --threads`** *INT*
  Parallel worker threads (auto-detected, capped at 7).

---

## 3. curate
Performs sequence quality control, length boundary filtering, MAFFT alignment, and phylogenetic tree construction (FastTree) for manual taxonomic inspection.

**Required arguments:**
* **`input_file`**
  Path to the uncurated BLAST output database (e.g., `BLAST_results/..._to_curate.fasta`).

**Optional arguments:**
* **`--min_length`** *INT*
  Minimum sequence length to retain. Shorter fragments are dropped prior to alignment (default: 150).
* **`--max_length`** *INT*
  Maximum sequence length to retain. Longer sequences are dropped prior to alignment.
* **`-N, --number_ns`** *INT*
  Maximum allowable ambiguous nucleotides (`N`s) per sequence (default: 0).
* **`-o, --old_database`** *FILE*
  Path to a previous version of the curated database.
* **`-M, --mafft_online`** *FILE*
  Path to a pre-aligned file if MAFFT was run externally/online.

---

## 4. complete
Fills missing taxonomic hierarchy levels, merges new sequences with older database versions, and finalizes the official FASTA file.

**Arguments:**
* **`-b, --blast_file`** *FILE*
  The new database BLAST file (`..._to_curate.fasta`).
* **`-c, --curated_file`** *FILE*
  The manually inspected, aligned FASTA file.
* **`-o, --old_database`** *FILE*
  Path to the existing master database to update/merge.
* **`-u, --updated_database`** *STRING*
  Filename for the finalized output database (default: `database.fasta`).

---

## 5. evaluate
Generates summary coverage dataframes, monophyletic/polyphyletic group diagnostics, and optional primer binding mismatch visual checks.

**Required arguments:**
* **`reference_database`**
  Path to the finalized reference FASTA database.
* **`monophyletic_group`**
  Path to the `..._post_curation_monophyletic_group.txt` file generated during the `complete` stage.

**Optional arguments:**
* **`-f, --forward_primer`** *STRING*
  Forward primer sequence to evaluate binding mismatches (5'-3').
* **`-r, --reverse_primer`** *STRING*
  Reverse primer sequence to evaluate binding mismatches (5'-3').

---

## 6. reformat
Reformats the finalized database headers and generates companion taxonomy files required by third-party classifiers.

**Required arguments:**
* **`reference_database`**
  Path to the finalized reference database.
* **`format`**
  Specifies the target classifier format. Must be one of:
  * `sintax` - SINTAX format (`>ACC;tax=d:...,p:...`)
  * `rdp` - Ribosomal Database Project format
  * `dadt` - DADA2 `assignTaxonomy` header format
  * `dads` - DADA2 `assignSpecies` header format
  * `idt` - IDTAXA format (`>Root;Domain;...`)
  * `qiime` - QIIME 2 FASTA and companion 2-column taxonomy `.txt` file

# TROUBLESHOOTING & NCBI RATE LIMITS
When processing high-volume taxa, NCBI Entrez may drop HTTP connections, resulting in `IncompleteRead(0 bytes read)` errors.

**Adaptive Chunking:**
EchoPipe automatically handles this by retrying the chunk once. If it fails again, it dynamically halves the batch size (e.g., 5000 → 2500 → 1250) to allow large species to pass without crashing the pipeline.

**Manual Intervention:**
If connection errors persist on your network, manually lower the thread count and batch size:
`echopipe create species.txt ref.fasta -T 2 -b 1000`

# EXAMPLES

**1. Create a reference template using a random subset of 50 species:**
`echopipe template species.txt -f ACACCGCCCG -r GTAYACTTACC -n 50 -e user@mail.com -a API_KEY`

**2. Mine target marker sequences with strict E-value filtering:**
`echopipe create species.txt ref_template.fasta -e user@mail.com -a API_KEY -E 20 -b 5000`

**3. Curate sequences by applying strict length boundaries:**
`echopipe curate BLAST_results/run_to_curate.fasta --min_length 200 --max_length 600 -N 0`

**4. Reformat the final database for QIIME 2:**
`echopipe reformat Database_2026.fasta qiime`

# AUTHOR
EchoPipe Development Team.

# SEE ALSO
mafft(1), fasttree(1), blastn(1)