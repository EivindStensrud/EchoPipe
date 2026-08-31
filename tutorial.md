# Tutorial: How to Use EchoPipe
## Test run on Scandinavian amphibians and batra primers

## Background
EchoPipe is an iterative, reproducible pipeline for creating, curating, evaluating, and reformatting reference databases for environmental DNA (eDNA) metabarcoding studies. The tool allows for creating reference databases for every given taxa, primer set, and genetic region. 

By utilizing an iterative workflow, computational workload, memory, and runtime are minimized while maximizing target taxonomic coverage. The entire pipeline is distributed as a Python package installable via `pip` with a unified command-line interface (`echopipe`).

---

## Installation & Requirements

### Mandatory Inputs & Dependencies
* **[NCBI API Key](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317)** *(Increases rate limits up to 10 requests/sec)*
* **Species list with scientific names** (CSV or TXT)
* **Email address** (for NCBI Entrez contact)
* **Forward and Reverse primer sequences** spanning the target region (5'-3' orientation)
* **Python ≥ 3.9**

### Step-by-Step Installation

1. **Create and activate a virtual environment** (using `conda`):
   ```bash
   conda create -n echopipe python=3.10 -y
   conda activate echopipe

Install EchoPipe via pip:

Bash
pip install git+[https://github.com/EivindStensrud/EchoPipe.git](https://github.com/EivindStensrud/EchoPipe.git)

Verify Installation:

Bash
echopipe --help

<img width="982" height="1024" alt="Figure_1_workflow" src="https://github.com/user-attachments/assets/d2088ad1-da3c-4aa4-b14c-b1fa9c9e3785" />

Figure 1. EchoPipe modular workflow for iterative database creation and curation. The pipeline architecture is divided into primary stages: Reference Template Generation (yellow), Sequence Retrieval and Database Creation (orange), Diagnostic Curation (dark orange), Database Completion/Evaluation (red), and Classifier Reformatting.

Stage 1: Generate Reference Template (template)
The template subcommand generates an initial, uncurated marker reference database using primer alignment matching.

Create an initial template reference (template)  
The template subcommand generates a template reference database.  
Full Options for echopipe templateUncurated template generation arguments:  
input_file: Txt or CSV file species names or a fasta file.  
-f, --forward: The forward primer used to find region of interest, (5'-3').  
-r, --reverse: The reverse primer used to find region of interest, (5'-3').  
-e, --email: Your email if NCBI needs to contact you.  
-a, --api_key: The user's NCBI API key.  
-q, --query: Custom query additions.  
-t, --threshold: The minimum length of a sequence (default: 150).  
-l, --length: The longest allowed sequence length (default: 22000).  
-m, --max: Number of sequences downloaded per species (default: 1).  
-p, --provided_sequences: Use a fasta file as reference template.  
-z, --longest_amplicon_size: Multiplier for median length (default: 2).  
-n, --random_subset: Number of random species to use.  
-sf, --subset_file: Path to a file containing a specific subset of species.  
-T, --threads: Number of parallel threads to use (default: auto-detected, max 7).  

Arguments used to finish the curated reference template database:  
-C, --Complete: Completes the reference template database.   
input_file_species: Txt or CSV file species names.  

1. Build Initial Unaligned Template
Bash
echopipe template Norwegian_swedish_amphibians_250214.csv \
  -f ACACCGCCCGTCACCCT \
  -r GTAYACTTACCATGTTACGACTT \
  -l 1000 \
  -e your.email \
  -a your.api.key \
  -t 40 \
  -m 5 \
  -q "AND 12S"  


2. Inspect and Clean Template Alignment
Open Reference_template_creation/aligned_sequences_to_curate.fasta in an MSA viewer (e.g., Jalview).

Color alignment by nucleotide to view conserved regions.

Verify that sequences overlap both forward and reverse primer binding sites.

Remove obviously dubious sequences or non-target genomic fragments.

Save the file without modifying its filename or relative directory path.

3. Complete Reference Template
Run the completion flag (-C) to build reference_template_database.fasta:  

Bash
echopipe template unique_Norwegian_swedish_amphibians_250214.csv \
  -e your.email \
  -a your.api.key \
  -C  


Stage 2: Create Reference Database (create)
The create subcommand mines NCBI for candidate sequences, extracts local marker regions via local BLAST against the template database, and handles network disruptions automatically.

Key CLI Options for echopipe create
input_file: A txt file or CSV with a list of species names.  
input_database: Path to the input reference database fasta file.  
-e, --email: User's email address.  
-a, --api_key: User's NCBI API key.  
-s, --sort: Sort by length (Not recommended).  
-c, --maxcount: Maximum accession numbers per species (default: 10000).  
-l, --maxlength: Longest allowed sequence length (default: 22000).  
-z, --ampliconsize: Minimum size an amplicon may be (default: 50).  
-m, --mitochondria: Search targets mitochondrial sequences.  
-r, --ribosomal: Search for mitochondrial 12S ribosomal DNA.  
-q, --query: Custom NCBI search term.  
-b, --batch_size: Batch size for downloading sequences (default: 5000).  
-t, --taxid: Use last saved taxid list.  
-E, --evalue: E-value for BLAST  (default: 20, indicating 5e-20). Increase for longer markers; keep lower for shorter markers to avoid introducing non-target gene regions.  
-R, --repeat: Repeat curation on previously downloaded sequences.  
-T, --threads: Number of parallel threads to use (default: auto-detected, max 7).  

Tip (E-Value Selection): Keep the E-value lower (less strict, ex. 15) for shorter markers to avoid introducing non-target gene regions. Increase the E-value (stricter, ex. 40) for longer markers to avoid introducing non-target gene regions.

Standard Database Creation Command
Bash
echopipe create unique_Norwegian_swedish_amphibians_250214.csv reference_template_database.fasta \
  -e your.email \
  -a your.api.key \
  -E 20 \
  -b 5000 \
  -T 4


Output Database Header Format
Dereplicated sequences track abundance using an appended sequence counter:

Plaintext
>gb|Accession_number|NCBI_taxonomy|counter
Sequence_data

Example:
>gb|KJ858774.1|Eukaryota;Chordata;Amphibia;Anura;Alytidae;Alytes;Alytes_obstetricans|2
ACACCGCCCGTCACCCTCCTCAACTAACTCAACCCCCTAACTAAAAGCTAACTGGTTAACAAGAAGAGGCAAGTCGTAACATGGTAAGTATA
Stage 3: Curate Database (curate)
The curate subcommand aligns extracted marker regions, applies strict sequence length boundary controls, builds FastTree gene trees, and flags potential misannotations or paralogs.

CLI Options for echopipe curate
input_file: Database to revise.  
-o, --old_database: The previous version of database.  
-N, --number_ns: Number of N's and ambiguous nucleotides allowed (default: 0).  
-M, --mafft_online: Path to MAFFT online alignment file.  
--min_length: Minimum sequence length to keep (default: 150).  
--max_length: Maximum sequence length to keep.  

Execution Command with Target Bounds
Bash
echopipe curate BLAST_results/2026-03-17_1_to_curate.fasta \
  --min_length 200 \
  --max_length 600 \
  -N 0

Diagnostic Curation Output Files
Database_curation/{date}/{date}_{run_number}_aligned.fasta: Aligned dataset (* marking new sequences).

Database_curation/{date}/{date}_{run_number}_tree_string.newick: Gene tree for visualization in iTOL or TreeViewer.

Database_curation/{date}/{date}_{run_number}_duplicate_sequences.txt: Identical sequences shared across different species.

Database_curation/{date}/{date}_{run_number}_sequences_to_delete.txt: User text file where misannotated accession numbers are placed to exclude them from the final database.

CLI Options for echopipe complete
-b, --blast_file: New database BLAST file.  
-c, --curated_file: Curated aligned FASTA.  
-o, --old_database: Existing master database.  
-u, --updated_database: Output database (default: database.fasta).  

Stage 4: Complete Official Database (complete)
Once misannotated accession numbers are added to sequences_to_delete.txt, finalize your official reference database.

Initial Database Finalization
Bash
echopipe complete \
  -b BLAST_results/2026-03-17_1_to_curate.fasta \
  -c Database_curation/2026-03-17_1/2026-03-17_1_aligned.fasta \
  -u Amphibian_2026-03-17_1.fasta


Updating an Existing Database Version
Merging new NCBI data into a previous EchoPipe database version preserves historical counters:

Bash
echopipe complete \
  -b BLAST_results/2026-03-17_2_to_curate.fasta \
  -u Amphibian_2026-03-17_2.fasta \
  -o Amphibian_2026-03-17_1.fasta

Generated Visualizations
Stage 5: Evaluate Database (evaluate)
The evaluate subcommand tests database completeness, taxid resolution, GC-content distribution, monophyly, and primer binding compatibility.

CLI Options for echopipe evaluate
reference_database: Path to the reference database.  
monophyletic_group: Path to the monophyletic groups text file.  
-f, --forward_primer: The forward primer sequence to check (5'-3').  
-r, --reverse_primer: The reverse primer sequence to check (5'-3').  

Bash
echopipe evaluate Amphibian_2026-03-17_1.fasta \
  Database_curation/2026-03-17_1/Curated_content/2026-03-17_1_post_curation_monophyletic_group.txt \
  -f ACACCGCCCGTCACCCT \
  -r GTAYACTTACCATGTTACGACTT
Primer Mismatch & Proportions Output
Upon completion, evaluate provides a soft suggestion to reformat the database for downstream taxonomic assignment classifiers.

Stage 6: Reformat for Taxonomic Classifiers (reformat)
The reformat subcommand converts official EchoPipe reference databases into header formats compatible with major bioinformatics tools.

Supported Formats
Format Code	Target Classifier	Example Header Output
sintax	SINTAX	>ACC;tax=d:Eukaryota,p:Chordata,c:Amphibia...
rdp	RDP Classifier	>ACC\troot;Eukaryota;Chordata;Amphibia...
dadt	DADA2 assignTaxonomy	>Eukaryota;Chordata;Amphibia;Anura;Alytidae;Alytes
dads	DADA2 assignSpecies	>ACC Alytes obstetricans
idt	IDTAXA	>Root;Eukaryota;Chordata;Amphibia;Anura;Alytidae...
qiime	QIIME 2	>ACC (+ companion 2-column taxonomy .txt file)
Example Reformat Commands
Bash
# Reformat for QIIME 2 (generates .fasta and .txt taxonomy mapping file)
echopipe reformat Amphibian_2026-03-17_1.fasta qiime

# Reformat for DADA2 assignTaxonomy
echopipe reformat Amphibian_2026-03-17_1.fasta dadt

# Reformat for SINTAX
echopipe reformat Amphibian_2026-03-17_1.fasta sintax


# Troubleshooting & Performance Tips
Handling NCBI Connection Issues (IncompleteRead(0 bytes read))
High-volume species (e.g., model organisms like Aquarana catesbeiana) can cause NCBI's Entrez server to drop HTTP streams.

EchoPipe includes adaptive download logic:

Initial Run: Attempts download using your specified -b / --batch_size (default: 5000).

First Retry: Re-attempts download with the same batch size after a 10-second pause.

Subsequent Retries: Dynamically halves the batch size (2500 → 1250 → 250) to bypass server timeouts.

If persistent network failures occur, manual overrides can be applied:

Bash
# Lower worker threads (-T) and batch size (-b) for unstable connections
echopipe create species.txt ref_template.fasta -e user@mail.com -a API_KEY -T 2 -b 1000