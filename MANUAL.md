# EchoPipe Manual

EchoPipe is a complete pipeline for reference database creation and curation. The script `echopipe.py` is divided into six main subcommands: `template`, `create`, `curate`, `complete`, `evaluate`, and `reformat`.

---
## 1. `template`
Generate a template reference database.

**Positional Arguments:**
* **`input_file`**: Txt or CSV file species names or a fasta file.
* **`input_file_species`**: Txt or CSV file species names (used specifically when completing the template with the `-C` flag).

**Optional Arguments:**
* **`-f, --forward`**: The forward primer used to find region of interest, (5'-3').
* **`-r, --reverse`**: The reverse primer used to find region of interest, (5'-3').
* **`-e, --email`**: Your email if NCBI needs to contact you (default: `email@email.com`).
* **`-a, --api_key`**: The user's NCBI API key (default: `api_key`).
* **`-q, --query`**: Custom query additions.
* **`-t, --threshold`**: The minimum length of a sequence (default: `150`).
* **`-l, --length`**: The longest allowed sequence length (default: `22000`).
* **`-m, --max`**: Number of sequences downloaded per species (default: `1`).
* **`-p, --provided_sequences`**: Use a fasta file as reference template.
* **`-z, --longest_amplicon_size`**: Multiplier for median length (default: `2.0`).
* **`-n, --random_subset`**: Number of random species to use.
* **`-sf, --subset_file`**: Path to a file containing a specific subset of species, typically a taxonomically diverse group of the target taxa.
* **`-T, --threads`**: Number of parallel threads to use (default: auto-detected, max 7).
* **`--rate_limit`**: Max requests per second for NCBI API (default: 7.0 with API key, 2.5 without).
* **`-C, --Complete`**: Completes the reference template database.

### Example Usage:
**Initial run:**
```
bash
python echopipe.py template species_list.csv -f GTCGGTAAAACTCGTGCCAGC -r CATAGTGGGGTATCTAATCCCAGTTTG -e email@email.com -a your_api_key
```

**Initial run, completion of template after manual curation:**
```  
bash 
python echopipe.py template -C unique_species_list.csv
```
---
## 2. `create`
Mine NCBI for reference sequences and creates a BLAST-ready database.

**Positional Arguments:**
* **`input_file`**: A txt file or CSV with a list of species names.
* **`input_database`**: Path to the input reference database fasta file.

**Optional Arguments:**
* **`-e, --email`**: User's email address (default: `email@email.com`).
* **`-a, --api_key`**: User's NCBI API key (default: `api_key`).
* **`-s, --sort`**: Sort by length (Not recommended).
* **`-c, --maxcount`**: Maximum accession numbers per species (default: `10000`).
* **`-l, --maxlength`**: Longest allowed sequence length (default: `22000`).
* **`-z, --ampliconsize`**: Minimum size an amplicon may be (default: `50`).
* **`-m, --mitochondria`**: Search targets mitochondrial sequences.
* **`-r, --ribosomal`**: Search for mitochondrial 12S ribosomal DNA.
* **`-q, --query`**: Custom NCBI search term.
* **`-b, --batch_size`**: Batch size for downloading sequences (default: `5000`).
* **`-t, --taxid`**: Use last saved taxid list.
* **`-E, --evalue`**: E-value for BLAST (default: `20`). Increase for longer markers; keep lower for shorter markers to avoid introducing non-target gene regions.
* **`-R, --repeat`**: Repeat curation on previously downloaded sequences.
* **`-T, --threads`**: Number of parallel threads to use (default: auto-detected, max 7).
* **`--rate_limit`**: Max requests per second for NCBI API (default: 7.0 with API key, 2.5 without).

**Initial run, create the raw reference database:**

```
bash 
python echopipe.py create species_list.csv reference_template_database.fasta
```
---
## 3. `curate`
Align and generate trees for manual curation.

**Positional Arguments:**
* **`input_file`**: Database to revise.

**Optional Arguments:**
* **`-o, --old_database`**: The previous version of database.
* **`-N, --number_ns`**: Number of N's and ambiguous nucleotides allowed (default: `0`).
* **`-M, --mafft_online`**: Path to MAFFT online alignment file.
* **`--min_length`**: Minimum sequence length to keep (default: `150`).
* **`--max_length`**: Maximum sequence length to keep (default: infinity).
* **`-f, --forward-primer`**: Forward primer sequence.
* **`-r, --reverse-primer`**: Reverse primer sequence.

```
bash 
python echopipe.py curate BLAST_results/{date}_{run_number}_to_curate.fasta --min_length 150 and --max_length 250
```
---
## 4. `complete`
Filter, merge, and finalize the database.

**Optional Arguments:**
* **`-b, --blast_file`**: New database BLAST file.
* **`-c, --curated_file`**: Curated aligned FASTA.
* **`-o, --old_database`**: Existing master database.
* **`-u, --updated_database`**: Name of your new database, (default: `database_{date}_{run_number}.fasta`).

```
bash
python echopipe.py complete -b BLAST_results/{date}_{run_number}_to_curate.fasta -c Database_curation/{date}_{run_number}/{date}_{run_number}_aligned.fasta -u Database_name_{date}_{run_number}.fasta
```
---
## 5. `evaluate`
Evaluate the database (GC content, primers).

**Positional Arguments:**
* **`reference_database`**: Path to the reference database.
* **`monophyletic_group`**: Path to the monophyletic groups text file.

**Optional Arguments:**
* **`-f, --forward_primer`**: The forward primer sequence to check (5'-3').
* **`-r, --reverse_primer`**: The reverse primer sequence to check (5'-3').
```
bash
python echopipe.py evaluate Database_name_{date}_{run_number}.fasta Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_post_curation_monophyletic_group.txt
```
---

## 6. `reformat`
Reformat database headers for popular taxonomic classifiers.

**Positional Arguments:**
* **`reference_database`**: The reference database for which the header format will be changed.
* **`format`**: Write in one of the available formats: 
  * `sintax` = SINTAX
  * `rdp` = Ribosomal Database Project (RDP)
  * `dadt` = DADA2 assignTaxonomy
  * `dads` = DADA2 assignSpecies
  * `idt` = IDTAXA
  * `qiime` = QIIME 2
 
## 7. `updating the database`
To update an existing reference database with newly available sequences from NCBI or to expand the sequence coverage, run the create command again.

If a species exceeded the previous download limit, you can run create with an increased maximum count (-c), add a specific query filter (e.g., --query "12s"), or run it at a later time when new accession numbers are published. New accession numbers will be automatically detected and processed against your master database.

Example Usage (Updating with a higher count, a custom query and the old reference database as template)
```
bash
python echopipe.py create unique_species_list.csv Database_name_{date}_{run_number}.fasta -c 20000 --query "12s" -e email@email.com -a your_api_key
```
