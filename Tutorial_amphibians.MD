# Tutorial on how to use the code base EchoPipe
## Test run on Scandinavian amphibians and batra primer
## Background
EchoPipe is an iterative, reproducible pipeline for creating, curating and evaluating reference databases for environmental DNA metabarcoding studies.  
By being an interative workflow, the workload, memory- and time usage is reduced, while coverage of the species of interest is increased.  
The entire pipeline is written in Python, and is distributed through Conda, which eases and streamlines the installation and usage. 
For installation, see README.md

## Requirements
### There are few required inputs for the workflow:  
** [NCBI API key](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) **  
** Species list with scientific name **  
** Email **  
** Forward and reverse primer sequence spanning the region of interest (both in 5'-3') **  
** [Conda version 24.5.0](https://anaconda.org/anaconda/conda/files?version=24.5.0) **  

Use the following code to ensure you have the correct conda version.

```
conda --version
Expected output: conda 24.5.0 
``` 
If the results is not conda 24.5.0, please update your conda version.

![image](https://github.com/user-attachments/assets/28efafb2-4204-4685-9075-5410de525229)
Flowchart of the EchoPipe's workflow. Each color represents a module of the workflow and the corresponding script. Light-yellow: EchoPipe – Reference template database creator. Yellow: EchoPipe – Database creation. Orange: EchoPipe – Database curation. Red: EchoPipe – Database completion. Illustration is made with InkScape.

## Preparation
This tutorial assumes that the EchoPipe scripts is located in the RAWDIR, in this example the Final_version path. After each program is finished, a recommended code line is presented to continue the workflow.

```
RAWDIR=Final_version
cd $RAWDIR
```

Download the conda environment for WSL, create it and activate it, using the following code.

```
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/refs/heads/main/environment_240930.yml
conda env create -f environment_240930.yml
conda activate EchoPipe 

```

## Create an initial template reference.
### Options for Echopipe_reference_template_creator.py  
options:
  -h, --help            show this help message and exit
  -n RANDOM_SUBSET, --random_subset RANDOM_SUBSET
                        Number of random species to use for template creation (e.g., 50).
  -sf SUBSET_FILE, --subset_file SUBSET_FILE
                        Path to a file containing a specific subset of species to use.

Arguments used to generate an uncurated reference template.
The required arguments are input_file, -f (--forward), -r (--reverse), -e (--email) and -a (--api_key):
  input_file            Txt or CSV file species names or a fasta file that is to be converted into the reference template database (-p is then required).
  -f FORWARD, --forward FORWARD
                        The forward primer used to find region of interest, (5'-3' orientation).
  -r REVERSE, --reverse REVERSE
                        The reverse primer used to find region of interest, (5'-3' orientation).
  -e EMAIL, --email EMAIL
                        Your email if NCBI needs to contact you.
  -a API_KEY, --api_key API_KEY
                        The user's NCBI API key, allows for faster downloads.
  -q QUERY, --query QUERY
                        The search result will include the user input search term(s). Example, limit the search to 12s region: -q "AND 12s" followed by the search term. To exclude a term write "NOT 12s".
  -t THRESHOLD, --threshold THRESHOLD
                        The minimum length of a sequence, including the primer regions. Any sequence shorter than this is discarded. Default cutoff is set to 150 bases.
  -l LENGTH, --length LENGTH
                        The longest allowed sequence length for template creation Default is set to 20 000. WARNING: The longer the sequence the more computational power is required to align the sequences.
  -m MAX, --max MAX     The number of sequences that are downloaded per species. Increasing the number may increase coverage while increasing computational power. Recommended total should not exceed 500.
                        Default = 1.
  -p, --provided_sequences
                        Use if a fasta file is provided to be used as a reference template.
  -z LONGEST_AMPLICON_SIZE, --longest_amplicon_size LONGEST_AMPLICON_SIZE
                        Advanced setting: multiplier for the median length of downloaded sequences to discard abnormally long sequences after trimming. Default is set to 2*median length.

Argument used to finish the curated reference template database.:
  -C, --Complete        Completes the reference template database.
  input_file_species    Txt or CSV file species names or a fasta file that is to be converted into the reference template database (-p is then required).

```
python Echopipe_reference_template.py Norwegian_swedish_amphibians_250214.csv -f ACACCGCCCGTCACCCT -r GTAYACTTACCATGTTACGACTT -l 1000 -e your.email -a your.api.key -t 40 -m 5 -q "AND 12S"  

```
Inspect aligned_sequences_to_curate.fasta with MSA-tool (Ex. Jalview).  
Remove clear dubious or wrongly annotated sequences from the aligned_sequences_to_curate.fasta file. 
No sequences should expand outside the primers. 
The most important is that sequences from only the correct marker region is included.  
Save the file with the same name, aligned_sequences_to_curate.fasta.

### Output files:
Log_files/
{Run_number}_log.txt, a log file with output of pipeline.  
  
Reference_template_creation/
aligned_sequences_to_curate.fasta, contains aligned versions of the reference sequences, use this for curation of reference template database.
filtered_aligned_sequences.fasta, contains aligned versions of the reference sequences (gaps removed).
non_approved_sequenes.txt, contains headers of sequences not accepted by the search terms.  
preformated_sequences.fasta, contains all non-aligned sequences analyzed.  
temp_split_file.fasta, contains batch-aligned sequences (temp files).  

duplicate_taxid_entires.txt, contains dupliacte species entries.  
unique_{species_list}.csv, contains a list of unique and validates species names.

### Example:
Open the aligned_sequences_to_curate.fasta, which is aligned using MAFFT (linsi) in Jalview, choose Colour by nucleotide to visualize sequence similarity.  
Scroll until you find a conserved region which overlaps with both primers, this ensure that the sequences in the reference template database is of the correct marker region.  
Keep in mind that gaps may be introduced as the alignment may contain unwanted sequences from other marker regions.  

<img width="952" alt="image" src="https://github.com/user-attachments/assets/8b8a3a89-e216-4122-8ac2-b3f48cf793f9">

Due to the reorder algorithm in the MSA, similiar sequences are stored together, and sequences from other marker regions are normally found at the bottom or top of the MSA.  
Highlighting and deleting (press delete) all sequences which are not similar to the primers in the file, is an effective way to remove sequences from other marker regions.  
  
Then save the file, and do not alter the path nor the file name.  





### Example:
Open the aligned_sequences_to_curate.fasta, which is aligned using MAFFT (linsi) in Jalview, choose Colour by nucleotide to visualize sequence similarity.  
Scroll until you find a conserved region which overlaps with both primers, this ensure that the sequences in the reference template database is of the correct marker region.  
Keep in mind that gaps may be introduced as the alignment also contain potential reference sequences from other marker regions.  

### Complete the creation and curation of the reference template database.
A file called reference_template_database.fasta, is then created in the $RAWDIR.  
The temp files is removed.  
The log file is updated.  

```
cd $RAWDIR
python Echopipe_reference_template.py unique_Norwegian_swedish_amphibians_250214.csv -e your.email -a your.api.key -C
```

## Create a reference database  
### Use the script Echopipe_database_creation.py  
This script is iterative, and could be re-run by changing the search options.
By dereplicating 100% identical sequences within a species, a counter is introduced, which reduce the file size, and help with curation of the database.
The format of the database is as follows:
```
>gb|Accession_number|NCBI_taxonomy|counter
Fasta_sequence

Example:

>gb|KU244636.1|Eukaryota;Chlorophyta;Mamiellophyceae;Mamiellales;Mamiellaceae;Micromonas;Micromonas_bravo|5
CCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAG
CTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTT
TCAACTTCCTGTAGAGGACGCGCTCTGGCTTCACGGCTGGACGCGGAGTCTACGTGGTTA
CTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAA
CACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAA
CAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGA
ACTTCTGCGAAAGCATTTGCCAAGGATGTTTTCATT
```
### Options for Echopipe_reference_template_database_creator.py  

options:
  -h, --help            show this help message and exit

Arguments used to generate the uncurated reference database.
The required arguments are input_file, input_database, -e (--email) and -a (--api_key):
  input_file            A txt file or CSV with a list of species names.
  input_database        Path to the input reference database fasta file.
  -e EMAIL, --email EMAIL
                        User's email address associated with the NCBI API key.
  -a API_KEY, --api_key API_KEY
                        User's NCBI API key.
  -s, --sort            Sort by length, if implemented, targets longer sequences (Not recommended).
  -c MAXCOUNT, --maxcount MAXCOUNT
                        Maximum number of accession numbers downloaded per species, default = 10 000.
  -l MAXLENGTH, --maxlength MAXLENGTH
                        Longest allowed sequence length for an accession number to be analysed, default = 22 000.
  -z AMPLICONSIZE, --ampliconsize AMPLICONSIZE
                        Minimum size an amplicon may be in order to accepted. Consider adjusting this based on the marker region.
                        Default = 50.
  -m, --mitochondria    Search targets mitochondrial sequences.
  -r, --ribosomal       Search result will only return accession numbers that have been annotated with mitochondrial 12S ribosomal
                        DNA.
  -q QUERY, --query QUERY
                        Search result will include the user input search term(s). Example, limit the search to 12S region: -q "AND
                        12S". To exclude a term write "NOT 12S".
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Amount of sequences that can be downloaded simultaneously. Default = 5000. Only need to consider using this
                        with values lower than the default and only use it in conjunction with --maxlength and expecting to
                        download very large sequences (size of chromosomes).
  -t, --taxid           Use last saved taxid list. Can be used to change NCBI taxonomy (Not recommended)
  -E EVALUE, --evalue EVALUE
                        The E-value used for BLAST. The higher the value entered, the more stringent the BLAST becomes. Default =
                        12, meaning: 10^-12.
  -R, --repeat          No new sequences are downloaded. The curation step needs to be ran first. Instead, all previously
                        downloaded sequences are BLASTed against the input database. Note: this options resets the file containing
                        analysed accession numbers and the output will be curated from the start.AATCAAGAACGAAAGT
```

### If this is the initial round:

```
python Echopipe_database_creation.py unique_Norwegian_swedish_amphibians_250214.csv reference_template_database.fasta -e your.email -a your.api.key
```
The scripts gives live updates, and tells the users what species are being processed.  
The scripts first sees if there are any new un-analyzed accession numbers on NCBI, and downloads if there are any.  
If there where any new un-analyzed accession numbers, the scripts then starts to download those.  
To check if the sequences related to the accession numbers are from the wanted marker region, a BLAST alignment against the initial reference template database is conducted.  
The longest match which are higher than the set E-value is kept, and are found in the $RAWDIR/BLAST_results/{date}_{run_number}_to_curate.fasta  

## Curate database
To ensure that the sequences are correctly annotated, Echopipe_database_curation.py, generates several tools to detect false-annotations.  
Example:
```
python Echopipe_database_curation.py BLAST_results/2024-09-30_1_to_curate.fasta
