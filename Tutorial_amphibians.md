# Tutorial on how to use the code base EchoPipe
## Test run on Scandinavian amphibians and batra primer
## Background
EchoPipe is an iterative, reproducible pipeline for creating, curating and evaluating reference databases for environmental DNA metabarcoding studies.  
By being an interative workflow, the workload, memory- and time usage is reduced, while coverage of the species of interest is increased.  
The entire pipeline is written in Python, and is distributed through Conda, which eases and streamlines the installation and usage. 
For installation, see README.md

## Requirements
### There are few required inputs for the workflow:  
**[NCBI API key](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317)**    
**Species list with scientific name**  
**Email**  
**Forward and reverse primer sequence spanning the region of interest (both in 5'-3')**  
**[Conda version 24.5.0](https://anaconda.org/anaconda/conda/files?version=24.5.0)**  

Use the following code to ensure you have the correct conda version.

```
conda --version
Expected output: conda 24.5.0 
``` 
If the results is not conda 24.5.0, please update your conda version.

<img width="982" height="1024" alt="Figure_1_workflow" src="https://github.com/user-attachments/assets/d2088ad1-da3c-4aa4-b14c-b1fa9c9e3785" />
  
Flowchart of the EchoPipe's workflow. Each color represents a module of the workflow and the corresponding script. Figure 1. EchoPipe modular workflow for iterative database creation and curation. The pipeline architecture is divided into four primary stages: Reference Template Generation (light yellow), Sequence Retrieval and Database Creation (light orange), Diagnostic Curation (orange), and Database Completion/Evaluation (red). The iterative feedback loop (bottom center) enables users to continuously update curated databases with novel accessions from public repositories without duplicating previous computational efforts.

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

Argument used to finish the curated reference template:  
  -C, --Complete        Completes the reference template.  
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

<img width="526" height="398" alt="Reference_template" src="https://github.com/user-attachments/assets/e8fa1a18-c8b5-4fa0-a81d-0ede996f87c6" />    
  
Visualization of the MSA alignment of the curated reference template file using Scandinavian amphibian species list and batra primers (Valentini et al. 2016)  

Due to the reorder algorithm in the MSA, similiar sequences are stored together, and sequences from other marker regions are normally found at the bottom or top of the MSA.  
Highlighting and deleting (press delete) all sequences which are not similar to the primers in the file, is an effective way to remove sequences from other marker regions.  
  
Then save the file, and do not alter the path nor the file name.  


### Complete the creation and curation of the reference template.  
A file called reference_template_database.fasta, is then created in the $RAWDIR.  
The temp files is removed.  
The log file is updated.  

```
cd $RAWDIR
python Echopipe_reference_template.py unique_Norwegian_swedish_amphibians_250214.csv -e your.email -a your.api.key -C
```

## Create a reference database  
### Use the script Echopipe_database_creation.py  
### Options for Echopipe_database_creation.py   
-a, Api.key from NCBI (required).  
-e, Email address for NCBI to contact you (required).  
-s, Sort, if included, longer sequences are targeted (optional).  
-c, Max count, The maximum number of accession numbers downloaded per species. Default = 10 000.  
-q, Query search word, target different marker regions using NCBI search term(s), ex. -q "AND 18s".  
-l, Maximum length of amplicon retrived from NCBI, default=22 000 (optional).  
-t, Taxid, use this if last run of taxid wants to be used. Only works if this is not the first run (optional).  
-b, Batch size, The amount of sequences that can be downloaded simultaneously. Default = 5000 (optional).  
-E, E-value, The E-value used for BLAST. The higher the value entered, the more stringent the BLAST becomes. Default = 5, (10^-5) (optional).  
-R, Repeat, No new sequences are downloaded. Instead, all previously downloaded sequences are BLASTed against the new input database. Note: this options resets the file containing analysed accession numbers and the output will be curated from the start (optional).  
-h, Display help.  


This script is iterative, and could be re-run by changing the search options.  
By dereplicating 100% identical sequences within a species, a counter is introduced, which reduce the file size, and help with curation of the database.  
The format of the database is as follows:  
```
>gb|Accession_number|NCBI_taxonomy|counter  
Fasta_sequence  

Example: 
gb|KJ858774.1|Eukaryota;Chordata;Amphibia;Anura;Alytidae;Alytes;Alytes_obstetricans|2  
ACACCGCCCGTCACCCTCCTCAACTAACTCAACCCCCTAACTAAAAGCTAACTGGTTAACAAGAAGAGGCAA
GTCGTAACATGGTAAGTATA
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

### If older version exist: 
We recommend to use the older version of the completed reference database as the input file, this may increase the retrival rate of potential reference sequences. 
  
Explore the {date}_over_max.csv file to see if some species have more unanalyzed sequences.  
If you want to increase from the default 10 000 sequences, you can use the option, -c 20 000, to increase to 20 000 analyzed sequences.  
Other options are to include or exclude search words, the default is no search words, which retrieves all DNA entries on NCBI.  
Including search words will reduce the number of potential reference sequences analyzed, and therefore recommended to avoid.  
If there are way more entries than 10 000, you may want to limit the search by including the marker region of interest, by using -q "AND 12s", or -q "AND mitochondrial".   
Creating optimal search terms may be challenging, and species group specific, explore we recommend to explore with different terms.  

Example:
This increases the number of sequences downloaded per species from 10 000 to 20 000.
```
python Echopipe_database_creation.py unique_Norwegian_swedish_amphibians_250214.csv Name.of.excisting.database.fasta -e your.email -a your.api.key -c 20000
```



There are several log files create when Echopipe_database_creation.py is ran.  
$RAWDIR/Log_files/{date}_{run_number}_log.txt, Is a log file storing the setttings and run output.  
$RAWDIR/Log_files/{date}_{run_number}_taxid_collection.txt, Stores NCBI species name, taxid, name from species list, and NCBI taxonomy.  
$RAWDIR/Log_files/duplicate_species_entries.txt, If two or more species name from the species input file share the same NCBI species name, they are stored here.  
$RAWDIR/Log_files/species_not_found.txt, Species in the species input file not found in the NCBI taxonomy are stored here.  

For each species, information about sequences obtained and analysed are stored per run are created.  
$RAWDIR/Species/{Species_name}/{Species_name}_{date}_{run_number}_accession_numbers.txt, stores analyzed accession numbers each run.  
$RAWDIR/Species/{Species_name}/{Species_name}_{date}_{run_number}_fasta_blast_results.txt, stores alignment data.  
$RAWDIR/Species/{Species_name}/{Species_name}_{date}_{run_number}_fasta_blast_results_check.fasta, stores potential new reference sequences.  
$RAWDIR/Species/{Species_name}/{Species_name}_analysed.accession_numbers.txt, stores all analyzed accession numbers from earlier runs.  
$RAWDIR/Species/{Species_name}/{Species_name}_new_accession_numbers_{date}_{run_number}.txt, stores new accession numbers each run, which have not been analyzed before.  
$RAWDIR/Species/{Species_name}/{Species_name}_new_seqs_{date}_{run_number}.fasta, stores the sequences from the new accession numbers, which have not been analyzed before.  

For each species, information about new accepted reference sequences are created in $RAWDIR/BLAST_results/.
$RAWDIR/BLAST_results/{date}_{run_number}_species_with_results.txt, stores what species accepted potential new reference sequences are found from.
$RAWDIR/BLAST_results/{date}_{run_number}_to_curate.fasta, stores what sequences which needs to be manually curated.

$RAWDIR/Template_databases/reference_template_database/, The BLAST compatible version of the reference database.

$RAWDIR/Database_curation/{date}/{date}_{run_number}_duplicate_species_entries.txt, If two or more species name from the species input file share the same NCBI species name, they are stored here.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_species_not_found.txt, If a species is not found, it is stored here.  

$RAWDIR/taxid_collection.txt, Stores NCBI species name, taxid, name from species list, and NCBI taxonomy.
$RAWDIR/{date}_{run_number}_over_max.csv, Stores information of what species have more max analyzed accession numbers (-c).  



## Curate database
To ensure that the sequences are correctly annotated, Echopipe_database_curation.py, generates several tools to detect false-annotations.  
Example:
```
python Echopipe_database_curation.py BLAST_results/2026-03-17_1_to_curate.fasta
```

The location and output of the curation tool.
$RAWDIR/Database_curation/{date}/{date}_{run_number}_aligned.fasta, Shows new sequences and old sequences aligned. New sequences are marked with *, note that first run all are marked.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_concatenated_file.fasta, Shows new sequences and old sequences not aligned. New sequences are marked with *, note that first run all are marked.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_dataframe.csv, Shows a brief summary if the species forms a monophyletic group.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_dubious_sequences.txt, Creates an empty files where the user may add weird or dubios sequences.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_sequences_to_delete.txt, An empty file where dubious sequenes which turned out to be false annotated should be stored. This file is used to delete sequences from the final reference database.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_duplicate_sequences.txt, Shows what species have sequences shared with other species.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_paraphyletic_group.txt, If sequences from a species creates a paraphyletic group, the species name is stored here.  
$RAWDIR/Database_curation/{date}/{date}_{run_number}_tree_string.newick, A newick file used to visualized the gene tree in programs such as iTOL and TreeViewer to detect false-annotated sequences.  

To effectively and reproducible detect and remove false-annotated sequences, we recommend to use a combination of the duplicate sequences, monophyletic and paraphyletic species name, combined with Web-Blast on both the aligned sequence aswell as accession number. We will highlight the importance, and benefits of using the counter when determining if the sequence is false annotated or not.

Example of a visualization of the Newick file.
<img width="629" height="203" alt="gene_tree_curation" src="https://github.com/user-attachments/assets/235d19ae-46f2-4417-b403-da18a0e1849c" />  

The genetree showcase a potential misannotation of Rana dalmatina with accession number MT872667.1.  

## Create the official database
### If this is the initial round:
Remember to store the accession numbers of the sequences you want to remove from the potential reference database in the $RAWDIR/Database_curation/{date}/{date}_{run_number}_sequences_to_delete.txt.
Congratulations, you are about to create your own reference database, choose a suitable name for it.

Example:
```
python Echopipe_database_completion.py -b BLAST_results/2026-03-17_1_to_curate.fasta -c Database_curation/2026-03-17_1/2026-03-17_1_aligned.fasta -u Amphibian_2026-03-17_1.fasta
```
### The following files are generated:
$RAWDIR/Amphibian_2026-03-17_1.fasta, the official reference database.
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_curated_tree_string.newick, A newick file used to visualized the gene tree in programs such as iTOL and TreeViewer to inspect intra- and interspecific variation within the marker region.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_duplicate_sequences.txt, Shows what species have sequences shared with other species.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_histogram_sequence_lengths.png, Shows and histogram of the length of the marker regions.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_histogram_sequences_per_species.png, Shows a histogram of the number of sequences per species.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_post_curation_monophyletic_group.txt, Shows what species forms a monophyletic group.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_post_curation_paraphyletic_group.txt, Shows what species forms a parphyletic group.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_aligned_Amphibian_2026-03-17_1.fasta, a MAFFT linsi alignment of the reference database.  
  

<img width="1000" height="600" alt="2026-03-17_1_histogram_sequence_lengths" src="https://github.com/user-attachments/assets/09dbd3c1-69d5-4b56-ad88-9b4a43f9dbdc" />  
  
Histogram over the distribution of sequence length in the database.    
<img width="1000" height="600" alt="2026-03-17_1_histogram_sequences_per_species" src="https://github.com/user-attachments/assets/427921ce-faf6-402f-95a6-428980b76d2a" />  
  
Histogram over the distribution of number of sequences per species in the database.    

### If older version exist:
Remember to store the accession numbers of the sequences you want to remove from the potential reference database in the $RAWDIR/Database_curation/{date}/{date}_{run_number}_sequences_to_delete.txt.  
Congratulations, you are about to update your reference database.  
This merges the counter of the earlier version(s) of your database.  
The second input is the new name of your database, and the third is the old version.  
```
python Echopipe_database_completion.py BLAST_results/2026-03-17_2_to_curate.fasta Amphibian_2026-03-17_2.fastaa --old_database Amphibian_2026-03-17_1.fasta  
```

## Evaluate database  
To futher evaluate the database, several tools are provided.  
The script include a monophyletic test, identicial sequences and more.  
```
python Echopipe_additional_evaluation.py Amphibian_2026-03-17_1.fasta Database_curation/2026-03-17_1/Curated_content/2026-03-17_1_post_curation_monophyletic_group.txt -f ACACCGCCCGTCACCCT -r GTAYACTTACCATGTTACGACTT  
```

Several dataframes are created and stored in $RAWDIR/Evaluation/Amphibian_2026-03-17_1/
$RAWDIR/Evaluation/Amphibian_2026-03-17_1/{Database_name}/{Database_name}_GC_content_histogram.png, Displays GC-content of the accepted reference sequences.  
$RAWDIR/Evaluation/Amphibian_2026-03-17_1/{Database_name}/{Database_name}_sequence_summary.csv, Gives a summary of the sequence.  
$RAWDIR/Evaluation/Amphibian_2026-03-17_1/{Database_name}/{Database_name}_evaluation_species_summary.csv, Gives a summary of the different species.  
$RAWDIR/Evaluation/Amphibian_2026-03-17_1/{Database_name}/{Database_name}_evaluation_family_summary.csv, Gives a summary of the different families.  

**Primer-mismatch visualization:**  

<img width="1200" height="800" alt="Forward_primer_nucleotide_proportions" src="https://github.com/user-attachments/assets/7c106843-b99d-4657-b403-0fd91b1b5796" />  

Visualization of primer-mismatch of the forward primer.  
  
<img width="1200" height="800" alt="Reverse_primer_nucleotide_proportions" src="https://github.com/user-attachments/assets/3e37237e-9f68-480e-b717-9a01760bd2e0" />  

Visualization of primer-mismatch of the reverse primer.

Additional primer amplification estimation is also provided.
