# Tutorial on how to use the code base EchoPipe
## Test run on Micromonas and Piredda primers
## Background
EchoPipe is an iterative, reproducible pipeline for creating, curating and evaluating reference databases for environmental DNA metabarcoding studies.  
By being an interative workflow, the workload, memory- and time usage is reduced, while coverage of the species of interest is increased.  
The entire pipeline is written in Python, and is distributed through Conda, which eases and streamlines the installation and usage. (For installation, see the Preparation step)  

## Requirements
### There are few required inputs for the workflow.
NCBI API key, needs to be obtained.  
Species list with scientific name.  
Email which NCBI may contact you.  
Forward and reverse primer sequence spanning the region of interest (both in 5'-3').  
Conda version 24.5.0

Check conda version with the following commands:

```
conda --version
conda 24.5.0 
``` 
If the results is not conda 24.5.0, please update your conda version.

![image](https://github.com/user-attachments/assets/28efafb2-4204-4685-9075-5410de525229)
Flowchart of the EchoPipe's workflow. Each color represents a module of the workflow and the corresponding script. Light-yellow: EchoPipe – Reference template database creator. Yellow: EchoPipe – Database creation. Orange: EchoPipe – Database curation. Red: EchoPipe – Database completion. Illustration is made with InkScape.

## Preparation
This tutorial assumes that the EchoPipe scripts is located in the RAWDIR, in this example the EchoPipe_Official/240927_PR2 path.

```
RAWDIR=EchoPipe_Official/240927_PR2
cd $RAWDIR
```
Download the conda environment for WSL, create it and activate it, using the following code.

```
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/refs/heads/main/environment_240930.yml
conda env create -f environment_240930.yml
conda activate EchoPipe 

```
## Create an initial template reference database.
### Options for Echopipe_reference_template_database_creator.py  
-f, forward primer sequence (5'-3') (required).  
-r, reverse primer sequence (5'-3') (required).  
-a, api.key from NCBI (required).  
-e, email address for NCBI to contact you (required).  
-q, query search word, target different marker regions using NCBI search term(s), ex. -q "AND 18s", default = NA (optional).  
-l, Maximum length of amplicon retrived from NCBI (optional).  
-t, minimum lengths of a sequences retrieved, default = 150 (optional).  
-m, maximum number of sequences per species retrieved, default = 1 (optional).  
-p, provided sequences, if a pre-existing fasta file is used as a reference template (optional).  
-z, longest amplicon length, Advanced setting. A value multiplied with the median length of the downloaded and trimmed sequences to discard sequences long sequences post trimming. Default is set to 2*median length species list (Optional).  
-C, Complete database, run after curating the database (shown below).  
-h, display help.


```
python Echopipe_reference_template_database_creator.py -f CCAGCASCYGCGGTAATTCC -r ACTTTCGTTCTTGATYRATGA -l 35000 -e your.email -t 200 -m 100 -a your.NCBIApiKey micromonas_species_list.txt
cd $RAWDIR/Reference_template_creation
```
    
    
Inspect aligned_sequences_to_curate.fasta with MSA-tool (Ex. Jalview).  
Remove clear dubious or wrongly annotated sequences from the aligned_sequences_to_curate.fasta file. 
No sequences should expand outside the primers.  
The most important is that sequences from only the correct marker region is included.  
Save the file with the same name, aligned_sequences_to_curate.fasta.  

### Output files:
aligned_sequences_to_curate.fasta, contains aligned versions of the reference sequences, use this for curation of reference template database.
filtered_aligned_sequences.fasta, contains aligned versions of the reference sequences (gaps removed).
non_approved_sequenes.txt, contains headers of sequences not accepted by the search terms.  
preformated_sequences.fasta, contains all non-aligned sequences analyzed.  
reference_template_database_inputs.txt, contains all the search terms used.   

### Example:
Open the aligned_sequences_to_curate.fasta, which is aligned using MAFFT (linsi) in Jalview, choose Colour by nucleotide to visualize sequence similarity.  
Scroll until you find a conserved region which overlaps with both primers, this ensure that the sequences in the reference template database is of the correct marker region.  
Keep in mind that gaps may be introduced as the alignment also contain potential reference sequences from other marker regions.  

<img width="952" alt="image" src="https://github.com/user-attachments/assets/8b8a3a89-e216-4122-8ac2-b3f48cf793f9">

Due to the reorder algorithm in the MSA, similiar sequences are stored together, and sequences from other marker regions are normally found at the bottom of the MSA.  
Highlighting and deleting (press delete) all sequences which are not similar to the primers in the file, is an effective way to remove sequences from other marker regions.  
Notice the similarity between the Forward primer and the non-selected sequences (white background color).  
Then save the file, and do not alter the path nor the file name.  
  
  
<img width="950" alt="image" src="https://github.com/user-attachments/assets/4d228472-4bd4-4788-831c-3e7cf8dd96f9">
  


### Complete the creation and curation of the reference template database.
A file called reference_template_database.fasta, is then created in the $RAWDIR.

```
cd $RAWDIR
python Echopipe_reference_template_database_creator.py -C 
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
ACTTCTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGT
```

### Options for Echopipe_database_creation.py   
-a, Api.key from NCBI (required).  
-e, Email address for NCBI to contact you (required).  
-s, Sort, if included, longer sequences are targeted (optional).  
-c, Max count, The maximum number of accession numbers downloaded per species. Default = 10 000.  
-q, Query search word, target different marker regions using NCBI search term(s), ex. -q "AND 18s".  
-l, Maximum length of amplicon retrived from NCBI, default=22 000 (optional).  
-t, Taxid, use this if last run of taxid wants to be used. Only works if this is not the first run (optional).  
-m, Mitochondria, the search targets mitochondrial sequences (optional).  
-r, Ribosomal, the search targets mitochondrial 12S ribosomal DNA (optional).  
-p, Provided sequences, if a pre-existing fasta file is used as a reference template (optional).  
-b, Batch size, The amount of sequences that can be downloaded simultaneously. Default = 5000 (optional).
-E, E-value, The E-value used for BLAST. The higher the value entered, the more stringent the BLAST becomes. Default = 5, (10^-5) (optional).
-R, Repeat, No new sequences are downloaded. Instead, all previously downloaded sequences are BLASTed against the new input database. Note: this options resets the file containing analysed accession numbers and the output will be curated from the start (optional).  
-h, Display help.

### If this is the initial round:

```
python Echopipe_database_creation.py micromonas_species_list.txt reference_template_database.fasta -e your.email -a your.NCBIApiKey
```
The scripts gives live updates, and tells the users what species are being processed.  
The scripts first sees if there are any new un-analyzed accession numbers on NCBI, and downloads if there are any.  
If there where any new un-analyzed accession numbers, the scripts then starts to download those.  
To check if the sequences related to the accession numbers are from the wanted marker region, a BLAST alignment against the initial reference template database is conducted.  
The longest match which are higher than the set E-value is kept, and are found in the $RAWDIR/BLAST_results/{date}_{run_number}_to_curate.fasta  

### Example:
Accession number download  
<img width="947" alt="image" src="https://github.com/user-attachments/assets/c55a90c2-d62f-4c89-9813-7a10fa6bc4aa">
  
Download of new un-analyzed accession numbers and a subsequent BLAST against the initial reference template database.  
  
<img width="940" alt="image" src="https://github.com/user-attachments/assets/0e29cc50-799a-4ed9-9e75-21764666e61d">


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
python Echopipe_database_creation.py micromonas_species_list.txt Micromas_24_09_30.fasta -e your.email -a your.NCBIApiKey -c 20000
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
python Echopipe_database_curation.py BLAST_results/2024-09-30_1_to_curate.fasta
```
  
<img width="533" alt="image" src="https://github.com/user-attachments/assets/05018847-e26e-4377-bbd6-df488c013d0b">
  
There location and usage of the produced tools.
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

<img width="553" alt="image" src="https://github.com/user-attachments/assets/1f37df6b-8615-4f32-9016-c688ab74af25">


<img width="485" alt="image" src="https://github.com/user-attachments/assets/9be211f6-3832-4be3-8888-cf37c5852da5">



## Create the official database
### If this is the initial round:
Remember to store the accession numbers of the sequences you want to remove from the potential reference database in the $RAWDIR/Database_curation/{date}/{date}_{run_number}_sequences_to_delete.txt.
Congratulations, you are about to create your own reference database, choose a suitable name for it.
```
python Echopipe_database_completion.py BLAST_results/2024-09-30_1_to_curate.fasta Micromas_24_09_30.fasta
```
  
<img width="944" alt="image" src="https://github.com/user-attachments/assets/4be3393f-bb1a-4d21-a4b5-aaf370122eac">  

  
### The following files are generated:
$RAWDIR/Micromas_24_09_30.fasta, the official reference database.
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_curated_tree_string.newick, A newick file used to visualized the gene tree in programs such as iTOL and TreeViewer to inspect intra- and interspecific variation within the marker region.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_duplicate_sequences.txt, Shows what species have sequences shared with other species.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_histogram_sequence_lengths.png, Shows and histogram of the length of the marker regions.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_histogram_sequences_per_species.png, Shows a histogram of the number of sequences per species.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_post_curation_monophyletic_group.txt, Shows what species forms a monophyletic group.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_post_curation_paraphyletic_group.txt, Shows what species forms a parphyletic group.  
$RAWDIR/Database_curation/{date}_{run_number}/Curated_content/{date}_{run_number}_aligned_Micromas_24_10_02.fasta, a MAFFT linsi alignment of the reference database.  
  
![2024-10-02_1_histogram_sequence_lengths](https://github.com/user-attachments/assets/2ae94f63-c253-4432-9cd8-739988e3d9b0)  

![2024-10-02_1_histogram_sequences_per_species](https://github.com/user-attachments/assets/00f32619-40fc-42d4-8905-282f47fb7604)  

   

### If older version exist:
Remember to store the accession numbers of the sequences you want to remove from the potential reference database in the $RAWDIR/Database_curation/{date}/{date}_{run_number}_sequences_to_delete.txt.
Congratulations, you are about to update your reference database.
This merges the counter of the earlier version(s) of your database.
The second input is the new name of your database, and the third is the old version.
```
python Echopipe_database_completion.py BLAST_results/2024-09-30_2_to_curate.fasta Micromas_24_10_02.fasta --old_database Micromas_24_09_30.fasta
```


## Evaluate database
To futher evaluate the database, several tools are provided.
The script include a monophyletic test, identicial sequences and more.
```
python Echopipe_additional_evaluation.py Micromas_24_09_30.fasta Database_curation/2024-09-30_1/Curated_content/2024-09-30_1_post_curation_monophyletic_group.txt
```

Several dataframes are created and stored in $RAWDIR/Evaluation/Micromas_24_10_02/{Database_name}/
$RAWDIR/Evaluation/Micromas_24_10_02/{Database_name}/{Database_name}_GC_content_histogram.png, Displays GC-content of the accepted reference sequences.  
$RAWDIR/Evaluation/Micromas_24_10_02/{Database_name}/{Database_name}_sequence_summary.csv, Gives a summary of the sequence.  
$RAWDIR/Evaluation/Micromas_24_10_02/{Database_name}/{Database_name}_evaluation_species_summary.csv, Gives a summary of the different species.  
$RAWDIR/Evaluation/Micromas_24_10_02/{Database_name}/{Database_name}_evaluation_family_summary.csv, Gives a summary of the different families.  

GC-content:  

![Micromas_24_10_02_GC_content_histogram](https://github.com/user-attachments/assets/c935d21e-4eae-4e2e-b96c-0d7827084433)  
  
  
Sequence summary:  
  
  
<img width="772" alt="image" src="https://github.com/user-attachments/assets/ec296794-08d2-49a3-bafa-c3a50b3834f8">  
  
  
Species summary:  
  
<img width="736" alt="image" src="https://github.com/user-attachments/assets/bef10dcb-824b-4563-9c73-2950fe5b4d0b">  
  
  
Family summary:  
  
<img width="452" alt="image" src="https://github.com/user-attachments/assets/965ad752-3cf5-4113-b06e-6607ee023f81">  
  
  
## Reformat database
To ensure the database is usable for your annotation method, a reformat script can be run, which creates a new copy of the database.
Available annotation methods are: 

sintax = SINTAX  
rdp = Ribosomal Database Project (RDP)  
dadt = DADA2 assignTaxonomy  
dads = DADA2 assignSpecies  
idt = IDTAXA  
qiime = QIIME 2  

Example:  
Create a DADA2 assignSpecies compatible database  
  
```
python Echopipe_reformat.py Micromas_24_09_30.fasta dads 
```
The files are then stored in $RAWDIR, with the annotation extention.
Example: 
Micromas_24_10_02_dads.fasta
