# Tutorial on how to use the code base EchoPipe
## Test run on micromonas and piredda primers

## Preparation
This tutorial assumes that the EchoPipe scripts is located in the RAWDIR, in this example the EchoPipe_Official/240927_PR2 path

```
RAWDIR=EchoPipe_Official/240927_PR2
cd $RAWDIR
```
Download the conda environment from WSL, create it and activate it.

```
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/refs/heads/main/environment_240930.yml
conda env create -f environment_240930.yml
conda activate EchoPipe 

```
## Create an initial template reference database.
### Options for Echopipe_reference_template_database_creator.py  
-f, forward primer sequence (required).  
-r, reverse primer sequence (required).  
-a, api.key from NCBI (required).  
-e, email address for NCBI to contact you (required).  
-q, query search word, target different marker regions using NCBI search term(s), ex. -q "AND 18s.  
-l, Maximum length of amplicon retrived from NCBI (optional).  
-t, minimum lengths of a sequences retrieved (optional).  
-m, maximum number of sequences per species retrieved (optional).  
-p, provided sequences, if a pre-existing fasta file is used as a reference template (optional).  
-z, longest amplicon length, Advanced setting. A value multiplied with the median length of the downloaded and trimmed sequences to discard sequences long sequences post trimming. Default is set to 2*median length species list.  
-C, Complete database, run after curating the database (shown below).  
-h, display help.


```
python Echopipe_reference_template_database_creator.py -f CCAGCASCYGCGGTAATTCC -r ACTTTCGTTCTTGATYRATGA -l 35000 -e eivisten@uio.no -t 200 -m 100 -a 86856fd70b62350246a86289c63698621607 micromonas_species_list.txt
cd $RAWDIR/Reference_template_creation
```
    
    
Inspect aligned_sequences_to_curate.fasta, located in $RAWDIR/Reference_template_creation. Visualize in MSA-tool, ex. Jalview.   
Remove clear dubious or wrongly annotated sequences from the aligned_sequences_to_curate.fasta file.  
The most important is that sequences from only the correct marker region is included.  
Save the file with the same name, aligned_sequences_to_curate.fasta.  
Additional, all fasta files are stored in the filtered_aligned_sequences.fasta, are fi

### Example:
Open the aligned_sequences_to_curate.fasta, which is aligned using MAFFT (linsi) in Jalview, choose Colour by nucleotide to visualize sequence similarity.  
Scroll until you find a conserved region which overlaps with both primers, this ensure that the sequences in the reference template database is of the correct marker region.  
Keep in mind that gaps may be introduced as the alignment also contain potential reference sequences from other marker regions.  

<img width="952" alt="image" src="https://github.com/user-attachments/assets/8b8a3a89-e216-4122-8ac2-b3f48cf793f9">

Due to the reorder algorithm in the MSA, dissimilar sequences sorts similiar sequences together, and sequences from other marker regions is visualized at the bottom of the MSA.  
Highlighting and deleting (press delete) all sequences which are not similar to the primers in the file, is an effective way to remove sequences from other marker regions.  
Notice the similarity between the Forward primer and the non-selected sequences (white background color).  
Then save the file, and do not alter the path nor the file name.  

<img width="950" alt="image" src="https://github.com/user-attachments/assets/4d228472-4bd4-4788-831c-3e7cf8dd96f9">

Additional, some other files are created.  
filtered_aligned_sequences.fasta contains aligned versions of the reference sequences (gaps removed).
non_approved_sequenes.txt contains headers of sequences not accepted by the search terms.  
preformated_sequences.fasta contains all non-aligned sequences analyzed.  
reference_template_database_inputs.txt contains all the search terms used.  

### Complete the creation and curation of the reference template database.
A file called reference_template_database.fasta, is then created in the $RAWDIR.

```
cd $RAWDIR
python Echopipe_reference_template_database_creator.py -C 
```

## Create an initial reference database  
## Use the script Echopipe_database_creation.py  
This script is iterative, and could be re-run by changing the search options.
### Options for Echopipe_database_creation.py   
-a, api.key from NCBI (required).  
-e, email address for NCBI to contact you (required).  
-q, query search word, target different marker regions using NCBI search term(s), ex. -q "AND 18s.  
-l, Maximum length of amplicon retrived from NCBI, default=22 000 (optional).  
-t, taxid, use this if last run of taxid wants to be used. Only works if this is not the first run (optional).  
-m, maximum number of sequences per species retrieved (optional).  
-p, provided sequences, if a pre-existing fasta file is used as a reference template (optional).  
-z, longest amplicon length, Advanced setting. A value multiplied with the median length of the downloaded and trimmed sequences to discard sequences long sequences post trimming. Default is set to 2*median length species list.  
-h, display help.

If this is the initial round:

```
python Echopipe_database_creation.py micromonas_species_list.txt reference_template_database.fasta -e eivisten@uio.no -a 86856fd70b62350246a86289c63698621607
```

### If older version exist: 
-c 20000, downloads up to 20 000 sequences from each species.
Required input file is then the old version of the database, and the species name. Increase the search by altering the -c or -l.

```
python Echopipe_database_creation.py micromonas_species_list.txt Micromas_24_09_30.fasta -e your.email -a your.NCBIApiKey -c 20000
```

## Curate database
Use data created during this step to detect potential falsely annotated sequences.
Accession number of dubious sequences may be copied to {date}.dubious_sequences.txt file.
If the dubious sequences is deemed falsely annotated, add the accession number in the {date}.sequences_to_delete.txt file.

```
python Echopipe_database_curation.py BLAST_results/2024-09-30_1_to_curate.fasta
```

# Create the official database, give an official name, this case. Micromonas_24_09_31.fasta
## If this is the initial round:
```
python Echopipe_database_completion.py BLAST_results/2024-09-30_1_to_curate.fasta Micromas_24_09_30.fasta
```

## If older version exist:
```
python Echopipe_database_completion.py BLAST_results/2024-09-30_2_to_curate.fasta Micromas_24_09_30.fasta --old_database Micromas_24_09_30.fasta
```
## Evaluate database
To evaluate the database, several tools are provided, and output is created.
The script include a monophyletic test, identicial sequences and more.
```
python Echopipe_additional_evaluation.py Micromas_24_09_30.fasta Database_curation/2024-09-30_1/Curated_content/2024-09-30_1_post_curation_monophyletic_group.txt
```

## Reformat database
To ensure the database is usable for your annotation method, a reformat script can be run, which creates a new copy of the database.
```
Echopipe_reformat.py Micromas_24_09_30.fasta dads 
```
