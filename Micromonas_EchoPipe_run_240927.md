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

## Create an initial reference database, required input are:
Primer pairs, -f (forward) and -r (reverse) primer sequence.  
-l, Maximum length of amplicon retrived from NCBI.  
-t, minimum lengths of a sequences retrieved.  
-m, maximum number of sequences per species retrieved.  
-a, api.key from NCBI.  
-e, email address for NCBI to contact you.  
species list  

```
python Echopipe_reference_template_database_creator.py -f CCAGCASCYGCGGTAATTCC -r ACTTTCGTTCTTGATYRATGA -l 35000 -e eivisten@uio.no -t 200 -m 100 -a 86856fd70b62350246a86289c63698621607 micromonas_species_list.txt
cd $RAWDIR/Reference_template_creation

```

Inspect aligned_sequences_to_curate.fasta, located in $RAWDIR/Reference_template_creation. Visualize in MSA-tool, ex. Jalview. 
Remove clear dubious or wrongly annotated sequences from the aligned_sequences_to_curate.fasta file. The most important is that sequences from only the correct marker region is included.
Save the file with the same name, aligned_sequences_to_curate.fasta.

## Complete the creation and curation of the reference template database.
## A file called reference_template_database.fasta, is then created in the $RAWDIR.

```
cd $RAWDIR

python Echopipe_reference_template_database_creator.py -C 
```

## Use the script Echopipe_database_creation.py. Which can be iterated, by using the previously version of the database as a new reference_template_database. See below for example.
If this is the initial round:

```
python Echopipe_database_creation.py micromonas_species_list.txt reference_template_database.fasta -e eivisten@uio.no -a 86856fd70b62350246a86289c63698621607
```

### If older version exist: 
-c 20000, downloads up to 20 000 sequences from each species.
Required input file is then the old version of the database, and the species name. Increase the search by altering the -c or -l.

```
python Echopipe_database_creation.py micromonas_species_list.txt amphibase_240514.fasta -e your.email -a your.NCBIApiKey -c 20000
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
