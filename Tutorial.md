# Tutorial on how to use the code base EchoPipe

## Preparation
This tutorial assumes that the EchoPipe scripts is located in the RAWDIR.
Find location of where you want the output to be stored, example.

```
RAWDIR=/EchoPipe/Amphibians
cd $RAWDIR
```

Activate the EchoPipe conda environment, more information in the README.md file.

```
conda activate EchoPipe 
```

## Create initial reference template database
Run Echopipe_reference_template_database_creator.py script to create a initial reference template database.

```
python Echopipe_reference_template_database_creator.py -f ACACCGCCCGTCACCCT -r GTAYACTTACCATGTTACGACTT -e your.email -t 50 -m 5 -a your.NCBIApiKey All_amphibian_names.txt

cd $RAWDIR/Reference_template_creation
```

Inspect aligned_sequences_to_curate.fasta, located in $RAWDIR/Reference_template_creation. Visualize in MSA-tool, ex. Jalview. 
Remove dubios or wrongly annotated sequences from the aligned_sequences_to_curate.fasta file.

Create a official reference_template_database.fasta, which will be used to find similiar sequences from NCBI.
```
cd $RAWDIR
python Echopipe_reference_template_database_creator.py -C 
```
## Create
Use the script Echopipe_database_creation.py. Which can be iterated, by using the previously version of the database as a new reference_template_database. See below for example.
```
python Echopipe_database_creation.py All_amphibian_names.txt reference_template_database.fasta -e your.email -a your.NCBIApiKey
```

## Curate
Use data created during this step to detect potential falsely annotated sequences.
Accession number of dubios sequences may be copied to {date}.dubious_sequences.txt file.
If the dubious sequences is deemed falsely annotated, add the accession number in the {date}.sequences_to_delete.txt file.

```
python Echopipe_database_curation.py BLAST_results/2024-05-13_1_to_curate.fasta #Name of file is decided on run date.
```

## Complete the database, and merge with older version of database. 
Give the database an official name, and congratulation, the database is now ready to use for species annotation.

```
python Echopipe_database_completion.py BLAST_results/2024-05-13_1_to_curate.fasta amphibase_240514.fasta 
```
## Evaluation
To evaluate the database, several tools are provided, and output is created.
The script include a monophyletic test, identicial sequences and more.

```
python Echopipe_additional_evaluation.py amphibase_240513.fasta Database_curation/2024-05-13_1/Curated_content/2024-05-13_1_post_curation_monophyletic_group.txt

```
## Reformat
To ensure the database is usable for your annotation method, a reformat script can be run, which creates a new copy of the database.
```
Echopipe_reformat.py amphibase_240513.fasta dads 
```