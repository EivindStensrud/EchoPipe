#!/usr/bin/env python3

# This program is intended to be used for curating reference databases, and is thus performed after that program has been executed.
# I do not take credit for most of the following code as it has been provided to me by other creative sources.

# NOTE: MAFFT must be installed! Was downloaded from: https://mafft.cbrc.jp/alignment/software/linux.html
#       The downloaded file, for Ubuntu, had the .deb extension and was installed in bash with: sudo dpkg -i mafft_<version>.deb
# NOTE: FastTree must be installed! In base: conda install bioconda::fasttree
# NOTE: ETE3 must be installed! In bash: pip install ete3


from Bio.Align.Applications import MafftCommandline # Let's us utilize MAFFT in Biopython as long as it has been installed.
from Bio.Phylo.Applications import FastTreeCommandline # Let's us utilize FastTree in Biopython as long as it has been installed.
from Bio import SeqIO # Let's us parse through sequences easily to detect any unvalid nucleotides (i.e., anything other than ACTG).
from ete3 import PhyloTree # Used to make a phylogenetic tree out of our tree string (newick file).
import pandas as pd # If this program is run together with database creation this import is redundant. Enables us to have neat dataframes.
import time
import os
import re
import argparse # Allows the program to be run from the command line.

###
###
###

parser = argparse.ArgumentParser(
    prog="EchoPipe - Database curation",
    description="A tool for curating reference databases by creating phylogenetic trees to be reviewed.",
    epilog="Version 1.0")

parser.add_argument('input_file', type=str,
    help="The database that is to be reviewed. Expects a .fasta file from BLAST_results/.")
parser.add_argument('-o', '--old_database', type=str, default="",
    help="The most recent database that is to be updated in this iteration. Note: The reference template database is not intended for this.")

args = parser.parse_args()

new_potential_sequences = os.path.abspath(args.input_file) # New sequences that are to be added to the existing ones and manually evaluated.
most_recent_database = os.path.abspath(args.old_database) # Input by the user! Does not exist on the first run, but is used on subsequent runs after curation and evaluation, where dubious sequences should be removed.

program_timer = time.time()


###
###
###

# Filters away sequences that has other bases than ACGT (could, for example, be N). Creates a temporary output_file.
def filter_sequences(input_file, output_file):
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)  # Turns the sequence into a string.
            if all(base in "ACGT" for base in sequence):  # Check if all bases are A, C, G, or T.
                out_handle.write(f">{record.id}\n")
                out_handle.write(str(record.seq) + "\n")
                #SeqIO.write(record, out_handle, "fasta") # This could be used as a single line, but causes issues downstreams with the formatting.
            else:
                print(f"The accession number {record.id.split('|')[1]} included bases other than ACGT in its sequence and is therefore omitted.")

# Add * to new_potential_sequences in order to be able to distinguish new sequences from previously curated ones when evaluating the tree.
def add_symbol_to_new_sequence(new_potential_sequences):
    new_sequence_marked = [] # When the function is invoked an empty list called new_sequence_marked is created.
    with open(new_potential_sequences, 'r') as new_sequence: # Reads through the file...
        for line in new_sequence: # ... line by line...
            if line.startswith('>'): # ... and if a line starts with '>', then a * and newline is added to the end. This means each header receives a * at the end so that we may identify them more easily.
                updated_line = line.strip() + "*\n"
                new_sequence_marked.append(updated_line)
            else: # The unaltered line (i.e., the sequence) is added to the list as well.
                new_sequence_marked.append(line)
    return new_sequence_marked

# Creates a dictionary with duplicate sequences and their headers.
def find_duplicates(concatenated_fasta_file):
    sequence_dict = {}

    with open(concatenated_fasta_file, "r") as file:
        current_header = None
        current_sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    if current_sequence in sequence_dict:
                        sequence_dict[current_sequence].append(current_header)
                    else:
                        sequence_dict[current_sequence] = [current_header]
                current_header = line
                current_sequence = ""
            else:
                current_sequence += line
        if current_header and current_sequence:
            if current_sequence in sequence_dict:
                sequence_dict[current_sequence].append(current_header)
            else:
                sequence_dict[current_sequence] = [current_header]
    duplicate_sequences = {sequence: headers for sequence, headers in sequence_dict.items() if len(headers) > 1}
    return duplicate_sequences

def get_species_name(leaf):
    # Should account for subspecies names.
    name_parts = leaf.name
    parts = re.split(pattern, name_parts)
    parts = [part for part in parts if part]
    return '_'.join(parts[1].split('_')[7:-1])

def get_family_name(leaf):
    # Assuming the family identifier has index 6.
    name_parts = leaf.name
    parts = re.split(pattern, name_parts)
    parts = [part for part in parts if part]
    return parts[1].split('_')[5]

def is_monophyletic(tree, species_leaves):
    # Find the common ancestor of the specified species leaves.
    common_ancestor = tree.get_common_ancestor(species_leaves)
    
    # Check if the common ancestor is an internal node (non-terminal).
    if set(common_ancestor.get_leaf_names()) == set(leaf.name for leaf in species_leaves):
        return True
    elif len(species_leaves) == 1 and species_leaves[0].is_leaf():
        # Special case: species with only one sequence
        return True
    else:
        return False

###
###
###

database_curation = "Database_curation/" # Directory where files associated with the database curation goes.
run_name = os.path.basename(args.input_file).rsplit("_", maxsplit=2)[0] # The name of this run will be the same as the one that generated the input file.
working_directory = database_curation + run_name


os.makedirs(working_directory, exist_ok=True)
os.chdir(working_directory)

temp_filtered_file = "temp_sequence_file.fasta"
filter_sequences(new_potential_sequences, temp_filtered_file)

concatenated_file = f"{run_name}_concatenated_file.fasta" # Fasta file containing the new, potential sequences and the previous ones.
maffted_file = f"{run_name}_aligned.fasta" # Name of the initial output file. Used downstreams. NOTE: fasta format!
tree_string = f"{run_name}_tree_string.newick" # Tree file that is generated and used to make our phylogenetic tree. NOTE: newick format!

###
###
###

new_sequence_marked_list = add_symbol_to_new_sequence(temp_filtered_file) # Calls the function, which adds a * at the end of each header in the file. This modification is stored in the list new_sequence_marked.

# Concatenate most_recent_database and new_sequence_marked, i.e., put the previous sequences together with the new ones, into a list.
try: # Use a 'try' block since most_recent_database may not be defined (if it is the first run).
    with open(most_recent_database, "r") as first_file: # Open and reads the file's content ...
        for line in first_file:
            new_sequence_marked_list.append(line) # ... and appends it to the list new_sequence_marked_list.
except: # In the case of a name error occurring, which would happen if most_recent_database is not defined, the rest of the program carries on.
    pass

new_sequence_marked_list = [line.replace('|', '_').replace(';', '_') for line in new_sequence_marked_list] # The list, now containing lines from most_recent_database and filtered new_potential_sequences, have two characters swapped out in each of its lines.

with open(concatenated_file, "w") as second_file: # Opens a new file, concatenated_file to write onto.
    for concatenated_line in new_sequence_marked_list:
        second_file.write(concatenated_line) # Each concatenated_line from new_sequence_marked_list is written on to the file 'concatenated_file.fasta'.



# Create the MAFFT commandline object with the "linsi" algorithm.
mafft_cline = MafftCommandline(
    input=concatenated_file,
    localpair=True, # This gives us "linsi".
    maxiterate=1000, # Set in order to not make it go on for too long.
    thread=4, # Speeds up the process by running more processors at once. Number can be changed based on a computer's capabilities.
    reorder=True # Makes it more presentable.
)

# Run MAFFT
mafft_cline(stdout=maffted_file)

print("Alignment complete.\n") # Optional message that is displayed once this process has finished.

duplicate_sequences = find_duplicates(concatenated_file)
if duplicate_sequences:
    duplicate_sequences_file = f"{run_name}_duplicate_sequences.txt"
    with open(duplicate_sequences_file, "w") as file:
        for sequence, headers in duplicate_sequences.items():
            for header in headers:
                file.write(header + "\n")
            file.write(sequence + "\n\n")
    print(f"Duplicate sequences are noted in {duplicate_sequences_file}.\nCan be used for assist for viewing the gene tree and curating the database.\nPay attention to species names.\n")


# Run FastTree to generate a Newick tree file.
fasttree_cline = FastTreeCommandline(
    nt=True,
    input=maffted_file,
    out=tree_string
)

fasttree_cline()

#from Bio import Phylo # Phylo is used for things such as reading, manipulating and visualizing phylogenetic trees.
# Parse the output tree file using Biopython.
# tree = Phylo.read(tree_string, "newick") # Optional. Requires Phylo to be imported.

# Print the tree
# Phylo.draw_ascii(tree) # Optional. Requires Phylo to be imported.

tree = PhyloTree(tree_string)

# Dictionary to store species and their associated leaves.
species_leaves_dict = {}
pattern = r'\.\d+_|[\d.]+'

# Populate the dictionary with species and their leaves.
for leaf in tree:
    species_name = get_species_name(leaf)
    if species_name not in species_leaves_dict:
        species_leaves_dict[species_name] = []
    species_leaves_dict[species_name].append(leaf)


# Check if each species forms a monophyletic group.
monophyletic_species = []
paraphyletic_species = []
for species_name, species_leaves in species_leaves_dict.items():
    result = is_monophyletic(tree, species_leaves)
    if result == True:
        monophyletic_species.append(species_name)
    else:
        paraphyletic_species.append(species_name)

table_data = {'Species': [], 'Family': [], 'Is Monophyletic': []}

# Check if each set of species forms a paraphyletic group and record family information.
for species_name, species_leaves in species_leaves_dict.items():
    is_monophyletic_result = is_monophyletic(tree, species_leaves)
    families = set(get_family_name(leaf) for leaf in species_leaves)

    for family in families:
        table_data['Species'].append(species_name)
        table_data['Family'].append(family)
        table_data['Is Monophyletic'].append(is_monophyletic_result)

# Create a DataFrame from the table data.
df = pd.DataFrame(table_data)

# Display the DataFrame. Optional.
# print(df)

if monophyletic_species:
    monophyletic_species = sorted(monophyletic_species)
    monophyletic_species.insert(0, "These species came out as monophyletic in the gene tree. Still, consider reviewing branch length to for any oddities.\n\n")
    with open(f"{run_name}_monophyletic_group.txt", "w") as file: # This creates a text file of monophyletic entries.
        for mono in monophyletic_species:
            file.write(f"{mono.replace('_', ' ')}\n")

if paraphyletic_species:
    paraphyletic_species = sorted(paraphyletic_species)
    paraphyletic_species.insert(0, "These species came out as paraphyletic in the gene tree. Consider reviewing them. Feel free to make a note in this document if you deem the entries to be acceptable.\n\n")
    with open(f"{run_name}_paraphyletic_group.txt", "w") as file: # This creates a text file of paraphyletic entries.
        for para in paraphyletic_species:
            file.write(f"{para.replace('_', ' ')}\n")
        
df.to_csv(f"{run_name}_dataframe.csv", index=False) # df (Dataframe) is saved as a .csv file. The index=False parameter ensures that the DataFrame index is not written to the file.


# Two files are created that are associated with the run. The user is meant to manually type in accession numbers that are dubious and those that should be deleted from the database (after curation!)
with open(f"{run_name}_dubious_sequences.txt", "w") as file:
    file.write("Write down accession numbers below - one per line - of those sequences with dubious results based on manual curation.\n")

with open(f"{run_name}_sequences_to_delete.txt", "w") as file:
    file.write("Write down accession numbers below - one per line - of the dubious sequences that remain dubious based on manual curation.\n")

os.remove(temp_filtered_file)

program_duration = round(time.time() - program_timer, 2)
print(f"\nIt took {program_duration} seconds from start to finish.\n")

os.chdir("../..")