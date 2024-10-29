#!/usr/bin/env python3

"""
This script will merge the content of two databases and provide a new output file where new sequences have been appended and known sequences are counted.
"""

import warnings
from Bio import BiopythonDeprecationWarning
# Suppress the specific Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from Bio.Align.Applications import MafftCommandline # Let's us utilize MAFFT in Biopython as long as it has been installed.
from Bio.Phylo.Applications import FastTreeCommandline # Let's us utilize FastTree in Biopython as long as it has been installed.
from Bio import SeqIO
from ete3 import PhyloTree
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import argparse # Allows the program to be run from the command line.

parser = argparse.ArgumentParser(
    prog="EchoPipe - Database completion",
    description="A script for updating the current database with new content.",
    epilog="Version 1.0")

parser.add_argument('new_content', type=str,
    help="The path to the file with new content for the database.")
parser.add_argument('updated_database', type=str,
    help="The name of the new, updated database.")
parser.add_argument('-o', '--old_database', type=str, default="",
    help="The path to the old, existing database. Note: The reference template database is not intended for this.")

args = parser.parse_args()

new_db_content = os.path.abspath(args.new_content)
updated_database = os.path.abspath(args.updated_database)
old_database = os.path.abspath(args.old_database) if args.old_database else ""

run_name = os.path.basename(new_db_content).rsplit("_", maxsplit=2)[0] # The name of this run will be the same as the one that generated the input file.

###
###
###

# Adds existing sequences to a dictionary, appends new sequences and updates the counter for known sequences. Dictionary can be written to a file by calling the function write_output(output_path, species_data).
def process_fasta_file(file_path, species_data):
    for record in SeqIO.parse(file_path, "fasta"):
        header = f">{record.id}"
        sequence = str(record.seq)
        counter = int(header.split('|')[-1])
        species_key = header.split('|')[2].split(';')[-1]
        header_without_counter = header.rsplit('|', 1)[0]

        # Check if the sequence is already seen for the species
        if species_key in species_data:
            if sequence in species_data[species_key]['sequences']:
                # Update the counter for the specific sequence match
                species_data[species_key]['sequences'][sequence]['counter'] += counter
            else:
                # Append the new sequence for the species
                species_data[species_key]['sequences'][sequence] = {'counter': counter, 'header_without_counter': header_without_counter}
        else:
            species_data[species_key] = {'sequences': {sequence: {'counter': counter, 'header_without_counter': header_without_counter}}, 'sequence_set': {sequence}}



# Writes the content of the dictionary species_data to the defined output_path. Writes the content alphabetically based on species' names.
def write_output(output_path, species_data):
    with open(output_path, 'w') as output_file:
        # Sort species keys alphabetically
        sorted_species_keys = sorted(species_data.keys())
        
        for species_key in sorted_species_keys:
            data = species_data[species_key]
            
            # Include all sequences for the species
            for seq, seq_info in data['sequences'].items():
                counter = seq_info['counter']
                header_without_counter = seq_info['header_without_counter']
                output_file.write(f"{header_without_counter}|{counter}\n{seq}\n")


# Checks which accession numbers are to be removed based on the user's manual curation.
def read_sequences_to_remove(sequences_to_remove):
    with open(sequences_to_remove, "r") as file:
        return {line.strip() for line in file}


# This function removes the accession numbers and their corresponding sequences.
def remove_sequences(input_fasta, accession_numbers):
    temp_lines = []
    with open(input_fasta, "r") as infile:
        skip = False
        for line in infile:
            if line.startswith('>'):
                accession_number = line.split('|')[1]
                if accession_number in accession_numbers:
                    print(f"{accession_number} and its sequence has been removed.")
                    skip = True
                else:
                    skip = False
            if not skip:
                temp_lines.append(line)

    with open(input_fasta, "w") as outfile:
        for line in temp_lines:
            outfile.write(line)


# Filters away sequences that has other bases than ACGT (could, for example, be N).
def filter_sequences_in_place(input_fasta):
    records = list(SeqIO.parse(input_fasta, "fasta"))

    filtered_records = [record for record in records if all(base in "ACGT" for base in str(record.seq).upper())]

    with open(input_fasta, "w") as out_handle:
        SeqIO.write(filtered_records, out_handle, "fasta")


# Function to parse the fasta file and collect species counters using Biopython
def parse_fasta(file_path):
    species_counters = defaultdict(int)
    
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.description
        parts = header.split('|')
        counter = int(parts[-1])
        species_name = parts[2].split(';')[-1]
        species_counters[species_name] += counter
                
    return species_counters


# Function to plot a histogram showing how many sequences there are per species.
def plot_histogram_sequence_per_species(species_counters):
    values = sorted(species_counters.values(), reverse=True)
    
    # Calculate the mean and median.
    mean_value = np.mean(values)
    median_value = np.median(values)
    
    # Find the index of the value closest to the mean and median.
    mean_index = np.abs(np.array(values) - mean_value).argmin()
    median_index = np.abs(np.array(values) - median_value).argmin()

    plt.figure(figsize=(10, 6))
    
    # Create a list of indices for x-axis
    indices = list(range(1, len(values) + 1))
    
    # Plot the histogram
    plt.bar(indices, values, width=0.9)
    
    # Add lines for mean and median, rounded to one decimal place.
    plt.axvline(x=mean_index + 1, color='b', linestyle='dotted', linewidth=1.5, label=f'Mean: {mean_value:.1f} sequences')
    plt.axvline(x=median_index + 1, color='r', linestyle='dashed', linewidth=1.5, label=f'Median: {median_value:.1f} sequences')
    
    # Labeling the axes
    plt.xlabel('Species')
    plt.ylabel('Number of sequences')
    plt.title('Number of sequences per species')
    plt.grid(axis='y', linestyle='--', linewidth=0.7)
    
    plt.legend()
    plt.tight_layout()
    
    # Save the figure as a PNG file.
    plt.savefig(f'{run_name}_histogram_sequences_per_species.png')


# Function to plot a histogram showing the lengths of sequences and how many there are.
def plot_histogram_sequence_lenghts(updated_database):
    sequence_lengths = [len(record.seq) for record in SeqIO.parse(updated_database, "fasta")]

    mean_length = np.mean(sequence_lengths)
    median_length = np.median(sequence_lengths)

    plt.figure(figsize=(10, 6))

    # Create the histogram with each bar representing one specific sequence length.
    plt.hist(sequence_lengths, bins=range(min(sequence_lengths), max(sequence_lengths) + 2), width=1)

    # Add lines for mean and median, rounded to one decimal place.
    plt.axvline(mean_length, color='b', linestyle='dotted', linewidth=1.5, label=f'Mean: {mean_length:.1f} bp')
    plt.axvline(median_length, color='r', linestyle='dashed', linewidth=1.5, label=f'Median: {median_length:.1f} bp')

    # Plot customization.
    plt.title('Distribution of sequence lengths')
    plt.xlabel('Sequence length (bp)')
    plt.ylabel('Number of sequences')
    plt.grid(axis='y', linestyle='--', linewidth=0.7)

    plt.legend()
    plt.tight_layout()

    # Save the figure as a PNG file.
    plt.savefig(f'{run_name}_histogram_sequence_lengths.png')


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

species_data = {} # Creates an empty dictionary wherein sequences, their counters and accession numbers are kept track of.

# Fills up the dictionary species_data.
if old_database:
    process_fasta_file(old_database, species_data)
process_fasta_file(new_db_content, species_data)

# Creates the new database.
write_output(updated_database, species_data)

# Removes unwanted sequences, if any, from the new database.
accession_numbers = read_sequences_to_remove(f"Database_curation/{run_name}/{run_name}_sequences_to_delete.txt")
remove_sequences(updated_database, accession_numbers)

# Removes sequences that has nucleotides other than ACGT in them, such as N.
filter_sequences_in_place(updated_database)

###
###
###

curated_directory = f"Database_curation/{run_name}/Curated_content"
os.makedirs(curated_directory, exist_ok=True)
os.chdir(curated_directory)

print(f"Following informative files are found in the directory {curated_directory}\n")

duplicate_sequences = find_duplicates(updated_database)
if duplicate_sequences:
    duplicate_sequences_file = f"{run_name}_duplicate_sequences.txt"
    with open(duplicate_sequences_file, "w") as file:
        for sequence, headers in duplicate_sequences.items():
            for header in headers:
                file.write(header + "\n")
            file.write(sequence + "\n\n")
    print(f"Duplicate sequences that persisted after curation are noted in {duplicate_sequences_file}.\n")

# Creates histograms.
species_counters = parse_fasta(updated_database)
plot_histogram_sequence_per_species(species_counters)
plot_histogram_sequence_lenghts(updated_database)

temp_file = "temporary_file.fasta"
maffted_database = f"aligned_{args.updated_database}"
tree_string = f"{run_name}_curated_tree_string.newick"

with open(updated_database, "r") as outfile:
    with open(temp_file, "w") as infile:
        for line in outfile:
            infile.write(line.replace("|", "_").replace(";", "_"))

# Runs MAFFT with the "linsi" algorithm.
mafft_cline = MafftCommandline(
    input=temp_file,
    localpair=True,
    maxiterate=1000,
    thread=4,
    reorder=True
)

# Run Mafft
mafft_cline(stdout=maffted_database)


# Run FastTree
fasttree_cline = FastTreeCommandline(
    input=maffted_database,
    out=tree_string
)

fasttree_cline()

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

monophyletic_species = sorted(monophyletic_species)
monophyletic_species.insert(0, "These species came out as monophyletic in the gene tree after curation.\n\n")
with open(f"{run_name}_post_curation_monophyletic_group.txt", "w") as file: # This creates a text file of monophyletic entries.
    for mono in monophyletic_species:
        file.write(f"{mono.replace('_', ' ')}\n")

paraphyletic_species = sorted(paraphyletic_species)
paraphyletic_species.insert(0, "These species came out as paraphyletic in the gene tree after curation.\n\n")
with open(f"{run_name}_post_curation_paraphyletic_group.txt", "w") as file: # This creates a text file of paraphyletic entries.
    for para in paraphyletic_species:
        file.write(f"{para.replace('_', ' ')}\n")

os.remove(temp_file)

os.chdir("../../..")

###
###
###

# Moves the old_database from its current location to a directory dedicated to keeping old databases.
if old_database:
    path_to_old_databases = "Old_databases/"
    os.makedirs(path_to_old_databases, exist_ok=True)
    os.rename(old_database, path_to_old_databases + args.old_database)

print(f"\n\nThe reference database {updated_database} has successfully been created!\n")