#!/usr/bin/env python3

# ---------------------------------------------------------------------------
# Copyright (c) 2024, Daniel Borg, Eivind Stensrud, Alexander Eiler
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#
# This script will merge the content of two databases and provide a new output file where new sequences have been appended and known sequences are counted.
# Read more on the GitHub repository: https://github.com/EivindStensrud/EchoPipe/tree/main
#
# ---------------------------------------------------------------------------

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
from datetime import datetime # So that we can get the current day's date.
import time # In order to add delays, so that servers are not put under much strain (and if one would like to see how long it takes for the code to execute).
import pathlib
import sys


###
###
###

def process_fasta_file(file_path, species_data, accessions_to_keep=None, log_file=None):
    for record in SeqIO.parse(file_path, "fasta"):
        
        # --- FILTERING LOGIC ---
        if '|' in record.id:
            accession_number = record.id.split('|')[1]
        elif '_' in record.id:
            accession_number = record.id.split('_')[1] 
        else:
            accession_number = record.id 

        if accessions_to_keep is not None and accession_number not in accessions_to_keep:
            print(f"Skipping unapproved sequence: {accession_number}")
            if log_file:
                with open(log_file, "a") as f:
                    f.write(f"{accession_number}\n") # Appends to _removed_accessions.txt
            continue 
        # -----------------------

        header = f">{record.id}"
        sequence = str(record.seq)
        
        if not all(base in "ACGT" for base in sequence.upper()):
            print(f"Skipping sequence with invalid bases (e.g., N): {accession_number}")
            if log_file:
                with open(log_file, "a") as f:
                    f.write(f"{accession_number}\n") # Appends to _removed_accessions.txt
            continue

        # --- NEW SAFE COUNTER & SPECIES EXTRACTION ---
        # 1. Remove the asterisk if it's there
        clean_header = header.rstrip('*')
        
        # 2. Extract the counter, species name, and base header based on the format
        if '|' in clean_header:
            counter_str = clean_header.split('|')[-1]
            species_name = clean_header.split('|')[2].split(';')[-1]
            header_without_counter = '|'.join(clean_header.split('|')[:-1])
        else:
            # If using the aligned file, everything is separated by underscores
            parts = clean_header.split('_')
            counter_str = parts[-1]
            # Genus and species are the two words right before the counter
            species_name = f"{parts[-3]}_{parts[-2]}"
            header_without_counter = '_'.join(parts[:-1])
            
        try:
            counter = int(counter_str)
        except ValueError:
            counter = 1 # Fallback just in case the format is totally unexpecte

        # --- MERGING LOGIC -N--
        if species_name in species_data:
            if sequence in species_data[species_name]['sequences']:
                species_data[species_name]['sequences'][sequence]['counter'] += counter
            else:
                species_data[species_name]['sequences'][sequence] = {
                    'counter': counter, 
                    'header_without_counter': header_without_counter
                }
        else:
            species_data[species_name] = {
                'sequences': {
                    sequence: {
                        'counter': counter, 
                        'header_without_counter': header_without_counter
                    }
                }
            }


def write_output(output_path, species_data):
# Writes the content of the dictionary species_data to the defined output_path. Writes the content alphabetically based on species' names.

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


def read_sequences_to_remove(file_path):
# Checks which accession numbers are to be removed based on the user's manual curation.
    to_remove = set()
    if file_path and os.path.exists(file_path):
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("Write down"): 
                    to_remove.add(line)
    return to_remove


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


def plot_histogram_sequence_per_species(species_counters, run_name):
# Function to plot a histogram showing how many sequences there are per species.
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


def plot_histogram_sequence_lenghts(updated_database, run_name):
# Function to plot a histogram showing the lengths of sequences and how many there are.

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


def find_duplicates(concatenated_fasta_file):
# Creates a dictionary with duplicate sequences and their headers.
    sequence_dict = {}
    with open(concatenated_fasta_file, "r") as file:
        current_header = None
        current_sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    # Remove dashes before storing in the dictionary
                    cleaned_sequence = current_sequence.replace("-", "")
                    if cleaned_sequence in sequence_dict:
                        sequence_dict[cleaned_sequence].append(current_header)
                    else:
                        sequence_dict[cleaned_sequence] = [current_header]
                current_header = line
                current_sequence = ""
            else:
                # Append line to current_sequence
                current_sequence += line
        
        # Handle the last sequence in the file
        if current_header and current_sequence:
            cleaned_sequence = current_sequence.replace("-", "")
            if cleaned_sequence in sequence_dict:
                sequence_dict[cleaned_sequence].append(current_header)
            else:
                sequence_dict[cleaned_sequence] = [current_header]
    # Identify duplicate sequences
    duplicate_sequences = {sequence: headers for sequence, headers in sequence_dict.items() if len(headers) > 1}
    return duplicate_sequences


def get_species_name(leaf, pattern):
    # Should account for subspecies names.
    name_parts = leaf.name
    parts = re.split(pattern, name_parts)
    parts = [part for part in parts if part]
    return '_'.join(parts[1].split('_')[7:-1])

def get_family_name(leaf, pattern):
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

def append_and_print_message(log_file, msg):
# A function that prints a message (msg) and appends it on to a file.

    with open(log_file, "a") as file:
        file.write(msg)
    print(msg)

def read_counter(counter_file, date=None):
#Read the counter value for the current date from the counter file.
    try:
        if os.path.exists(counter_file):
            with open(counter_file, "r") as f:
                data = f.read().strip()
                
                # Parse the counter values into a dictionary
                counters = dict(line.split(":") for line in data.splitlines())

                if date is not None:
                    # Return the counter for the given date, or 0 if not found
                    return int(counters.get(date, 0))
                else:
                    # If no date is provided, find the latest date
                    if counters:
                        latest_date = max(counters.keys(), key=lambda d: datetime.strptime(d, '%Y-%m-%d'))
                        return int(counters[latest_date])
                    else:
                        return 0
    except Exception as e:
        print(f"Error reading the counter file: {e}")
    
    # Default to counter 0 if file doesn't exist or error occurs
    return 0

def get_accession(header_line):
    header_line = header_line.strip()
    if header_line.startswith(">gb_"):
        # Aligned curated style: >gb_ACCESSION_...
        return header_line.split("_")[1]
    elif "|" in header_line:
        # BLAST style: >gb|ACCESSION|...
        return header_line.split("|")[1]
    else:
        # fallback
        return header_line.split()[0]


def get_log_file_name(logs_dir):
#Generate a unique log file name based on the latest date in the counter file and the counter.
    # Ensure the logs directory exists
    os.makedirs(logs_dir, exist_ok=True)

    # Path to the counter file
    counter_file = os.path.join(logs_dir, "run_counters.txt")
    
    # Find the latest date with an entry in the counter file
    latest_date = None
    if os.path.exists(counter_file):
        with open(counter_file, "r") as f:
            lines = f.read().strip().splitlines()
            if lines:
                latest_date = max(line.split(':')[0] for line in lines)
    
    if latest_date is None:
        # If there is no date in the file, use today's date
        latest_date = datetime.today().strftime('%Y-%m-%d')

    # Read the current counter for this date
    counter = read_counter(counter_file, latest_date)
    
    # Formulate the run name
    run_name = f"{latest_date}_{counter}"
    log_file_path = os.path.join(logs_dir, f"{run_name}_log.txt")
    log_file_path = os.path.abspath(log_file_path)

    
    return log_file_path, run_name

def extract_accession_numbers_from_aligned_file(fasta_file):
    accession_numbers = set()
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                accession = get_accession(line)
                accession_numbers.add(accession)
    return accession_numbers


def read_sequences_to_keep(file_path):
    to_keep = set()
    if file_path and os.path.exists(file_path):
        # Parses the aligned FASTA and grabs the accession numbers
        for record in SeqIO.parse(file_path, "fasta"):
            if '|' in record.id:
                accession = record.id.split('|')[1]
            elif '_' in record.id:
                accession = record.id.split('_')[1]
            else:
                accession = record.id
            to_keep.add(accession)
    return to_keep


# Join all elements in sys.argv with a space separator
def get_command_string():
    command_string = ' '.join(sys.argv)
    return command_string

def log_removed_accessions(accession_numbers, removed_accession_log, run_name):
    with open(removed_accession_log, "a") as file:
        for accession in accession_numbers:
            file.write(f"{accession.strip()}\n")  # Use .strip() to avoid issues with newline characters


###
###
###
def main():

    parser = argparse.ArgumentParser(
        prog="EchoPipe - Database completion",
        description="A script for updating the current database with new content.",
        epilog="Version 1.0")

    parser.add_argument('-b','--blast_file', type=str,
        help="The path to the new database BLAST file.")
    parser.add_argument('-c', '--curated_file', type=str,
        help="The path to the file with curated_alinged_fasta.")
    parser.add_argument('-u', '--updated_database', type=str,default="database.fasta",
        help="The name of the new, updated database.")
    parser.add_argument('-o', '--old_database', type=str, default="",
        help="The path to the old, existing database. Note: The reference template database is not intended for this.")

    args = parser.parse_args()

    new_db_content = os.path.abspath(args.blast_file)
    aligned_curated_file = os.path.abspath(args.curated_file)
    updated_database = os.path.abspath(args.updated_database)
    old_database = os.path.abspath(args.old_database) if args.old_database else ""
    thread_number = int(os.cpu_count() // 2 )

    command_string = get_command_string()

    logs_dir = "Log_files/"
    log_file, run_name = get_log_file_name(logs_dir)

    program_timer = time.time()

    curated_directory = f"Database_curation/{run_name}/Curated_content"

    os.makedirs(curated_directory, exist_ok=True)

    aligned_curated_file_format=f"Database_curation/{run_name}/{run_name}_aligned.fasta"
    removed_accession_log = f"Database_curation/{run_name}/{run_name}_removed_accessions.txt"

    append_and_print_message(log_file,f"\n\nThe command used to run the script was: python {command_string}\n\n")


    species_data = defaultdict(lambda: {'sequences': {}}) # Creates an empty dictionary wherein sequences, their counters and accession numbers are kept track of.

    accessions_to_delete = read_sequences_to_remove(aligned_curated_file)

    accessions_to_keep = read_sequences_to_keep(aligned_curated_file) if aligned_curated_file else None

    # 2. Process old database (-o) normally (we don't filter the master database)
    if old_database:
        print(f"Loading old database: {old_database}")
        process_fasta_file(old_database, species_data)

    # 3. Process the NEW sequences (-b) AND pass the 'keep' list to filter them
    if new_db_content:
        print(f"Loading and filtering new sequences: {new_db_content}")
        process_fasta_file(new_db_content, species_data, accessions_to_keep, removed_accession_log)
        
    # 4. Save the merged result to the updated database (-u)
    if updated_database:
        print(f"Writing updated database to: {updated_database}")
        write_output(updated_database, species_data)


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
    plot_histogram_sequence_per_species(species_counters,run_name)
    plot_histogram_sequence_lenghts(updated_database, run_name)

    temp_file = "temporary_file.fasta"
    temp_file = os.path.abspath(temp_file)
    maffted_database = f"aligned_{args.updated_database}"
    maffted_database = os.path.abspath(maffted_database)
    tree_string = f"{run_name}_curated_tree_string.newick"
    tree_string = os.path.abspath(tree_string)

    with open(updated_database, "r") as outfile:
        # Count the number of sequences
        n = sum(1 for line in outfile if line.startswith(">"))
        print(f"Number of sequences: {n}")  # Debug: Print the number of sequences

        outfile.seek(0) #start from first line in input file
        with open(temp_file, "w") as infile:
            for line in outfile:
                infile.write(line.replace("|", "_").replace(";", "_"))


    sequences_per_file_linsi = 1000
    sequences_per_file_auto = 10000

    #alignes using mafft linsi if less than 1000 fasta sequences in file. Otherwise using Mafft auto.
    if n <= sequences_per_file_linsi:
        # Create the MAFFT commandline object with the "linsi" algorithm.
        append_and_print_message(log_file,f"Align with MAFFT linsi using {n} sequences\n")
        mafft_cline = MafftCommandline(
            input=temp_file,
            localpair=True, # This gives us "linsi".
            maxiterate=100, # Set in order to not make it go on for too long.
            thread=thread_number, # Speeds up the process by running more processors at once. Number can be changed based on a computer's capabilities.
            reorder=True  # Makes it more presentable.
            )
    
    elif n <= sequences_per_file_auto:
        # Create the MAFFT commandline object with the "linsi" algorithm.
        append_and_print_message(log_file,f"Align with MAFFT auto using {n} sequences with {thread_number} threads\n")
        mafft_cline = MafftCommandline(
            input=temp_file,
            auto=True,
            thread=thread_number, # Speeds up the process by running more processors at once. Number can be changed based on a computer's capabilities.
            reorder=True  # Makes it more presentable.
            )
    
    else:
        append_and_print_message(log_file, f"Align with MAFFT auto using {n} sequences\n")
        mafft_cline = MafftCommandline(
            input=temp_file,
            auto=True,
            reorder=True

            )

    # Run MAFFT
    try:
        stdout, stderr = mafft_cline()
        with open(maffted_database, "w") as out_handle:
            out_handle.write(stdout)
        print("Alignment complete.\n")
    except Exception as e:
        print(f"MAFFT error: {e}, try with larger RAM, or upload to {temp_file} MAFFT Online: https://mafft.cbrc.jp/alignment/server/large.html?aug31 \n")
        print(f"With large alignment, we recommend using MAFFT online with default settings, with following exceptions: \n")
        print(f"UPPERCASE/lowercase: Same as input")
        print(f"Strategy: PartTree")


    duplicate_sequences = find_duplicates(maffted_database)
    if duplicate_sequences:
        duplicate_sequences_file = f"{run_name}_duplicate_sequences.txt"
        with open(duplicate_sequences_file, "w") as file:
            for sequence, headers in duplicate_sequences.items():
                for header in headers:
                    file.write(header + "\n")
                file.write(sequence + "\n\n")
        print(f"Duplicate sequences are noted in {duplicate_sequences_file}.\nCan be used for assist for viewing the gene tree and curating the database.\nPay attention to species names.\n")
    else:
        print("No duplicates sequenes found\n")

    #fasttree different approach

    sequences_per_file_fastree_normal = 10000
    if n >= sequences_per_file_fastree_normal:
        # Create the MAFFT commandline object with the "linsi" algorithm.
        append_and_print_message(log_file,f"Running FastTree with fastest settings: \ninput: {maffted_database} \noutput: {tree_string} with {n} sequences")
        fasttree_cline = FastTreeCommandline(
        nt=True,
        fastest=True,
        cat=20,
        input=maffted_database,
        out=tree_string)

    else:
        append_and_print_message(log_file,f"Running FastTree with standard settings: \ninput: {maffted_database} \noutput: {tree_string} with {n} sequences")
        fasttree_cline = FastTreeCommandline(
        nt=True,
        input=maffted_database,
        out=tree_string)


    if not os.path.exists(maffted_database):
        print(f"Error: The input file '{maffted_database}' does not exist.")
    else:
        try:
            stdout, stderr = fasttree_cline()
            if stderr:
                print("FastTree running:", stderr)
        except Exception as e:
            print(f"Error when running FastTree: {e}")


    tree = PhyloTree(tree_string)

    # Dictionary to store species and their associated leaves.
    species_leaves_dict = {}
    pattern = r'\.\d+_|[\d.]+'

    # Populate the dictionary with species and their leaves.
    for leaf in tree:
        species_name = get_species_name(leaf,pattern)
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
        families = set(get_family_name(leaf, pattern) for leaf in species_leaves)

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

    current_working_directory = os.getcwd()

    relative_path_database = os.path.relpath(updated_database, current_working_directory)

    end_timestamp = datetime.today().strftime("%H:%M:%S") # Creates a variable called end_timestamp with the format HH:MM:SS.
    program_duration = round(time.time() - program_timer, 2)
    append_and_print_message(log_file,
        f"\nIt took {program_duration} seconds from start to finish.\n"
        f"The program finished at {end_timestamp}.\n"
        f"The reference database {relative_path_database} has successfully been created!\n"
          "##############################################################################\n\n")




    print("\nRecommendation:\n")
    print("Conduct an additional evaluation of the database, see the tutorial and manual for more thorough recommendations and instructions:\n")
    print("https://github.com/EivindStensrud/EchoPipe/blob/main/Tutorial.md\n")
    print("Complete the database and conduct basic evaluation:\n")
    print(f"Code line:\n\033[32mpython Echopipe_additional_evaluation.py {relative_path_database} Database_curation/{run_name}/Curated_content/{run_name}_post_curation_monophyletic_group.txt\033[0m\n")
    print("To include evaluation of the primer binding sites, add the sequence of the forward and reverse primers using the -f (forward_primer_sequence) and -r (reverse_primer_sequence) flags:\n")
    print(f"Code line:\n\033[32mpython Echopipe_additional_evaluation.py {relative_path_database} Database_curation/{run_name}/Curated_content/{run_name}_post_curation_monophyletic_group.txt -f (forward_primer_sequence) -r (reverse_primer_sequence) \033[0m\n")


    os.chdir("../..")

if __name__ == "__main__":
    main()