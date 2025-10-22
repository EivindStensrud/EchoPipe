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
# This program is intended to be used for curating reference databases, and is thus performed after that program has been executed.
#
# ---------------------------------------------------------------------------



import warnings
from Bio import BiopythonDeprecationWarning
# Suppress the specific Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from Bio.Align.Applications import MafftCommandline # Let's us utilize MAFFT in Biopython as long as it has been installed.
from Bio.Phylo.Applications import FastTreeCommandline # Let's us utilize FastTree in Biopython as long as it has been installed.
from Bio import SeqIO, AlignIO # Let's us parse through sequences easily to detect any unvalid nucleotides (i.e., anything other than ACTG).
from ete3 import PhyloTree # Used to make a phylogenetic tree out of our tree string (newick file).
import pandas as pd # If this program is run together with database creation this import is redundant. Enables us to have neat dataframes.
import time
from datetime import datetime # So that we can get the current day's date.
import os
import re
import argparse # Allows the program to be run from the command line.
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed # Used to parallelize script
from collections import defaultdict


###
### defs
###

# Filters away sequences that has other bases than ACGT (could, for example, be N). Creates a temporary output_file.
# With the else statement, sequences can be accepted if they have less than the number_ambigious_nucleotides defined by the user.
def filter_sequences(input_file, output_file, number_ambigious_nucleotides):
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)  # Turns the sequence into a string.
            if all(base in "ACGT" for base in sequence):  # Check if all bases are A, C, G, or T.
                out_handle.write(f">{record.id}\n")
                out_handle.write(str(record.seq) + "\n")
                #SeqIO.write(record, out_handle, "fasta") # This could be used as a single line, but causes issues downstreams with the formatting.
            else: 
                non_acgt_count = sum(1 for base in sequence if base not in "ACGT")
                if non_acgt_count <= number_ambigious_nucleotides:
                    out_handle.write(f">{record.id}\n")
                    out_handle.write(str(record.seq) + "\n")

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

# A function that prints a message (msg) and appends it on to a file.
def append_and_print_message(log_file, msg):
    with open(log_file, "a") as file:
        file.write(msg)
    print(msg)


def read_counter(counter_file, date=None):
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

def get_log_file_name(logs_dir):
    """Generate a unique log file name based on the latest date in the counter file and the counter."""
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

# Join all elements in sys.argv with a space separator
def get_command_string():
    command_string = ' '.join(sys.argv)
    return command_string

def number_threads():
    executor = ThreadPoolExecutor()
    # report the number of worker threads chosen by default
    thread_number = executor._max_workers
    thread_number = min((thread_number / 2 ), 7) # max 7, as API allows for up to 10 searches per second. 
    return thread_number

def split_fasta(input_file, sequences_per_file):
    # Create a directory to store split files
    output_dir = "split_fasta"
    os.makedirs(output_dir, exist_ok=True)
    
    # Read all sequences from the input FASTA file
    sequences = list(SeqIO.parse(input_file, "fasta"))
    total_sequences = len(sequences)
    
    # Dynamically set max sequences per file if specified number is greater than total
    if sequences_per_file > total_sequences:
        sequences_per_file = total_sequences
    
    num_files = (total_sequences + sequences_per_file - 1) // sequences_per_file


def align_batches(input_dir, output_dir):
    """Aligns each batch of FASTA files using MAFFT."""
    os.makedirs(output_dir, exist_ok=True)
    
    for file in os.listdir(input_dir):
        if file.endswith(".fasta"):
            input_file = os.path.join(input_dir, file)
            output_file = os.path.join(output_dir, f"aligned_{file}")

            # Create the MAFFT command line object
            mafft_cline = MafftCommandline(input=input_file, thread=7, auto=True)

            # Run MAFFT and capture the output
            try:
                stdout, stderr = mafft_cline()
                with open(output_file, "w") as out_handle:
                    out_handle.write(stdout)
                print(f"Aligned {input_file} -> {output_file}")
            except Exception as e:
                print(f"Error aligning {input_file}: {e}")

def merge_alignments(input_dir, output_file):
    with open(output_file, "w") as outfile:
        for file in sorted(os.listdir(input_dir)):
            if file.startswith("aligned_") and file.endswith(".fasta"):
                with open(os.path.join(input_dir, file)) as infile:
                    outfile.write(infile.read())
                print(f"Merged {file} into {output_file}")

def assess_alignment(file_path):
    alignment = AlignIO.read(file_path, "fasta")
    
    # Check overall alignment length and number of sequences
    print(f"Number of sequences: {len(alignment)}")
    print(f"Length of alignment: {alignment.get_alignment_length()}")
    
    # Calculate gaps and conservation
    gap_count = 0
    conserved_count = 0
    total_positions = alignment.get_alignment_length()

    for i in range(total_positions):
        column = alignment[:, i]
        if '-' in column:  # Count gaps
            gap_count += 1
        if len(set(column)) == 1:  # Count conserved positions
            conserved_count += 1

    print(f"Total gaps: {gap_count}")
    print(f"Total conserved positions: {conserved_count}")


def filter_unique_fasta(file_path):
    # Step 1: Count occurrences of each header
    header_count = defaultdict(int)
    
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.id
        header_count[header] += 1

    # Step 2: Create a temporary list for unique records
    unique_records = []

    for record in SeqIO.parse(file_path, "fasta"):
        if header_count[record.id] == 1:  # Only keep unique records
            unique_records.append(record)

    # Step 3: Overwrite the original file with unique records
    with open(file_path, "w") as output_handle:
        SeqIO.write(unique_records, output_handle, "fasta")

    print(f"Replaced the original file '{file_path}' with unique sequences.")




###
###
###
def main():


    parser = argparse.ArgumentParser(
        prog="EchoPipe - Database curation [-h] [-o OLD_DB] [-N NUM_NS] INPUT_FILE",
        description="A tool for curating reference databases.",
        epilog="Version 1.0")

    parser.add_argument('input_file', type=str, metavar='INPUT_FILE',
        help="Database to revise. Expects a .fasta file from BLAST_results/.")
    parser.add_argument('-o', '--old_database', type=str, default= "", metavar='OLD_DB',
        help="The previous version of database. Note: The reference template database is not intended for this.")
    parser.add_argument('-N', '--number_ns', type=int, default=0, metavar='NUM_NS',
        help="Number of N's and ambigious nucleotides allowed in a reference sequence, default = 0.")
    parser.add_argument('-M', '--mafft_online', type=str, default= "", metavar='MAFFT_online',
        help="If MAFFT online is used to align sequences. Show where the alignment are.")

    args = parser.parse_args()

    command_string = get_command_string()

    new_potential_sequences = os.path.abspath(args.input_file) # New sequences that are to be added to the existing ones and manually evaluated.
    most_recent_database = os.path.abspath(args.old_database) # Input by the user! Does not exist on the first run, but is used on subsequent runs after curation and evaluation, where dubious sequences should be removed.
    MAFFT_online_alignment = os.path.abspath(args.mafft_online) if args.mafft_online else None # Mafft online alignment
    number_ambigious_nucleotides = args.number_ns
    program_timer = time.time()
    thread_number = number_threads()

    filter_unique_fasta(new_potential_sequences)

    logs_dir = "Log_files/"
    log_file, run_name = get_log_file_name(logs_dir)


    database_curation_dir = "Database_curation/" # Directory where files associated with the database curation goes.
    working_directory = database_curation_dir + run_name


    counter_file = f"{logs_dir}run_counters.txt"  # File to store counters

    # Get current counter for the date, starting with 1 if it's a new date
    i = read_counter(counter_file)

    os.makedirs(working_directory, exist_ok=True)
    os.chdir(working_directory)

    concatenated_file = f"{run_name}_concatenated_file.fasta" # Fasta file containing the new, potential sequences and the previous ones.
    maffted_file = f"{run_name}_aligned.fasta" # Name of the initial output file. Used downstreams. NOTE: fasta format!
    tree_string = f"{run_name}_tree_string.newick" # Tree file that is generated and used to make our phylogenetic tree. NOTE: newick format!

    append_and_print_message(log_file,f"\n\nThe command used to run the script was: python {command_string}\n\n")
    

    if MAFFT_online_alignment and os.path.exists(MAFFT_online_alignment):
        print(f"MAFFT online allignment found.\n")

        concatenated_file = MAFFT_online_alignment
        temp_filtered_file = "temp_sequence_file.fasta"
        if os.path.exists(temp_filtered_file):
            print(f"Temp_file exist")
            # Count sequences in the concatenated file
            fh = open(temp_filtered_file)
            n = sum(1 for line in fh if line.startswith(">"))
            fh.close()
        else:
            filter_sequences(new_potential_sequences, temp_filtered_file, number_ambigious_nucleotides)
            # Count sequences in the concatenated file
            fh = open(temp_filtered_file)
            n = sum(1 for line in fh if line.startswith(">"))
            fh.close()

    else:
        print(f"MAFFT allignment conducted locally\n")
        # Filter sequences from the input FASTA file
        temp_filtered_file = "temp_sequence_file.fasta"
        filter_sequences(new_potential_sequences, temp_filtered_file, number_ambigious_nucleotides)

        new_sequence_marked_list = add_symbol_to_new_sequence(temp_filtered_file) # Calls the function, which adds a * at the end of each header in the file. This modification is stored in the list new_sequence_marked.

        # Concatenate most_recent_database and new_sequence_marked, i.e., put the previous sequences together with the new ones, into a list.
        try: # Use a 'try' block since most_recent_database may not be defined (if it is the first run).
            with open(most_recent_database, "r") as first_file: # Open and reads the file's content ...
                for line in first_file:
                    new_sequence_marked_list.append(line) # ... and appends it to the list new_sequence_marked_list.
        except: # In the case of a name error occurring, which would happen if most_recent_database is not defined, the rest of the program carries on.
            pass

        new_sequence_marked_list = [line.replace('|', '_').replace(';', '_').replace(' ','') for line in new_sequence_marked_list] # The list, now containing lines from most_recent_database and filtered new_potential_sequences, have two characters swapped out in each of its lines.

        with open(concatenated_file, "w") as second_file: # Opens a new file, concatenated_file to write onto.
            for concatenated_line in new_sequence_marked_list:
                second_file.write(concatenated_line) # Each concatenated_line from new_sequence_marked_list is written on to the file 'concatenated_file.fasta'.

        

        # Count sequences in the concatenated file
        fh = open(concatenated_file)
        n = sum(1 for line in fh if line.startswith(">"))
        fh.close()

        sequences_per_file_linsi = 1000
        sequences_per_file_auto = 10000

        #alignes using mafft linsi if less than 1000 fasta sequences in file. Otherwise using Mafft auto.
        if n <= sequences_per_file_linsi:
            # Create the MAFFT commandline object with the "linsi" algorithm.
            append_and_print_message(log_file,f"Align with MAFFT linsi using {n} sequences")
            mafft_cline = MafftCommandline(
                input=concatenated_file,
                localpair=True, # This gives us "linsi".
                maxiterate=100, # Set in order to not make it go on for too long.
                thread=thread_number, # Speeds up the process by running more processors at once. Number can be changed based on a computer's capabilities.
                reorder=True  # Makes it more presentable.
                )
        
        elif n <= sequences_per_file_auto:
            # Create the MAFFT commandline object with the "linsi" algorithm.
            append_and_print_message(log_file,f"Align with MAFFT auto using {n} sequences with {thread_number} threads")
            mafft_cline = MafftCommandline(
                input=concatenated_file,
                auto=True,
                thread=thread_number, # Speeds up the process by running more processors at once. Number can be changed based on a computer's capabilities.
                reorder=True  # Makes it more presentable.
                )
        
        else:
            append_and_print_message(log_file, f"Align with MAFFT auto using {n} sequences")
            mafft_cline = MafftCommandline(
                input=concatenated_file,
                auto=True

                )

        # Run MAFFT
        try:
            stdout, stderr = mafft_cline()
            # Process output to ensure sequences are in one line
            formatted_output = []
            current_sequence = []
            
            for line in stdout.splitlines():
                if line.startswith(">"):  # If it's a header line
                    if current_sequence:  # If there is a previous sequence, add it
                        formatted_output.append("".join(current_sequence))  # Join and add it to output
                        current_sequence = []  # Reset for new sequence
                    formatted_output.append(line)  # Add the header line
                else:
                    current_sequence.append(line.strip())  # Collect sequence lines

            # Add the last sequence to output if it exists
            if current_sequence:
                formatted_output.append("".join(current_sequence))
            
            # Write formatted output to file
            with open(maffted_file, "w") as out_handle:
                out_handle.write("\n".join(formatted_output) + "\n")
            print("Alignment complete.\n")
        except Exception as e:
            print(f"MAFFT error: {e}, try with larger RAM, or upload to {concatenated_file} MAFFT Online: https://mafft.cbrc.jp/alignment/server/large.html?aug31 \n")
            print(f"With large alignment, we recommend using MAFFT online with default settings, with following exceptions: \n")
            print(f"UPPERCASE/lowercase: Same as input")
            print(f"Strategy: PartTree")

    


    duplicate_sequences = find_duplicates(concatenated_file)
    if duplicate_sequences:
        duplicate_sequences_file = f"{run_name}_duplicate_sequences.txt"
        with open(duplicate_sequences_file, "w") as file:
            for sequence, headers in duplicate_sequences.items():
                for header in headers:
                    file.write(header + "\n")
                file.write(sequence + "\n\n")
        print(f"Duplicate sequences are noted in {duplicate_sequences_file}.\nCan be used for assist for viewing the gene tree and curating the database.\nPay attention to species names.\n")


    #fasttree different approach

    sequences_per_file_fastree_normal = 10000
    if n >= sequences_per_file_fastree_normal:
        # Create the MAFFT commandline object with the "linsi" algorithm.
        append_and_print_message(log_file,f"Running FastTree with fastest settings, input: {maffted_file} and output: {tree_string} with {n} sequences")
        fasttree_cline = FastTreeCommandline(
        nt=True,
        fastest=True,
        cat=20,
        input=maffted_file,
        out=tree_string)
    
    else:
        append_and_print_message(log_file,f"Running FastTree with standard settings, input: {maffted_file} and output: {tree_string} with {n} sequences")
        fasttree_cline = FastTreeCommandline(
        nt=True,
        input=maffted_file,
        out=tree_string)


    if not os.path.exists(maffted_file):
        print(f"Error: The input file '{maffted_file}' does not exist.")
    else:
        try:
            stdout, stderr = fasttree_cline()
        except Exception as e:
            print(f"Error when running FastTree: {e}")

    
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
        species_name = get_species_name(leaf, pattern)
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

    end_timestamp = datetime.today().strftime("%H:%M:%S") # Creates a variable called end_timestamp with the format HH:MM:SS.
    program_duration = round(time.time() - program_timer, 2)
    append_and_print_message(log_file,
        f"\nEchopipe_database_creation.py\n"
        f"It took {program_duration} seconds from start to finish.\n"
        f"The program finished at {end_timestamp}.\n"
        "##############################################################################\n\n")

    print("\nRecommendation:\n")
    print("Conduct a manual curation step, see the tutorial and manual for more thorough recommendations and instructions:\n")
    print("https://github.com/EivindStensrud/EchoPipe/blob/main/Tutorial.md\n")
    print("To finalize the reference database, choose a database name, we recommend to include organisms group and date:\n")
    print(f"Code line: \033[32mpython Echopipe_database_completion.py -b BLAST_results/{run_name}_to_curate.fasta -c Database_curation/{run_name}/{run_name}_aligned.fasta -u Database_name_{run_name}.fasta\033[0m\n")


    os.chdir("../..")

if __name__ == "__main__":
    main()