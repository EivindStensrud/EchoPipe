#!/usr/bin/env python3

"""
This script will merge the content of two databases and provide a new output file where new sequences have been appended and known sequences are counted.
"""

import os
import argparse # Allows the program to be run from the command line.

parser = argparse.ArgumentParser(
    prog="ECHoPipe - Database updater",
    description="A script for updating the current database with new content.",
    epilog="Version 1.0")

parser.add_argument('new_content', type=str,
    help="The path to the file with new content for the database.")
parser.add_argument('updated_database', type=str,
    help="The name of the new, updated database.")
parser.add_argument('-o', '--old_database', type=str, default="",
    help="The path to the old, existing database. Note: The reference template database is not intended for this.")

args = parser.parse_args()

new_db_content = args.new_content
updated_database = args.updated_database
old_database = args.old_database

run_name = os.path.basename(new_db_content).rsplit("_", maxsplit=2)[0] # The name of this run will be the same as the one that generated the input file.

###
###
###

# Adds existing sequences to a dictionary, appends new sequences and updates the counter for known sequences. Dictionary can be written to a file by calling the function write_output(output_path, species_data).
def process_fasta_file(file_path, species_data):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith('>gb|'):
                header = lines[i].strip()
                sequence = lines[i + 1].strip()
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

                i += 2
            else:
                i += 1


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


# Moves the old_database from its current location to a directory dedicated to keeping old databases.
if old_database:
    path_to_old_databases = "Old_databases/"
    os.makedirs(path_to_old_databases, exist_ok=True)
    os.rename(old_database, path_to_old_databases + old_database)


# Removes unwanted sequences, if any, from the new database.
accession_numbers = read_sequences_to_remove(f"Database_curation/{run_name}/{run_name}_sequences_to_delete.txt")
remove_sequences(updated_database, accession_numbers)

print(f"\n\nThe reference database {updated_database} has successfully been created!\n")
