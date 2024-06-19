#!/usr/bin/env python3

### This program is intended for creating dataframes to further assist with the evaluation of the reference database.


from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import os
import argparse

###
###
###

parser = argparse.ArgumentParser(
    prog="EchoPipe - Additional evaluation",
    description="Provides three dataframes from the reference database.",
    epilog="Version 1.0")

parser.add_argument('reference_database', type=str,
	help="The path to the reference database that the dataframes are going to be based on.")
parser.add_argument('monophyletic_group', type=str,
	help="The path to the text file containing the monophyletic groups that was generated when the reference database was completed.")

args = parser.parse_args()

reference_database = os.path.abspath(args.reference_database)
monophyletic_group = os.path.abspath(args.monophyletic_group)
dataframe_name = os.path.splitext(args.reference_database)[0]
evaluation_directory = os.path.join("Evaluation", dataframe_name)

os.makedirs(evaluation_directory, exist_ok=True)
os.chdir(evaluation_directory)

###
###
###

def find_duplicates(reference_database):
    sequence_dict = {}

    with open(reference_database, "r") as file:
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

###
###
###

with open(monophyletic_group, "r") as file:
	monophyletic_list = file.readlines()
monophyletic_list = [line.strip() for line in monophyletic_list]

duplicate_sequences = find_duplicates(reference_database)


dataframe_one_dictionary = {'superkingdom': [],
                        'phylum': [],
                        'class': [],
                        'order': [],
                        'family': [],
                        'genus': [],
                        'species': [],
                        'monophyletic': [],
                        'species_sequences': [],
                        'total_count': [],
                        'duplicate_sequences': []
                       }

dataframe_two_dictionary = {'superkingdom': [],
                        'phylum': [],
                        'class': [],
                        'order': [],
                        'family': [],
                        'genus': [],
                        'species': [],
                        'accession_number': [],
                        'sequence': [],
                        'count': [],
                        'length': []
                        }

lineage_keys = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

for record in SeqIO.parse(reference_database, "fasta"):
    split_header = record.id.split("|")
    lineage = split_header[2].split(";")
    species = lineage[6]
    accession = split_header[1]
    counter = int(split_header[3])


    if species not in dataframe_one_dictionary['species']:
        duplicate_sequence_number = 0
        for entries in duplicate_sequences.values():
            for entry in entries:
                if lineage[6] in entry:
                    duplicate_sequence_number += 1

        for key, lineage_value in zip(lineage_keys, lineage):
            dataframe_one_dictionary[key].append(lineage_value)

        dataframe_one_dictionary['monophyletic'].append("Yes" if species.replace('_', ' ') in monophyletic_list else "No")
        dataframe_one_dictionary['species_sequences'].append(1)
        dataframe_one_dictionary['total_count'].append(counter)
        dataframe_one_dictionary['duplicate_sequences'].append(duplicate_sequence_number)
    else:
        index = dataframe_one_dictionary['species'].index(species)
        dataframe_one_dictionary['total_count'][index] += counter
        dataframe_one_dictionary['species_sequences'][index] += 1

    for key, lineage_value in zip(lineage_keys, lineage):
        dataframe_two_dictionary[key].append(lineage_value)

    dataframe_two_dictionary['accession_number'].append(accession)
    dataframe_two_dictionary['sequence'].append(record.seq)
    dataframe_two_dictionary['count'].append(counter)
    dataframe_two_dictionary['length'].append(len(record.seq))


df_one = pd.DataFrame(dataframe_one_dictionary)
df_one.to_csv(f"{dataframe_name}_evaluation_species_summary.csv", sep=";", index=False)

df_two = pd.DataFrame(dataframe_two_dictionary)
df_two.to_csv(f"{dataframe_name}_sequence_summary.csv", sep=";", index=False)


# Convert monophyletic from string to boolean
df_one['monophyletic'] = df_one['monophyletic'].map({'Yes': True, 'No': False})

# Number of unique species per family
species_per_family = df_one.groupby('family')['species'].nunique().reset_index(name='number_of_species')

# Number of monophyletic species per family
monophyletic_species_per_family = df_one[df_one['monophyletic']].groupby('family')['species'].nunique().reset_index(name='number_of_monophyletic_species')

# Merge all the results into a single DataFrame
result = species_per_family.merge(monophyletic_species_per_family, on='family', how='left')

# Fill NaN values with 0 in the Number of Monophyletic Species column (in case there are no monophyletic species)
result['number_of_monophyletic_species'].fillna(0, inplace=True)
result['number_of_monophyletic_species'] = result['number_of_monophyletic_species'].astype(int)

# Add a new column that gives the number of paraphyletic species based on the total number of species minus the number of monophyletic species.
result['number_of_paraphyletic_species'] = result['number_of_species'] - result['number_of_monophyletic_species']

df_three = result

reconstructed_dict = {}

for sequence, headers in duplicate_sequences.items():
    for header in headers:
        family = header.split(';')[-3]
        reconstructed_dict.setdefault(sequence, []).append(family)

# Initialize a dictionary to count duplicates for each family
duplicate_counts = defaultdict(int)

# Iterate through the dictionary and count duplicates for each family
for sequences, families in reconstructed_dict.items():
    unique_families = set(families)  # Use a set to count each family only once per sequence
    for family in unique_families:
        duplicate_counts[family] += 1

# Insert duplicate counts into the dataframe
df_three['total_duplicate_sequences'] = df_three['family'].map(duplicate_counts).fillna(0).astype(int)

df_three.to_csv(f"{dataframe_name}_evaluation_family_summary.csv", sep=";", index=False)

os.chdir("../..")
print(f"\nDataframes have been created.\nThe dataframes are found in {os.path.abspath(evaluation_directory)}")