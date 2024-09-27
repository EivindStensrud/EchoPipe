#!/usr/bin/env python3

### This program is intended for creating dataframes to further assist with the evaluation of the reference database.


from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import tempfile
import argparse


###
###
###

parser = argparse.ArgumentParser(
    prog="EchoPipe - Additional evaluation",
    description="Provides additional dataframes from the reference database, with the option to also view primer compatibility - either one or both directions at the same time.",
    epilog="Version 1.0")

parser.add_argument('reference_database', type=str,
    help="The path to the reference database that the dataframes are going to be based on.")
parser.add_argument('monophyletic_group', type=str,
    help="The path to the text file containing the monophyletic groups that was generated when the reference database was completed.")
parser.add_argument('-f', '--forward_primer', type=str, default=None,
    help="The forward primer that wishes to be checked against the reference database's sequence (5'-3' orientation).")
parser.add_argument('-r', '--reverse_primer', type=str, default=None,
    help="The reverse primer that wishes to be checked against the reference database's sequence (5'-3' orientation).")

args = parser.parse_args()

reference_database = os.path.abspath(args.reference_database)
monophyletic_group = os.path.abspath(args.monophyletic_group)
database_name = os.path.splitext(args.reference_database)[0]
evaluation_directory = os.path.join("Evaluation", database_name)
forward_primer = args.forward_primer
reverse_primer = args.reverse_primer

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

# Get a list of positions that is used for checking mismatches between the aligned primer and sequence.
def get_positions(sequence_input):
    return list(range(len(sequence_input)))

# Shows if there is a mismatch at a particular position.
def positions_to_check(mode, positions, aligned_primer, aligned_sequence, record_result):
    number_of_mismatches = 0
    length = len(positions)
    mismatch_details = []

    number_of_mismatches_mode = f"{mode}_number_of_mismatches"
    mismatches_mode = f"{mode}_mismatches"
    
    for pos in positions:
        if pos < length:
            sequence_char = aligned_sequence[pos]
            primer_char = aligned_primer[pos]
            if sequence_char != primer_char:
                number_of_mismatches += 1
                bad_position = length - pos 
                mismatch_details.append((bad_position, primer_char, sequence_char))
        else:
            print(f"Position {pos+1} is out of bounds for the aligned sequences.")

    record_result[number_of_mismatches_mode] = number_of_mismatches
    if mismatch_details:
        record_result[mismatches_mode] = mismatch_details
    else:
        record_result[mismatches_mode] = "NA"

# Align the primer with the sequence and check for mismatches.
def primer_alignment(mode, primer_sequence, sequence_to_align, record_result):
    aligner = PairwiseAligner()
    aligner.match = 1
    aligner.mismatch = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    aligner.end_gap_score = 0
    #aligner.target_end_gap_score = 0
    #aligner.query_end_gap_score = 0
    #aligner.target_internal_open_gap_score = -50 # Can be used to prevent gaps inside of the primer.
    aligner.mode = 'global'
    
    alignments = aligner.align(primer_sequence, sequence_to_align)

    primer_length = len(primer_sequence)

    max_score = aligner.match * primer_length
    score_cutoff = max_score * 0.5
    length_cutoff = primer_length * 1.3

    sequence_key = f"{mode}_sequence"
    
    best_alignment = None # Implemented for the purpose of ensuring that only one alignment per accession number is used. As it is now, it will omit alignments that share the same score and only keep one of them.
    best_score = -float('inf')  # Initialize with a very low score.

    for alignment in alignments:
        if alignment.score > score_cutoff and length_cutoff >= alignment.length:
            if alignment.score > best_score:
                best_score = alignment.score
                best_alignment = alignment

    if best_alignment:
        aligned_primer = alignment[0]
        
        # Find the first and last non-gap character in the primer alignment, and gets rid of them from the primer and the sequence.
        start = 0
        end = len(aligned_primer)
        while start < end and aligned_primer[start] == '-':
            start += 1
        while end > 0 and aligned_primer[end - 1] == '-':
            end -= 1
        aligned_primer = aligned_primer[start:end]
        aligned_sequence = alignment[1][start:end]

        record_result[sequence_key] = aligned_sequence
        
        positions = get_positions(aligned_primer)
        positions_to_check(mode, positions, aligned_primer, aligned_sequence, record_result)
    else:
        record_result[sequence_key] = "NA"
        record_result[f"{mode}_number_of_mismatches"] = "NA"
        record_result[f"{mode}_mismatches"] = "NA"

# Checks if there are any mismatches between the input primer and sequences, and the position from the 3' end (basically the inverse orientation).
def check_mismatches(mismatches, position_checks, mismatch_sets):
    for item in mismatches:
        if isinstance(item, tuple):
            position, primer_char, sequence_char = item
            if position in position_checks and (primer_char, sequence_char) in mismatch_sets:
                return True
    return False

# Checks if there are any gaps in either the primer or in the sequence. May be wise to introduce additional parameters to it in order to make it more strict.
def has_gap(mismatches):
    for item in mismatches:
        if isinstance(item, tuple):
            position, primer_char, sequence_char = item
            if primer_char == '-' or sequence_char == '-':
                return True
    return False

# Generates a diagram showing the nucleotide proportions. Stacked bars.
def nucleotide_proportions_diagram(df, mode): # Expected input here is "proportions_df".
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.fillna(0)

    df['Most_Common'] = df.idxmax(axis=1)

    plot_data = df.drop(columns='Most_Common')

    fix, ax = plt.subplots(figsize=(12, 8))

    plot_data.plot(kind='bar', stacked=True, ax=ax, colormap='coolwarm')

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Most_Common'], rotation=0, ha='center')

    plt.xlabel("Position 5' to 3'")
    plt.ylabel('Proportion')
    plt.title('Nucleotide proportions at each position')

    plt.legend(title='Nucleotide/Gaps', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()  # Automatically adjust subplot parameters to fit the figure area.

    file_name = "Reverse_primer_nucleotide_proportions.png"
    if mode == 'forward_sequence':
        file_name = "Forward_primer_nucleotide_proportions.png"
    plt.savefig(file_name)

# Generates a table showing the frequencies/occurrences of each kind of nucleotide at a given position relative to the consensus sequence and input primer.
def proportions_table(df, mode, padded_primer):
    file_name = "Reverse_primer_nucleotide_frequencies.csv"
    if mode == "forward_sequence":
        file_name = "Forward_primer_nucleotide_frequencies.csv"    
    df.to_csv(file_name, sep=";", index=True)

    # Open the CSV file and modify the top-left cell
    with open(file_name, 'r+') as file:
        lines = file.readlines()
    
        # Modify the first line, first column (top-left cell)
        lines[0] = 'Consensus_sequence_5_to_3' + ';' + lines[0][1:]
    
        # Move the cursor to the beginning of the file and write the modified content
        file.seek(0)
        file.writelines(lines)

    new_row = ';'.join(padded_primer)
    
    with open(file_name, 'r+') as file:
        content = file.read()
        file.seek(0, 0)
        file.write(new_row.rstrip('\r\n') + '\n' + content)

# Generates sequence alignments using the linsi algorithm with MAFFT.
def multiple_sequence_alignment(df, mode): # mode is either forward_sequence or reverse_sequence.
    # Convert DataFrame to FASTA format and write to a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as temp_input:
        for index, row in df.iterrows():
            if not row[mode] == "NA":
                temp_input.write(f">{row['species_name']}_{row['accession_number']}\n{row[mode]}\n")
            temp_input_path = temp_input.name
    
    # Define output file for the aligned sequences
    file_name = "Reverse_sequences_aligned.fasta"
    if mode == "forward_sequence":
        file_name = "Forward_sequences_aligned.fasta"
    output_file = file_name
    
    # Create MAFFT command line instance
    mafft_cline = MafftCommandline(
        input=temp_input_path,
        localpair=True,
        maxiterate=1000,
        thread=4,
        reorder=True
    )
    
    # Run MAFFT alignment
    stdout, stderr = mafft_cline()
    
    # Save the alignment output to the output file
    with open(output_file, "w") as handle:
        handle.write(stdout)
    
    # Clean up temporary files
    os.remove(temp_input_path)

### mode will either be "forward_sequence" or "reverse_sequence" when the function is invoked.
def omega_function(df, mode, primer): # Name it something appropriately.
    multiple_sequence_alignment(df, mode) # Creates the sequence alignments.
    df = df[mode]
    sequences = [sequence for sequence in df if sequence != "NA"]
    max_len = max(len(seq) for seq in sequences) # Gives us the longest sequence from the df. Include ones with gaps introduced.
    padded_sequences = [seq.ljust(max_len, '-') for seq in sequences] # Add gaps at the end of sequences that are shorter than max_len to ensure functionality downstreams.
    adjusted_primer = primer.ljust(max_len, '-')
    padded_primer = []
    padded_primer.append("Primer_used_5_to_3")
    for nucleotide in adjusted_primer:
        padded_primer.append(nucleotide)
    
    # Initialize a list to hold the frequency dictionaries for each position.
    frequencies = [] 

    for i in range(max_len):
        column = [seq[i] for seq in padded_sequences]
        freq = Counter(column)
        frequencies.append(freq)

    # Calculate proportions and most common nucleotide per position.
    proportions = []
    most_common_seq = []
    for freq in frequencies:
        total = sum(freq.values())
        proportion = {char: count / total for char, count in freq.items()}
        proportions.append(proportion)
        most_common_seq.append(freq.most_common()[0][0])

    # Convert proportions to a dataframe.
    proportions_df = pd.DataFrame(proportions).fillna(0) # Fill NaN with 0.
    nucleotide_order = ['A', 'C', 'G', 'T', '-']
    proportions_df = proportions_df[nucleotide_order]

    # Makes a csv file displaying the occurrences of each nucleotide at each position relative to the consenus primer.
    df_to_transpose = pd.DataFrame(frequencies).fillna(0).astype(int)
    df_to_transpose = df_to_transpose[nucleotide_order]

    df_transposed = df_to_transpose.T
    new_headers = [frequence.most_common()[0][0] for frequence in frequencies]
    df_transposed.columns = new_headers

    nucleotide_proportions_diagram(proportions_df, mode) # Creates the diagrams.
    proportions_table(df_transposed, mode, padded_primer) # Creates the tables.

###
###
###

with open(monophyletic_group, "r") as file:
    monophyletic_list = file.readlines()
monophyletic_list = [line.strip() for line in monophyletic_list]

duplicate_sequences = find_duplicates(reference_database)


# Dataframe for the species summary.
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
                        'duplicate_sequences': [],
                        'average_gc_content': [],
                        'lowest_gc_content': [],
                        'highest_gc_content': []
                       }

# Dataframe for the sequence summary.
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
                        'length': [],
                        'gc_content': []
                        }

lineage_keys = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

species_gc_content = {}
family_gc_content = {}

for record in SeqIO.parse(reference_database, "fasta"):
    split_header = record.id.split("|")
    lineage = split_header[2].split(";")
    species = lineage[6]
    family = lineage[4]
    accession = split_header[1]
    counter = int(split_header[3])
    gc_content = gc_fraction(record.seq) * 100

    if species not in species_gc_content:
        species_gc_content[species] = []
    species_gc_content[species].append(gc_content)

    if family not in family_gc_content:
        family_gc_content[family] = []
    family_gc_content[family].append(gc_content)


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
    dataframe_two_dictionary['gc_content'].append(round(gc_content, 2))

# Adds GC-content to the dataframe with species data.
for species, species_gc in species_gc_content.items():
    average_species_gc_content = round(sum(species_gc) / len(species_gc), 2)
    lowest_species_gc_content = round(min(species_gc), 2)
    highest_species_gc_content = round(max(species_gc), 2)
    dataframe_one_dictionary['average_gc_content'].append(average_species_gc_content)   
    dataframe_one_dictionary['lowest_gc_content'].append(lowest_species_gc_content)   
    dataframe_one_dictionary['highest_gc_content'].append(highest_species_gc_content)

df_one = pd.DataFrame(dataframe_one_dictionary)
df_one.to_csv(f"{database_name}_evaluation_species_summary.csv", sep=";", index=False)

df_two = pd.DataFrame(dataframe_two_dictionary)
df_two.to_csv(f"{database_name}_sequence_summary.csv", sep=";", index=False)


# Makes a histogram of the database's GC-content.
gc_contents = dataframe_two_dictionary['gc_content']
mean_gc = np.mean(gc_contents)
median_gc = np.median(gc_contents)

plt.figure(figsize=(10, 6))

plt.hist(gc_contents, bins=30)

plt.axvline(mean_gc, color='b', linestyle='dotted', linewidth=1.5, label=f'Mean: {mean_gc:.2f}%')
plt.axvline(median_gc, color='r', linestyle='dashed', linewidth=1.5, label=f'Median: {median_gc:.2f}%')

plt.xlabel('GC-content (%)')
plt.ylabel('Number of sequences')
plt.title("The reference database's GC-content distribution")
plt.grid(axis='y', linestyle='--', linewidth=0.7)

plt.legend()
plt.tight_layout()

plt.savefig(f"{database_name}_GC_content_histogram.png")

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

# The dataframe for family summary.
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

df_three['average_gc_content'] = None
df_three['lowest_gc_content'] = None
df_three['highest_gc_content'] = None


# Adds GC-content to the dataframe with family data.
for family, family_gc in family_gc_content.items():
    average_family_gc_content = round(sum(family_gc) / len(family_gc), 2)
    lowest_family_gc_content = round(min(family_gc), 2)
    highest_family_gc_content = round(max(family_gc), 2)
    df_three.loc[df_three['family'] == family, 'average_gc_content'] = average_family_gc_content
    df_three.loc[df_three['family'] == family, 'lowest_gc_content'] = lowest_family_gc_content
    df_three.loc[df_three['family'] == family, 'highest_gc_content'] = highest_family_gc_content

df_three.to_csv(f"{database_name}_evaluation_family_summary.csv", sep=";", index=False)

###
###
###

if forward_primer or reverse_primer: # This blocked is skipped entirely if neither -f nor -r is used.
    if forward_primer:
        forward_max = int(len(forward_primer) * 1.1)
    if reverse_primer:
        reverse_max = int(len(reverse_primer) * 1.1)

    primer_alignment_results = []

    for record in SeqIO.parse(reference_database, "fasta"):
        sequence = str(record.seq)
        if forward_primer:
            forward_sequence = sequence[:forward_max]
        else:
            forward_sequence = None
        if reverse_primer:
            reverse_sequence = str(Seq(sequence[-reverse_max:]).reverse_complement())
        else:
            reverse_sequence = None

        lineage = record.id.split("|")[2]
        species_name = lineage.split(";")[-1]
        accession_number = record.id.split("|")[1]
        counter = record.id.split("|")[3]

        record_result = {
            'species_lineage': lineage,
            'species_name': species_name,
            'accession_number': accession_number,
            'counter': counter
        }
        if forward_primer:
            record_result.update({
                'forward_sequence': None,
                'forward_number_of_mismatches': None,
                'forward_mismatches': None,
                'forward_sequence_status': None
            })
        if reverse_primer:
            record_result.update({
                'reverse_sequence': None,
                'reverse_number_of_mismatches': None,
                'reverse_mismatches': None,
                'reverse_sequence_status': None
            })

        if forward_primer:
            primer_alignment("forward", forward_primer, forward_sequence, record_result)
        if reverse_primer:
            primer_alignment("reverse", reverse_primer, reverse_sequence, record_result)

        primer_alignment_results.append(record_result)

    df = pd.DataFrame(primer_alignment_results)
    df.to_csv(f"{database_name}_primer_result_all_entries_info.csv", sep=";", index=False)

    # Define mismatch sets for the different rules (based on Stadhouders el al., 2010).
    purine_purine_mismatches = {('A', 'T'), ('A', 'G'), ('G', 'A'), ('G', 'C'), ('C', 'G')}
    pyrimidine_pyrmidine_mismatches = {('T', 'A'), ('T', 'C'), ('C', 'T')}
    purine_pyrimidine_mismatches = {('C', 'A'), ('A', 'C'), ('G', 'T'), ('T', 'G')}

    filtered_records_perfect = []
    filtered_records_ok = []
    filtered_records_bad = []
    filtered_records_no_data = []

    for record in primer_alignment_results:
        avoid_forward = False
        dodge_forward = False
        avoid_reverse = False
        dodge_reverse = False
        if forward_primer:
            forward_mismatch = record['forward_mismatches']
            forward_mismatches = record['forward_number_of_mismatches'] # The number of mismatches.
            forward_sequence = record['forward_sequence']
            avoid_forward = has_gap(forward_mismatch) # Check for gaps.
        if reverse_primer:
            reverse_mismatch = record['reverse_mismatches']
            reverse_mismatches = record['reverse_number_of_mismatches'] # The number of mismatches.
            reverse_sequence = record['reverse_sequence']
            avoid_reverse = has_gap(reverse_mismatch) # Check for gaps.

        all_bad_pos = range(1, 6)
        # First rule set - Standard real-time PCR (Taq DNA polymerase-based)
        if forward_primer:
            if check_mismatches(forward_mismatch, {1}, purine_purine_mismatches.union(pyrimidine_pyrmidine_mismatches)):
                avoid_forward = True
            elif check_mismatches(forward_mismatch, {2}, purine_purine_mismatches):
                avoid_forward = True
        if reverse_primer:
            if check_mismatches(reverse_mismatch, {1}, purine_purine_mismatches.union(pyrimidine_pyrmidine_mismatches)):
                avoid_reverse = True
            elif check_mismatches(reverse_mismatch, {2}, purine_purine_mismatches):
                avoid_reverse = True
        
        # Second rule set - Real-time RT-PCR using specific reverse primer (Taq DNA polymerase)
        if forward_primer:
            if check_mismatches(forward_mismatch, all_bad_pos, purine_purine_mismatches):
                avoid_forward = True
            elif check_mismatches(forward_mismatch, {1, 2}, pyrimidine_pyrmidine_mismatches):
                avoid_forward = True
            elif check_mismatches(forward_mismatch, {1}, purine_pyrimidine_mismatches):
                avoid_forward = True

        # Third rule set - rTth DNA polymerase-based real-time PCR using specific reverse primer
        if reverse_primer:
            if check_mismatches(reverse_mismatch, all_bad_pos, purine_purine_mismatches.union(pyrimidine_pyrmidine_mismatches).union(purine_pyrimidine_mismatches)):
                avoid_reverse = True


        # Check if the number of mismatches exceeds our threshold.
        if forward_primer:
            if type(forward_mismatches) == int:
                if forward_mismatches > 2:
                    avoid_forward = True
                elif 1 <= forward_mismatches <= 2:
                    dodge_forward = True
        if reverse_primer:
            if type(reverse_mismatches) == int:
                if reverse_mismatches > 2:
                    avoid_reverse = True
                elif 1 <= reverse_mismatches <= 2:
                    dodge_reverse = True

        if forward_primer:
            forward_sequence_status = (
                "Bad" if avoid_forward else
                "Ok" if dodge_forward else
                "NA" if forward_sequence == "NA" else
                "Perfect"
            )
            record['forward_sequence_status'] = forward_sequence_status

        if reverse_primer:
            reverse_sequence_status = (
                "Bad" if avoid_reverse else
                "Ok" if dodge_reverse else
                "NA" if reverse_sequence == "NA" else
                "Perfect"
            )
            record['reverse_sequence_status'] = reverse_sequence_status
            
        if forward_primer and reverse_primer: # This is triggered when both primers have been used as inputs.
            if avoid_forward or avoid_reverse:
                filtered_records_bad.append(record)
            elif dodge_forward or dodge_reverse:
                filtered_records_ok.append(record)
            elif forward_sequence == "NA" and reverse_sequence == "NA":
                filtered_records_no_data.append(record)
            else:
                filtered_records_perfect.append(record)
        elif forward_primer or reverse_primer: # If only one primer is used, this happens instead. Ensures functionality in both scenarios.
            if avoid_forward or avoid_reverse:
                filtered_records_bad.append(record)
            elif dodge_forward or dodge_reverse:
                filtered_records_ok.append(record)
            elif forward_sequence == "NA" or reverse_sequence == "NA":
                filtered_records_no_data.append(record)
            else:
                filtered_records_perfect.append(record)
            


    df_perfect = pd.DataFrame(filtered_records_perfect)
    df_perfect.to_csv(f"{database_name}_primer_result_perfect_entries.csv", sep=";", index=False)
    df_ok = pd.DataFrame(filtered_records_ok)
    df_ok.to_csv(f"{database_name}_primer_result_ok_entries.csv", sep=";", index=False)
    df_bad = pd.DataFrame(filtered_records_bad)
    df_bad.to_csv(f"{database_name}_primer_result_bad_entries.csv", sep=";", index=False)
    df_no_data = pd.DataFrame(filtered_records_no_data)
    df_no_data.to_csv(f"{database_name}_primer_result_no_data_entries.csv", sep=";", index=False)

if forward_primer:
    omega_function(df, "forward_sequence", forward_primer)
if reverse_primer:
    omega_function(df, "reverse_sequence", reverse_primer)

###
###
###

os.chdir("../..")
print(f"\nDataframes have been created.\nThe dataframes are found in {os.path.abspath(evaluation_directory)}")
if forward_primer or reverse_primer:
    print("Primer compatibility check has been performed. Please review the relevant files.")