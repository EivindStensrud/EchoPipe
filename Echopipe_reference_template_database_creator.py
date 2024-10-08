#!/usr/bin/env python3

# ---------------------------------------------------------------------------
# This script is for generating a template reference database.
#
# The user may either choose to download sequences based on a species list of their own, or they may provide a fasta file with sequences.
# ---------------------------------------------------------------------------

import warnings
from Bio import BiopythonDeprecationWarning
# Suppress the specific Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import reverse_complement
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from tqdm import tqdm
from urllib.error import HTTPError
import time
import os
import argparse
import statistics

###
###
###

parser = argparse.ArgumentParser(
    prog="EchoPipe - Reference template database creator",
    description="Use this script to create a reference template for downstream analysis in the EchoPipe workflow.",
    epilog="Version 1.0")

group_1 = parser.add_argument_group("Arguments used to generate an uncurated reference template.\nThe required arguments are input_file, -f (--forward), -r (--reverse), -e (--email) and -a (--api_key)")
group_1.add_argument('input_file', type=str, nargs='?',
    help="Txt or CSV file species names or a fasta file that is to be converted into the reference template database (-p is then required).")
group_1.add_argument('-f', '--forward',
    help="The forward primer used to find region of interest, (5'-3' orientation).")
group_1.add_argument('-r', '--reverse',
    help="The reverse primer used to find region of interest, (5'-3' orientation).")
group_1.add_argument('-e', '--email', type=str, default="eivisten@uio.no",
    help="Your email if NCBI needs to contact you.")
group_1.add_argument('-a', '--api_key', type=str, default="5538bd02e0d704e3416263ad4d51c74b7608",
    help="The user's NCBI API key, allows for faster downloads.")
group_1.add_argument('-q', '--query',
    help='The search result will include the user input search term(s). Example, limit the search to 12s region: -q "AND 12s" followed by the search term. To exclude a term write  "NOT 12s".')
group_1.add_argument('-t', '--threshold', type=int, default=150,
    help="The minimum length of a sequence, including the primer regions. Any sequence shorter than this is discarded. Default cutoff is set to 150 bases.")
group_1.add_argument('-l', '--length', type=int, default=22000,
    help="The longest allowed sequence length for template creation Default is set to 20 000. WARNING: The longer the sequence the more computational power is required to align the sequences.")
group_1.add_argument('-m', '--max', type=int, default=1,
    help="The number of sequences that are downloaded per species. Incresing the number may increase the coverage, while increasing the computational power. The total number of downloaded sequences is recommended to not exceed 500. Default = 1.")
group_1.add_argument('-p', '--provided_sequences', action='store_true', default="",
    help="Use if a fasta file is provided to be used as a reference template.")
group_1.add_argument('-z', '--longest_amplicon_size', type=float, default=2,
    help="Advanced setting. A value that is multiplied with the median length of the downloaded and trimmed sequences in order to discard sequences that are abnormally long despite the trimming process. Can be used if a more narrow range of sequence size is expected, or if it is expected to be much broader than the median. Default is set to 2*median length.")

group_2 = parser.add_argument_group("Argument used to finish the curated reference template database.")
group_2.add_argument('-C', '--Complete', action="store_true",
    help="Completes the reference template database.")
group_2.add_argument('input_file_species', type=str, nargs='?',
    help="Txt or CSV file species names or a fasta file that is to be converted into the reference template database (-p is then required).")


args = parser.parse_args()

working_directory = "Reference_template_creation/"
to_be_curated = "aligned_sequences_to_curate.fasta"
species_list_name = args.input_file
species_list_name_C = args.input_file_species


if not args.Complete:
    if not all ([args.input_file, args.forward, args.reverse, args.email, args.api_key]):
        parser.error("An input file with a list of species along with the following are required:\nForward primer (-f, --forward)\nReverse primer (-r, --reverse)\nEmail address (-e, --email)\nAPI key (-a, --api_key)\n\nThe flag to finish the database (-C, --Complete) require no other inputs.\n")
    else:
        input_file = os.path.abspath(args.input_file)
        forward_primer = args.forward.lower()
        reverse_primer = reverse_complement(args.reverse).lower()
        Entrez.email = args.email
        Entrez.api_key = args.api_key
        custom_query = f'{args.query}' if args.query else " "
        length_threshold = args.threshold
        longest_amplicon_size = args.longest_amplicon_size
        max_length = args.length
        retmax = args.max
        skip_download = args.provided_sequences


        ###
        ###
        ###

        try:
            os.makedirs(working_directory, exist_ok=True) # If the directory already exists then exist_ok=True ensures no error is raised and nothing new happens.
            os.chdir(working_directory)

            log_file = "reference_template_database_inputs.txt"
            with open(log_file, "w") as file:
                file.write(
                    "The script was run with the following settings:\n"
                    f"Input file: {args.input_file}\n"
                    f"Forward primer: 5'-{args.forward}-3'\n"
                    f"Reverse primer: 5'-{args.reverse}-3'\n"
                    f"Minimum sequence length: {args.threshold}\n"
                    f'Maximum amplicon cut off multiplier: {args.longest_amplicon_size}\n')
                if not skip_download:
                    file.write(
                        f"Email address: {args.email}\n"
                        f"API key: {args.api_key}\n"
                        f"The amount of sequences downloaded per species: {args.max}\n\n"
                        f"Custom query: {args.query}\n"
                        f'Search term: [Organism] AND ("{length_threshold}"[SLEN] : "{max_length}"[SLEN]) AND biomol_genomic[PROP] NOT "unverified" {args.query}')
                

            program_start = time.time()
            max_retries = 4

            ncbi_output_file = "preformated_sequences.fasta"

            if not skip_download:
                unique_names = set()
                with open(input_file, "r") as file:
                    for line in file:
                        unique_names.add(line.strip())
                with open(input_file, "w") as file:
                    for line in sorted(unique_names):
                        file.write(line + '\n')

                with open(input_file, "r") as file:
                    species_list = file.readlines()

                list_of_not_found = []
                with open(ncbi_output_file, "w") as file:
                    for idx, species in  tqdm(enumerate(species_list, start=1), total=len(species_list), desc="Going through the species list"):
                        retry_delay = 10
                        
                        for _ in range(max_retries):
                            try:
                                search = str(species.strip()) + f'[Organism] AND ("{length_threshold}"[SLEN] : "{max_length}"[SLEN]) AND biomol_genomic[PROP] NOT "unverified" {custom_query}'
                            
                                search_handle = Entrez.esearch(
                                    db="nucleotide",
                                    term=search,
                                    retmode="xml",
                                    retmax=retmax,
                                    usehistory="y"
                                )
                                search_record = Entrez.read(search_handle)
                                search_handle.close()

                                if int(search_record["Count"]) == 0:
                                    print(f'The search term: "{search}" did not return any results.')
                                    list_of_not_found.append(species)
                                else:    
                                    webenv = search_record["WebEnv"]
                                    query_key = search_record["QueryKey"]
                                
                                    fetch_handle = Entrez.efetch(
                                        db="nucleotide",
                                        rettype="fasta",
                                        retmode="text",
                                        retmax=retmax,
                                        webenv=webenv,
                                        query_key=query_key
                                    )
                                    fetch_record = fetch_handle.read()
                                    fetch_handle.close()
                                    file.write(fetch_record)
                                break
                            except HTTPError as e:
                                print(f"An error occurred:\n{e}")
                                error_message = e.read().decode()
                                print(error_message + "\n")
                            except Exception as e:
                                print(f"An error occurred:\n{e}")
                                error_message = e.read().decode()
                                print(error_message + "\n")

                if list_of_not_found:
                    with open("no_search_results.txt", "w") as file:
                        for line in list_of_not_found:
                            file.writelines(line)

            ###
            ###
            ###

            # Function to make batch_sizes dynamic, so that the size may be a bit more evenly distributed.
            def calculate_batches(total_records, min_size, max_size):
                # Calculate the minimum number of batches with the maximum batch size.
                num_batches = (total_records + max_size - 1) // max_size
                
                # Try to distribute the records as evenly as possible.
                base_size = total_records // num_batches
                remainder = total_records % num_batches # % = modulus.
                
                # Adjust the base_size and distribute the remainder.
                if base_size < min_size:
                    base_size = min_size
                    num_batches = (total_records + base_size - 1) // base_size
                    remainder = total_records % num_batches

                batch_sizes = [base_size] * num_batches
                for i in range(remainder):
                    batch_sizes[i] += 1

                return batch_sizes


            # Function to find the position of the first base in the forward primer.
            def find_5_end_fwd_position(aligned_seq):
                position_5_end_fwd = 0
                for char in aligned_seq:
                    if char != "-":
                        return position_5_end_fwd
                    position_5_end_fwd += 1
                return None


            # Function to find the position of the last occurring base in the reverse primer.
            def find_5_end_rev_position(aligned_seq):
                position_5_end_rev = len(aligned_seq)
                for char in reversed(aligned_seq):
                    if char != "-":
                        return position_5_end_rev
                    position_5_end_rev -= 1
                return None

            ###
            ###
            ###

            mafft_input_file = input_file if skip_download else ncbi_output_file
            output_prefix = "temp_split_file"

            min_batch_size = 35
            max_batch_size = 45

            records = list(SeqIO.parse(mafft_input_file, "fasta"))
            total_records = len(records)
            batch_sizes = calculate_batches(total_records, min_batch_size, max_batch_size)
            num_batches = len(batch_sizes)

            start_idx = 0
            for i, batch_size in enumerate(batch_sizes):
                end_idx = start_idx + batch_size
                output_file = f"{output_prefix}_{i + 1}.fasta"
                with open(output_file, "w") as file:
                    SeqIO.write(records[start_idx:end_idx], file, "fasta")
                    file.write(f">Forward_primer\n{forward_primer}\n")
                    file.write(f">Reverse_primer\n{reverse_primer}\n")
                start_idx = end_idx


            input_prefix = "temp_split_file_"
            maffted_file = "aligned_sequences.fasta"
            second_maffted_file = "filtered_aligned_sequences.fasta"

            not_accepted_sequences = []

            print("MAFFT under way. It may take some time depending on how long each sequence is.")


            mafft_max_retries = 2

            with open(maffted_file, "w") as file:
                for i in range(num_batches):
                    memory_saving_mode = False
                    mafft_retry_count = 0

                    input_file = f"{input_prefix}{i + 1}.fasta"
                    print(f"MAFFT is being performed on batch number {i + 1} out of {num_batches}.")

                    while True:
                        try:
                            mafft_cline = MafftCommandline(
                                input=input_file,
                                memsave=True if memory_saving_mode else False,
                            )
                            stdout, stderr = mafft_cline()

                            trimmed_and_filtered_sequences = []
                            fasta_string = StringIO(stdout)

                            alignment = AlignIO.read(fasta_string, "fasta")
                            for record in alignment:
                                aligned_seq = str(record.seq)
                                if record.id == "Forward_primer":
                                    start_pos_forward = find_5_end_fwd_position(aligned_seq)
                                elif record.id == "Reverse_primer":
                                    end_pos_reverse = find_5_end_rev_position(aligned_seq)

                            fasta_string.seek(0) # The string must be "reset" in order for it to be usable again.
                            for seq_record in SeqIO.parse(fasta_string, "fasta"):
                                fixed_sequence_length = len(seq_record.seq[start_pos_forward:end_pos_reverse].replace("-", ""))
                                if seq_record.id == "Forward_primer" or seq_record.id == "Reverse_primer":
                                    trimmed_and_filtered_sequences.append(f">{str(seq_record.description)}\n")
                                    trimmed_and_filtered_sequences.append(str(seq_record.seq.replace("-", "")) + "\n")
                                elif fixed_sequence_length > length_threshold:
                                    trimmed_and_filtered_sequences.append(f">{str(seq_record.description)}\n")
                                    trimmed_and_filtered_sequences.append(str(seq_record.seq[start_pos_forward:end_pos_reverse].replace("-", "")) + "\n")
                                else:
                                    #print(f"Sequence not approved. Consider manually looking up:\n{str(seq_record.description)}")
                                    not_accepted_sequences.append(seq_record.description)

                            for line in trimmed_and_filtered_sequences:
                                file.writelines(line)
                            os.remove(input_file)
                            break
                        except Exception as e:
                            print(f"An error occurred while processing batch {i + 1}: {e}")
                            if mafft_retry_count < mafft_max_retries:
                                mafft_retry_count += 1
                                memory_saving_mode = True
                                print("Retrying with memory-saving mode...")
                            else:
                                print(f"Maximum retry attempts reached. Skipping batch {i + 1}.")
                                break

            with open(second_maffted_file, "w") as file:
                mafft_cline = MafftCommandline(
                    input=maffted_file
                )
                stdout, stderr = mafft_cline()

                sequence_lengths = []
                second_trimmed_and_filtered_sequences = []
                second_fasta_string = StringIO(stdout)

                second_alignment = AlignIO.read(second_fasta_string, "fasta")
                for second_record in second_alignment:
                    second_aligned_seq = str(second_record.seq)
                    if second_record.id == "Forward_primer":
                        second_start_pos_forward = find_5_end_fwd_position(second_aligned_seq)
                    elif second_record.id == "Reverse_primer":
                        second_end_pos_reverse = find_5_end_rev_position(second_aligned_seq)
                    else:
                        sequence_length = len(second_aligned_seq.replace("-", ""))
                        sequence_lengths.append(sequence_length)

                median_length = statistics.median(sequence_lengths)

                second_fasta_string.seek(0)
                for second_seq_record in SeqIO.parse(second_fasta_string, "fasta"):
                    fixed_sequence_length = len(second_seq_record.seq[second_start_pos_forward:second_end_pos_reverse].replace("-", ""))
                    if second_seq_record.id == "Forward_primer" or second_seq_record.id == "Reverse_primer":
                        file.write(f">{str(second_seq_record.description)}\n")
                        file.write(str(second_seq_record.seq.replace("-", "")) + "\n")
                    elif fixed_sequence_length > length_threshold and fixed_sequence_length <= longest_amplicon_size*median_length:
                        file.write(f">{str(second_seq_record.description)}\n")
                        file.write(str(second_seq_record.seq[second_start_pos_forward:second_end_pos_reverse].replace("-", "")) + "\n")
                    else:
                        not_accepted_sequences.append(second_seq_record.description)

            os.remove(maffted_file)

            with open(to_be_curated, "w") as file:
                mafft_cline = MafftCommandline(
                    input=second_maffted_file,
                    localpair=True,
                    maxiterate=1000,
                    reorder=True
                )
                stdout, stderr = mafft_cline()
                file.write(stdout)


            if not_accepted_sequences:
                with open("non_approved_sequences.txt", "w") as file:
                    for line in not_accepted_sequences:
                        file.writelines(line + "\n")

            print("Reference template database created.")
            print(f"Make sure to review {working_directory}{to_be_curated} before finalizing it to a reference template database with -C (--Complete).\n")
            program_duration = round(time.time() - program_start, 2)
            print(f"It took {program_duration} seconds to run the program.")

        finally: # Is performed at the end; if the program finishes or if it crashes.
            os.chdir("..")

else: # This step is triggered with -C, --Complete. It removes the primer sequences from the file and makes sure to remove any gaps.
    with open("reference_template_database.fasta", "w") as file:
        for record in SeqIO.parse(working_directory + to_be_curated, "fasta"):
            if not "Forward_primer" in record.id and not "Reverse_primer" in record.id:
                file.write(f">{record.description}\n")
                file.write(f"{str(record.seq.replace('-', ''))}\n")

    print("reference_template_database.fasta has been created and ready for use with Echopipe_database_creation.py")

if vars(args).get('Complete'):
    print("\nRecommendation:\n")
    print(f"Next code line: python Echopipe_database_creation.py {args.input_file} reference_template_database.fasta -e {args.email} -a {args.api_key} \n")

else:
    print("\nRecommendation:\n")
    print("First remove sequences arising from other gene regions, see tutorial for detailed explanation:\n")
    print(f"Next code line: python Echopipe_reference_template_database_creator.py -e {args.email} -a {args.api_key} {args.input_file} -C \n")

