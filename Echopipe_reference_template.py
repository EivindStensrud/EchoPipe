#!/usr/bin/env python3

# ---------------------------------------------------------------------------
# Copyright (c) 2024, Daniel Borg, Eivind Stensrud, Alexander Eiler
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# This script is for generating a template reference database for the EchoPipe pipeline.
# Read more on the GitHub repository: https://github.com/EivindStensrud/EchoPipe/tree/main
#
# The user may either choose to download sequences based on a species list of their own, or they may provide a fasta file with sequences.
# ---------------------------------------------------------------------------


import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter("ignore", BiopythonDeprecationWarning)

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
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import sys
import glob
import random

###
### Functions
###

def find_5_end_fwd_position(aligned_seq):
# Function to find the position of the first base in the forward primer.
    position_5_end_fwd = 0
    for char in aligned_seq:
        if char != "-":
            return position_5_end_fwd
        position_5_end_fwd += 1
    return None

def find_5_end_rev_position(aligned_seq):
# Function to find the position of the last occurring base in the reverse primer.
    position_5_end_rev = len(aligned_seq)
    for char in reversed(aligned_seq):
        if char != "-":
            return position_5_end_rev
        position_5_end_rev -= 1
    return None

def calculate_batches(total_records, min_size, max_size):
    num_batches = (total_records + max_size - 1) // max_size
    base_size = total_records // num_batches
    remainder = total_records % num_batches
    
    if base_size < min_size:
        base_size = min_size
        num_batches = (total_records + base_size - 1) // base_size
        remainder = total_records % num_batches

    batch_sizes = [base_size] * num_batches
    for i in range(remainder):
        batch_sizes[i] += 1

    return batch_sizes

def append_and_print_message(log_file, msg):
    with open(log_file, "a") as file:
        file.write(msg)
    print(msg)

def get_command_string():
    command_string = ' '.join(sys.argv)
    return command_string

def number_threads():
    #number of threads
    executor = ThreadPoolExecutor()
    # report the number of worker threads chosen by default
    thread_number = executor._max_workers
    thread_number = min((thread_number / 2 ), 7) # max seven, as API allows for up to 8 searches per second.
    return thread_number

def fetch_taxid(species, email, api_key):
    """Fetches the unique NCBI TaxID for a given species name with retries."""
    Entrez.email = email
    Entrez.api_key = api_key
    
    max_retries = 3
    retry_delay = 10  # Start with 10 seconds, doubling on each failure
    
    for attempt in range(max_retries):
        try:
            # Perform the search
            handle = Entrez.esearch(db="Taxonomy", term=species.strip())
            record = Entrez.read(handle)
            handle.close()
            
            if record["IdList"]:
                return record["IdList"][0]
            else:
                return None  # No TaxID found (not an error, just no result)
                
        except HTTPError as e:
            # Retry on rate limits (429) or server errors (5xx)
            if e.code == 429 or 500 <= e.code < 600:
                if attempt < max_retries - 1:
                    print(f"\033[33mHTTP Error {e.code} for '{species}'. Retrying in {retry_delay}s...\033[0m")
                    time.sleep(retry_delay)
                    retry_delay *= 2
                    continue
            # Do not retry on 404 (Not Found) or 400 (Bad Request)
            return None
            
        except Exception as e:
            # Catch other connection/parsing errors
            if attempt < max_retries - 1:
                print(f"\033[33mError fetching '{species}': {e}. Retrying in {retry_delay}s...\033[0m")
                time.sleep(retry_delay)
                retry_delay *= 2
                continue
            return None

    return None

def fetch_species_data(species, email, api_key, length_threshold, max_length, retmax, custom_query, retry_delay=10, max_retries=4):
    attempt = 0
    while attempt < max_retries:
        try:
            Entrez.email = email
            Entrez.api_key = api_key
            
            search_term = f"{species.strip()}[Organism] AND (\"{length_threshold}\"[SLEN] : \"{max_length}\"[SLEN]) AND biomol_genomic[PROP] NOT \"unverified\" {custom_query}"
            
            search_handle = Entrez.esearch(
                db="nucleotide",
                term=search_term,
                retmode="xml",
                retmax=retmax,
                usehistory="y"
            )
            search_record = Entrez.read(search_handle)
            search_handle.close()

            if int(search_record["Count"]) == 0:
                print(f"\033[31mThe search term: {search_term} did not return any results.\033[0m", flush=True)
                return None
            
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

            print(f"\033[32mResults found for: {search_term}\033[0m", flush=True)
            return fetch_record

        except HTTPError as e:
            error_message = e.read().decode()
            print(f"An HTTP error occurred:\n\033[31m{error_message}\033[0m", flush=True)
        except Exception as e:
            error_message = str(e)
            print(f"An error occurred:\n\033[31m{error_message}\033[0m", flush=True)

        attempt += 1
        print(f"Retrying {species.strip()} in {retry_delay} seconds...", flush=True)
        time.sleep(retry_delay)

    print(f"\033[31mMax retries exceeded for {species.strip()}. Moving to next species.\033[0m", flush=True)
    return None

def download_species_data(species_list, email, api_key, length_threshold, max_length, retmax, custom_query, thread_number):
    list_of_not_found = []
    with open("preformated_sequences.fasta", "w") as file, ThreadPoolExecutor(max_workers=thread_number) as executor:
        futures = {}

        #progress bar
        with tqdm(total=len(species_list), desc="Downloading species data") as pbar:
            for species in species_list:
                time.sleep(random.uniform(0.1, 0.4)) 
                
                future = executor.submit(
                    fetch_species_data, species, email, api_key, length_threshold,
                    max_length, retmax, custom_query
                )
                futures[future] = species
                pbar.update(1)

            for future in as_completed(futures):
                fetched_data = future.result()
                if fetched_data:
                    file.write(fetched_data)
                else:
                    list_of_not_found.append(futures[future])
            


def mafft_alignment(input_file, output_file):
    try:
        mafft_cline = MafftCommandline(input=input_file, auto = True)
        stdout, stderr = mafft_cline()
        
        with open(output_file, "w") as file:
            file.write(stdout)
        
        return f"Alignment completed for {input_file}"
    except Exception as e:
        return f"Error in alignment for {input_file}: {str(e)}"

def parallel_mafft(files_to_align, output_prefix, thread_number):
    results = []
    num_workers = thread_number

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(mafft_alignment, input_file, f"{output_prefix}_{i}.fasta"): input_file for i, input_file in enumerate(files_to_align)}
        
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(result)

    return results

def run_mafft_parallel(files_to_align, length_threshold, longest_amplicon_size, aligned_mafft, thread_number):
    all_filtered_sequences = []
    all_not_accepted_sequences = []
    num_workers = thread_number


    print("Running MAFFT alignment.")
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Map the process_file function to the list of files
        results = executor.map(process_file, files_to_align, [length_threshold]*len(files_to_align), [longest_amplicon_size]*len(files_to_align))

        for filtered_sequences, not_accepted_sequences in results:
            all_filtered_sequences.extend(filtered_sequences)
            all_not_accepted_sequences.extend(not_accepted_sequences)

    # Write all the filtered sequences to the file
    with open(aligned_mafft, "w") as file:
        for description, sequence in all_filtered_sequences:
            file.write(f">{description}\n")
            file.write(sequence + "\n")

    return all_not_accepted_sequences

def process_file(aligned_file, length_threshold, longest_amplicon_size):
    not_accepted_sequences = []
    mafft_cline = MafftCommandline(input=aligned_file, auto=True)
    stdout, stderr = mafft_cline()

    sequence_lengths = []
    second_fasta_string = StringIO(stdout)
    second_alignment = AlignIO.read(second_fasta_string, "fasta")

    second_start_pos_forward = None
    second_end_pos_reverse = None

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

    filtered_sequences = []
    second_fasta_string.seek(0)
    for second_seq_record in SeqIO.parse(second_fasta_string, "fasta"):
        fixed_sequence_length = len(second_seq_record.seq[second_start_pos_forward:second_end_pos_reverse].replace("-", ""))
        if second_seq_record.id == "Forward_primer" or second_seq_record.id == "Reverse_primer":
            filtered_sequences.append((second_seq_record.description, str(second_seq_record.seq.replace("-", ""))))
        elif fixed_sequence_length > length_threshold and fixed_sequence_length <= longest_amplicon_size * median_length:
            filtered_sequences.append((second_seq_record.description, str(second_seq_record.seq[second_start_pos_forward:second_end_pos_reverse].replace("-", ""))))
        else:
            not_accepted_sequences.append(second_seq_record.description)

    return filtered_sequences, not_accepted_sequences

def main():
    parser = argparse.ArgumentParser(
        prog="EchoPipe - Reference template database creator",
        description="Use this script to create a reference template for downstream analysis in the EchoPipe workflow.",
        epilog="Version 1.0")

    group_1 = parser.add_argument_group("Arguments used to generate an uncurated reference template.\nThe required arguments are input_file, -f (--forward), -r (--reverse), -e (--email) and -a (--api_key)")
    group_1.add_argument('input_file', type=str, nargs='?',
        help="Txt or CSV file species names or a fasta file that is to be converted into the reference template database (-p is then required).")
    group_1.add_argument('-f', '--forward', type=str, default="",
        help="The forward primer used to find region of interest, (5'-3' orientation).")
    group_1.add_argument('-r', '--reverse', type=str, default="",
        help="The reverse primer used to find region of interest, (5'-3' orientation).")
    group_1.add_argument('-e', '--email', type=str, default="email@email.com",
        help="Your email if NCBI needs to contact you.")
    group_1.add_argument('-a', '--api_key', type=str, default="api_key",
        help="The user's NCBI API key, allows for faster downloads.")
    group_1.add_argument('-q', '--query',
        help='The search result will include the user input search term(s). Example, limit the search to 12s region: -q "AND 12s" followed by the search term. To exclude a term write  "NOT 12s".')
    group_1.add_argument('-t', '--threshold', type=int, default=150,
        help="The minimum length of a sequence, including the primer regions. Any sequence shorter than this is discarded. Default cutoff is set to 150 bases.")
    group_1.add_argument('-l', '--length', type=int, default=22000,
        help="The longest allowed sequence length for template creation Default is set to 20 000. WARNING: The longer the sequence the more computational power is required to align the sequences.")
    group_1.add_argument('-m', '--max', type=int, default=1,
        help="The number of sequences that are downloaded per species. Increasing the number may increase coverage while increasing computational power. Recommended total should not exceed 500. Default = 1.")
    group_1.add_argument('-p', '--provided_sequences', action='store_true', default='',
        help="Use if a fasta file is provided to be used as a reference template.")
    group_1.add_argument('-z', '--longest_amplicon_size', type=float, default=2,
        help="Advanced setting: multiplier for the median length of downloaded sequences to discard abnormally long sequences after trimming. Default is set to 2*median length.")

    group_2 = parser.add_argument_group("Argument used to finish the curated reference template database.")
    group_2.add_argument('-C', '--Complete', action="store_true",
        help="Completes the reference template database.")
    group_2.add_argument('input_file_species', type=str, nargs='?',
        help="Txt or CSV file species names or a fasta file that is to be converted into the reference template database (-p is then required).")
    parser.add_argument("-n", "--random_subset", type=int, help="Number of random species to use for template creation (e.g., 50).")
    parser.add_argument("-sf", "--subset_file", type=str, help="Path to a file containing a specific subset of species to use.")

    args = parser.parse_args()

    program_start = time.time()
    local_time_startime = time.localtime(program_start)
    formatted_time_startime = time.strftime("%Y-%m-%d %H:%M:%S", local_time_startime)

    if args.input_file:
        input_file = os.path.basename(args.input_file)
    elif args.input_file_species:
        input_file = os.path.basename(args.input_file_species)
    else:
        print("Error: Both ", input_file, " and ", input_file_species," are missing. Please provide at least.")
        sys.exit(1)
    if args.subset_file:
        args.subset_file = os.path.abspath(args.subset_file)

    working_directory = "Reference_template_creation/"
    to_be_curated = "aligned_sequences_to_curate.fasta"
    species_list_name = args.input_file
    species_list_name_C = args.input_file_species
    thread_number = number_threads()
    command_string = get_command_string()
    output_prefix = "temp_split_file"
    forward_primer = args.forward.lower()
    forward_primer = forward_primer.replace('i', 'n')
    reverse_primer = reverse_complement(args.reverse).lower()
    reverse_primer = reverse_primer.replace('i', 'n')


    if not args.Complete:
        if not all ([args.input_file, args.forward, args.reverse, args.email, args.api_key]):
            parser.error("An input file with a list of species along with the following are required:\nForward primer (-f, --forward)\nReverse primer (-r, --reverse)\nEmail address (-e, --email)\nAPI key (-a, --api_key)\n\nThe flag to finish the database (-C, --Complete) requires no other inputs.\n")
        else:
            input_file = os.path.abspath(args.input_file)
            Entrez.email = args.email
            Entrez.api_key = args.api_key
            custom_query = f"{args.query}" if args.query else " "
            length_threshold = args.threshold
            longest_amplicon_size = args.longest_amplicon_size
            max_length = args.length
            retmax = args.max
            skip_download = args.provided_sequences

            try:
                os.makedirs(working_directory, exist_ok=True)
                os.chdir(working_directory)

                logs = "../Log_files/"
                os.makedirs(logs, exist_ok=True)
                DATE = datetime.today().strftime("%Y-%m-%d")
                i = 1
                run_name = f"{DATE}_{i}"
                log_file = f"{logs}{run_name}_log.txt"
                append_and_print_message(log_file,f"\n\nThe command used to run the script was: python {command_string}\n")
                print(f"##############################################################################\n")


                with open(log_file, "w") as file:
                    file.write(
                        "-----------------------------------------------------------------------------------------------------------------------------\n"
                        "Thanks for using EchoPipe, an iterative pipeline to create, curate and evaluate your reference database for environmental DNA studies\n"
                        "For more information see the GitHub repository: https://github.com/EivindStensrud/EchoPipe/tree/main\n\n"
                        f"The command used to run the script was: python {command_string}\n \n"
                        "Reference_template_creation was conducted with the following settings:\n"
                        f"Initated at: {formatted_time_startime}\n"
                        f"Input file: {args.input_file}\n"
                        f"Forward primer: 5'-{args.forward}-3'\n"
                        f"Reverse primer: 5'-{args.reverse}-3'\n"
                        f"Minimum sequence length: {args.threshold}\n"
                        f"Maximum sequence length: {args.length}\n"
                        f"Maximum amplicon cutoff multiplier: {args.longest_amplicon_size}\n")
                    if not skip_download:
                        file.write(
                            f"Email address: {args.email}\n"
                            f"API key: {args.api_key}\n"
                            f"The amount of sequences downloaded per species: {args.max}\n\n"
                            f"Custom query: {args.query}\n"
                            f'Search term: [Organism] AND ("{length_threshold}"[SLEN] : "{max_length}"[SLEN]) AND biomol_genomic[PROP] NOT "unverified" {args.query}')
                
                # Existing unique name handling
                raw_names = set()
                with open(input_file, "r") as file:
                    for line in file:
                        if line.strip():
                            raw_names.add(line.strip())

                # Define lists to track results
                filtered_species_list = []
                duplicates_found = []
                seen_taxids = set()

                with ThreadPoolExecutor(max_workers=thread_number) as executor:
                    # Submit all tasks to the pool
                    future_to_name = {
                        executor.submit(fetch_taxid, name, args.email, args.api_key): name 
                        for name in raw_names
                    }
                    
                    time.sleep(random.uniform(0.1, 0.4)) 


                    # Process results as they complete
                    for future in tqdm(as_completed(future_to_name), total=len(raw_names), desc="Filtering TaxIDs"):
                        name = future_to_name[future]
                        try:
                            taxid = future.result()
                            
                            if taxid:
                                if taxid not in seen_taxids:
                                    seen_taxids.add(taxid)
                                    filtered_species_list.append(name)
                                else:
                                    msg = f"Skipping '{name}': Duplicate TaxID ({taxid}) already represented.\n"
                                    append_and_print_message(log_file, msg)
                                    duplicates_found.append(f"{name};{taxid}")
                            else:
                                msg = f"Warning: Could not resolve TaxID for '{name}'. Skipping.\n"
                                append_and_print_message(log_file, msg)
                                
                        except Exception as exc:
                            # Catch any unexpected errors from the thread itself
                            msg = f"Generated an exception for '{name}': {exc}\n"
                            append_and_print_message(log_file, msg)

                # Save a dedicated file for duplicate records
                if duplicates_found:
                    dup_file_name = f"duplicate_taxid_entries.txt"
                    dup_file_path =f"../{dup_file_name}"
                    with open(dup_file_path, "w") as dup_file:
                        dup_file.write("Species;TaxID\n")
                        for entry in duplicates_found:
                            dup_file.write(entry + "\n")
                    print(f"Duplicate records saved to: {dup_file_name}")

                base_name = os.path.basename(input_file)
                clean_output_file = f"unique_{base_name}"
                clean_output_file_path =f"../{clean_output_file}"
                
                with open(clean_output_file_path, "w") as file:
                    for name in filtered_species_list:
                        file.write(name + "\n")
                        
                print(f"\n\033[32mSuccess! Cleaned species list saved to: {clean_output_file}\033[0m")
                print(f"Please use THIS file for the next step.\n")

                species_for_template = filtered_species_list # Default = Use all

                if args.subset_file:
                    subset_names = set()
                    if os.path.exists(args.subset_file):
                        with open(args.subset_file, "r") as f:
                            for line in f:
                                if line.strip():
                                    subset_names.add(line.strip())
                        # Only use names that are in BOTH the subset file and our valid list
                        species_for_template = [s for s in filtered_species_list if s in subset_names]
                        print(f"\nSubset Mode: Using {len(species_for_template)} species from provided list for template.")
                    else:
                        print(f"\n\033[31mError: Subset file '{args.subset_file}' not found. Using full list.\033[0m")

                elif args.random_subset is not None:
                    # Check if requested number is valid
                    if 0 < args.random_subset < len(filtered_species_list):
                        species_for_template = random.sample(filtered_species_list, args.random_subset)
                        print(f"\nSubset Mode: Randomly selected {len(species_for_template)} species for template creation.")
                        
                        # Save the random selection to a file
                        base_name = os.path.basename(input_file).replace(".csv", "").replace(".txt", "")
                        subset_filename = f"subset_{args.random_subset}_{base_name}.txt"
                        subset_filename_path = f"../{subset_filename}"
                        
                        with open(subset_filename_path, "w") as f:
                            for s in species_for_template:
                                f.write(s + "\n")
                        print(f"List of selected species saved to: {subset_filename}")
                        
                    else:
                        print(f"\nSubset Mode: Requested number ({args.random_subset}) >= total species. Using full list.")


                ncbi_output_file = "preformated_sequences.fasta"

                # Use the SUBSET list for downloading
                if not skip_download:
                    download_species_data(species_for_template, args.email, args.api_key, 
                                          length_threshold, max_length, retmax, 
                                          custom_query, thread_number)
                mafft_input_file = input_file if skip_download else ncbi_output_file
                min_batch_size = 45
                max_batch_size = 55
                records = list(SeqIO.parse(mafft_input_file, "fasta"))
                total_records = len(records)
                batch_sizes = calculate_batches(total_records, min_batch_size, max_batch_size)
                num_batches = len(batch_sizes)

                start_idx = 0
                files_to_align = []
                for i, batch_size in enumerate(batch_sizes):
                    end_idx = start_idx + batch_size
                    output_file = f"{output_prefix}_{i + 1}.fasta"
                    output_file = os.path.abspath(output_file)
                    files_to_align.append(output_file)
                    with open(output_file, "w") as file:
                        SeqIO.write(records[start_idx:end_idx], file, "fasta")
                        file.write(f">Forward_primer\n{forward_primer}\n")
                        file.write(f">Reverse_primer\n{reverse_primer}\n")
                    start_idx = end_idx

                # Maximum number of attempts
                max_attempts_mafft = 5
                attempt_mafft = 0
                thread_number_mafft = thread_number

                while attempt_mafft < max_attempts_mafft:
                    try:
                        # Attempt to run MAFFT in parallel
                        not_accepted_sequences = run_mafft_parallel(
                            files_to_align, 
                            length_threshold, 
                            longest_amplicon_size, 
                            "filtered_aligned_sequences.fasta", 
                            thread_number_mafft
                        )
                        
                        # If successful, break out of the loop
                        break

                    except Exception as e:
                        print(f"Attempt {attempt_mafft + 1} failed with thread number {thread_number_mafft}: {e}")
                        
                        # Reduce thread number if it has not reduced to 1 already
                        if thread_number_mafft > 1:
                            thread_number_mafft -= 1  # Decrease by two
                            print(f"Reducing thread number to {thread_number_mafft} and will try again.")
                        else:
                            print("Minimum thread number reached. Cannot reduce further.")
                            break  # Exit the loop if no more threads can be reduced

                    attempt_mafft += 1  # Increment the attempt_mafft counter

                if attempt_mafft == max_attempts_mafft:
                    try:
                        thread_one = int(1)
                        not_accepted_sequences = run_mafft_parallel(
                            files_to_align, 
                            length_threshold, 
                            longest_amplicon_size, 
                            "filtered_aligned_sequences.fasta", 
                            thread_one
                        )
                        print("Max attempts reached. MAFFT alignment failed, retries with half batch size.")
                    except Exception as e:
                        print("Max attempts reached. MAFFT alignment failed, failed the retry with half batch size.")

                else:
                    print("MAFFT alignment completed.")

            
                print("Running MAFFT on sequences within the marker region.")
                with open(to_be_curated, "w") as file:
                    mafft_cline = MafftCommandline(
                        input="filtered_aligned_sequences.fasta",
                        localpair=True,
                        maxiterate=10,
                        reorder=True,
                        thread=int(thread_number)
                    )
                    stdout, stderr = mafft_cline()
                    file.write(stdout)

                if not_accepted_sequences:
                    with open("non_approved_sequences.txt", "w") as file:
                        for line in not_accepted_sequences:
                            file.writelines(line + "\n")

                program_end = time.time()
                program_duration = round(program_end - program_start, 2)

                local_time_endtime = time.localtime(program_end)
                formatted_time_endtime = time.strftime("%Y-%m-%d %H:%M:%S", local_time_endtime)

                append_and_print_message(log_file,
                f"\nA draft of the reference template database has been created.\n"
                f"Make sure to review {working_directory}{to_be_curated} before finalizing it to a reference template database with -C (--Complete).\n"
                f"Finished at: {formatted_time_endtime}.\n\n"
                f"It took {program_duration} seconds to run the program.\n\n"
                "##############################################################################\n")

            finally:
                
                os.chdir("..")

    else:
        logs = "Log_files/"
        os.makedirs(logs, exist_ok=True)
        DATE = datetime.today().strftime("%Y-%m-%d")
        i = 1
        run_name = f"{DATE}_{i}"
        log_file = f"{logs}{run_name}_log.txt"

        append_and_print_message(log_file,f"\n\nThe command used to run the script was: python {command_string}\n")
        print(f"##############################################################################\n")

        with open("reference_template_database.fasta", "w") as file:
            for record in SeqIO.parse(working_directory + to_be_curated, "fasta"):
                if not "Forward_primer" in record.id and not "Reverse_primer" in record.id:
                    file.write(f">{record.description}\n")
                    file.write(f'{str(record.seq.replace("-", ""))}\n')

        delete_temp_files = glob.glob(f"Reference_template_creation/{output_prefix}*")  # Finds all files starting with 'example'
        for file in delete_temp_files:
            try:
                os.remove(file)
            except Exception:
                pass


        program_end = time.time()
        program_duration = round(program_end - program_start, 2)

        local_time_endtime = time.localtime(program_end)
        formatted_time_endtime = time.strftime("%Y-%m-%d %H:%M:%S", local_time_endtime)

        append_and_print_message(log_file,
            f"\nReference_template_creation was completed by using the -C option:\n"
            f"Initiated at: {formatted_time_startime}\n"
            f"Finished at: {formatted_time_endtime}\n\n"
            f"It took {program_duration} seconds to run the program.\n\n"
            "##############################################################################\n")
        print("reference_template.fasta has been created and is ready for use with Echopipe_database_creation.py.\n")

    if vars(args).get("Complete"):
        print("Recommendation:\n")
        print(f"Next code line:\n\033[32mpython Echopipe_database_creation.py {args.input_file} reference_template_database.fasta -e {args.email} -a {args.api_key} \033[0m \n")
       
    else:
        print("Recommendation:\n")
        print("First remove sequences arising from other gene regions; see tutorial for detailed explanation: https://github.com/EivindStensrud/EchoPipe/tree/main\n")
        print(f"Next code line:\n\033[32mpython Echopipe_reference_template.py {clean_output_file} -e {args.email} -a {args.api_key} -C \033[0m \n")


if __name__ == "__main__":
    main()
