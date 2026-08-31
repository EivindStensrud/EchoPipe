#!/usr/bin/env python3

# ---------------------------------------------------------------------------
# Copyright (c) 2026, Eivind Stensrud, Daniel Borg, Alexander Eiler
#
# Distributed under the terms of the Modified BSD License.
# The full license is in the file LICENSE, distributed with this software.
# ---------------------------------------------------------------------------

import argparse
import glob
import os
import random
import re
import statistics
import sys
import tempfile
import time
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from io import StringIO
from Bio.SeqUtils import gc_fraction
import json

# Warning Suppression
warnings.simplefilter("ignore", FutureWarning)
try:
    from Bio import BiopythonWarning, BiopythonDeprecationWarning
    warnings.simplefilter("ignore", BiopythonWarning)
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
except ImportError:
    pass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, Entrez, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align.Applications import MafftCommandline
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio.Phylo.Applications import FastTreeCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Seq import reverse_complement
from ete3 import PhyloTree
from tqdm import tqdm
from urllib.error import HTTPError


# =============================================================================
# SHARED UTILITIES (All Pipeline Helpers)
# =============================================================================

def save_config(email, api_key, forward="", reverse="", config_file="Log_files/config.json"):
    """Saves user credentials and primers to avoid typing them repeatedly."""
    os.makedirs(os.path.dirname(config_file), exist_ok=True)
    with open(config_file, 'w') as f:
        json.dump({
            "email": email, 
            "api_key": api_key,
            "forward": forward,
            "reverse": reverse
        }, f)
        
def load_config(config_file="Log_files/config.json"):
    """Loads user credentials and primers if they exist."""
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            return json.load(f)
    return {"email": "email@email.com", "api_key": "api_key", "forward": "", "reverse": ""}

def get_command_string():
    """Returns the exact command line string used by the user."""
    return ' '.join(sys.argv)

def append_and_print_message(log_file, msg):
    """Prints a message to the console and appends it to the unified log file."""
    if log_file:
        with open(log_file, "a") as file:
            file.write(msg)
    print(msg)

def read_counter(counter_file, date=None):
    """Reads the run counter for the current date."""
    try:
        if os.path.exists(counter_file):
            with open(counter_file, "r") as f:
                data = f.read().strip()
                counters = dict(line.split(":") for line in data.splitlines())
                if date is not None:
                    return int(counters.get(date, 0))
                elif counters:
                    latest_date = max(counters.keys(), key=lambda d: datetime.strptime(d, '%Y-%m-%d'))
                    return int(counters[latest_date])
    except Exception as e:
        print(f"Error reading the counter file: {e}")
    return 0

def update_counter(counter_file, date, counter):
    """Updates the counters"""
    os.makedirs(os.path.dirname(counter_file), exist_ok=True)
    counters = {}
    if os.path.exists(counter_file):
        with open(counter_file, "r") as f:
            counters = dict(line.split(':') for line in f.read().splitlines())
    counters[date] = str(counter)
    with open(counter_file, "w") as f:
        for dt, cnt in counters.items(): 
            f.write(f"{dt}:{cnt}\n")

def get_log_file_name(logs_dir="Log_files/"):
    """Generates a standardized log file name based on the date and counter."""
    os.makedirs(logs_dir, exist_ok=True)
    counter_file = os.path.join(logs_dir, "run_counters.txt")
    
    latest_date = None
    if os.path.exists(counter_file):
        with open(counter_file, "r") as f:
            lines = f.read().strip().splitlines()
            if lines:
                latest_date = max(line.split(':')[0] for line in lines)
    
    if latest_date is None:
        latest_date = datetime.today().strftime('%Y-%m-%d')

    counter = read_counter(counter_file, latest_date)
    if counter == 0:
        counter = 1

    run_name = f"{latest_date}_{counter}"
    log_file_path = os.path.abspath(os.path.join(logs_dir, f"{run_name}_log.txt"))
    
    return log_file_path, run_name

def log_start(log_file, task_name):
    """Prints and logs the start of a task."""
    msg = f"\n{'='*80}\nStarting {task_name} at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    append_and_print_message(log_file, msg)

def log_end(log_file, start_time, next_step_msg):
    """Prints and logs the end of a task, including the duration and next step."""
    duration = round(time.time() - start_time, 2)
    msg = (f"\nTask finished at {datetime.now().strftime('%H:%M:%S')}. "
           f"Duration: {duration} seconds.\n"
           f"{'-'*80}\n\n"
           f"Recommendation:\n{next_step_msg}\n")
    append_and_print_message(log_file, msg)

def get_accession(header_line):
    """Safely extracts the accession number from standard or aligned FASTA headers."""
    header_line = header_line.strip()
    if header_line.startswith(">gb_"):
        return header_line.split("_")[1]
    elif "|" in header_line:
        return header_line.split("|")[1]
    else:
        return header_line.split()[0].replace('>', '')

def find_duplicates(fasta_file):
    """Identifies identical sequences in a FASTA file and groups their headers."""
    sequence_dict = defaultdict(list)
    current_header = None
    current_sequence = []
    
    with open(fasta_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    cleaned_sequence = "".join(current_sequence).replace("-", "")
                    sequence_dict[cleaned_sequence].append(current_header)
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
                
        if current_header and current_sequence:
            cleaned_sequence = "".join(current_sequence).replace("-", "")
            sequence_dict[cleaned_sequence].append(current_header)
            
    return {seq: headers for seq, headers in sequence_dict.items() if len(headers) > 1}

def get_species_name(leaf, pattern=r'\.\d+_|[\d.]+'):
    name_parts = re.split(pattern, leaf.name)
    parts = [p for p in name_parts if p]
    return '_'.join(parts[1].split('_')[7:-1])

def analysis_completed(finished_accession_numbers, analysed_accession_numbers, reuse_sequences):
    """Updates the analysed accession numbers list after BLAST."""
    if reuse_sequences:
        with open(analysed_accession_numbers, "w") as file:
            file.truncate(0)

    with open(finished_accession_numbers, "r") as file:
        finished_accession_content = set(file.readlines())
    with open(analysed_accession_numbers, "r") as file:
        previous_accession_content = set(file.readlines())

    new_additions = finished_accession_content.difference(previous_accession_content)

    with open(analysed_accession_numbers, "a") as file:
        for line in new_additions:
            file.write(line)

    with open(analysed_accession_numbers, "r") as file:
        sorted_accession_numbers = file.readlines()
        sorted_accession_numbers.sort()

    with open(analysed_accession_numbers, "w") as file:
        file.writelines(sorted_accession_numbers)

def concatenate_files(input_list, concatenated_file):
    """Concatenates multiple files into one."""
    with open(concatenated_file, "w") as outfile:
        for previous_file in input_list:
            with open(previous_file, "r") as infile:
                for line in infile:
                    outfile.write(line)

def columnize(blast_output):
    """Converts BLAST TSV output to a pandas DataFrame."""
    data_frame = pd.read_table(blast_output, header=None)
    data_frame.columns = [
        "qseqid", "qseq", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    return data_frame

def write_fasta(data, blast_filename):
    """Creates a FASTA file from processed BLAST results."""
    with open(blast_filename, "w") as file_conn:
        for _, row in data.iterrows():
            header = f'>{row["qseqid"]}|{row["identical_seqs"]}\n'
            sequence = str(row["qseq"])
            file_conn.write(header + sequence + '\n')

def create_directory(name_of_directory):
    """Safely creates a directory if it does not already exist."""
    os.makedirs(name_of_directory, exist_ok=True)

def does_file_exist(file_to_check_path):
    """Creates an empty file if the path does not exist."""
    if not os.path.isfile(file_to_check_path):
        with open(file_to_check_path, "w") as file:
            pass

def update_last_line(file_path, progress_info, idx):
    """Updates the last line of a file for ongoing progress tracking."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    if lines and idx != 1:
        lines[-1] = f"{progress_info}\n"
    else:
        lines.append(f"{progress_info}\n")
        
    with open(file_path, 'w') as file:
        file.writelines(lines)

def suggest_parallelism_fix(log_file=None):
    """Logs and prints helpful troubleshooting advice for NCBI connection issues."""
    msg = (
        "\n" + "="*80 + "\n"
        "[NCBI CONNECTION & THREAD NOTICE]\n"
        "Frequent 'IncompleteRead' or network timeout errors occur when NCBI throttles\n"
        "concurrent API requests or drops heavy HTTP download payloads.\n\n"
        "To resolve this, reduce parallel load by running the command with:\n\n"
        "  1. LOWER THREAD COUNT (-T / --threads):\n"
        "     Reduce parallel worker threads (e.g., -T 2 or -T 3) to prevent rate-limiting.\n\n"
        "  2. LOWER DOWNLOAD BATCH SIZE (-b / --batch_size):\n"
        "     Pass a smaller batch size (e.g., -b 1000 or -b 2000, default: 5000).\n\n"
        "Example Command:\n"
        "  python echopipe.py create <species.txt> <template.fasta> -T 2 -b 1000\n"
        + "="*80 + "\n"
    )
    append_and_print_message(log_file, msg)


def repeated_failures(consecutive_fail_counter, log_file, check_for_fail, type_of_download):
    """Monitors consecutive errors and pauses or kills the program if thresholds are met."""
    if consecutive_fail_counter == 2:
        msg = "\nTwo consecutive download errors encountered. Pausing for 10 minutes before reattempting.\n"
        append_and_print_message(log_file, msg)
        time.sleep(600)
    elif consecutive_fail_counter == 3:
        msg = "\nErrors persist. Shutting down the program."
        if check_for_fail:
            with open(log_file, "a") as file:
                file.write(f'\n{type_of_download} finished with errors:\n')
                file.write('\n'.join(check_for_fail) + '\n')
        append_and_print_message(log_file, msg)
        
        # Display the troubleshooting guidance before exiting
        suggest_parallelism_fix(log_file)
        
        sys.exit()

def loop_finished(check_for_fail, log_file, type_of_download, loop_duration):
    """Appends summary of a finished download loop to the log file."""
    with open(log_file, "a") as file:
        if check_for_fail:
            file.write(f'\n{type_of_download} finished with notes:\n')
            file.write('\n'.join(check_for_fail) + '\n')
        else:
            file.write(f'\n{type_of_download} finished! No errors detected.\n'
                       f'It took {loop_duration} seconds to finish.\n\n'
                       '--------------------------------------------------------------------------\n\n\n')

def write_lines_from_list(list_with_lines, lines_to_file_path, log_file=None, msg=None):
    """Writes a list of lines to a file, optionally logging a completion message."""
    with open(lines_to_file_path, "w") as file:
        file.writelines(list_with_lines)
    if msg is not None:
        append_and_print_message(log_file, msg)

def error_occurrence(error, species, retry_delay):
    """Handles runtime errors, prints diagnostics, and implements exponential backoff."""
    print(f"An error occurred: {error}")
    print(f"Error occurred for {species}")
    print(f"Initiating a {retry_delay} seconds pause before retrying.")
    time.sleep(retry_delay)
    return error, species

def find_5_end_fwd_position(aligned_seq):
    """Function to find the position of the first base in the forward primer."""
    position_5_end_fwd = 0
    for char in aligned_seq:
        if char != "-":
            return position_5_end_fwd
        position_5_end_fwd += 1
    return None

def find_5_end_rev_position(aligned_seq):
    """Function to find the position of the last occurring base in the reverse primer."""
    position_5_end_rev = len(aligned_seq)
    for char in reversed(aligned_seq):
        if char != "-":
            return position_5_end_rev
        position_5_end_rev -= 1
    return None

def calculate_batches(total_records, min_size, max_size):
    """Calculates optimal batch sizes for parallel processing."""
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

def number_threads(user_threads=None):
    """
    Determines the number of threads.
    If user_threads is provided, it uses that value but caps it at 7.
    Otherwise, auto-calculates thread count capped at 7.
    """
    if user_threads is not None and user_threads > 0:
        return min(int(user_threads), 7)

    executor = ThreadPoolExecutor()
    thread_number = executor._max_workers
    # Max 7, as NCBI API allows for up to 10 requests per second.
    thread_number = min((thread_number / 2), 7) 
    return max(1, int(thread_number))

def process_species_taxonomy(species, filtered_species_collection, taxonomic_ranks, max_retries):
    """Download and process species taxonomy data from NCBI."""
    retry_delay = 10
    consecutive_fail_counter = 0

    for i in range(max_retries):
        try:
            # Performs an esearch where we are interested in the taxids from the species in the list.
            taxid_handle = Entrez.esearch( 
                db="Taxonomy", 
                term=species,
                retmode="xml",
                usehistory="y"
            )
            taxid_record = Entrez.read(taxid_handle)
            taxid_handle.close()
            taxid = taxid_record["IdList"]

            time.sleep(0.5)

            # Check if search came out negative
            if taxid_record["Count"] == str(0) or "ErrorList" in taxid_record:  
                return {'species_not_found': species}
            # Check for duplicate entry
            elif any(taxid[0] == line.split(';')[1].strip() for line in filtered_species_collection):  
                return {'duplicate_species': species}
            else:
                webenv = taxid_record["WebEnv"]
                query_key = taxid_record["QueryKey"]

                fetch_handle = Entrez.efetch(
                    db="Taxonomy",
                    retmode="xml",
                    webenv=webenv,
                    query_key=query_key
                )
                fetch_taxonomy_data = Entrez.read(fetch_handle)
                fetch_handle.close()

                scientific_name = fetch_taxonomy_data[0]["ScientificName"]
                LineageEx = fetch_taxonomy_data[0]["LineageEx"]
                taxa_rank_list = []

                for taxonomic_rank in taxonomic_ranks:
                    taxonomic_rank_exist = False 

                    for item in LineageEx: 
                        if item.get("Rank") == taxonomic_rank: 
                            taxonomic_rank_exist = True 
                            tax_name = (item["ScientificName"])
                            cleaned_taxa_name = re.match(r'([^:;(\s]+)', tax_name)
                            if cleaned_taxa_name:
                                cleaned_taxa_name = cleaned_taxa_name.group(0).strip()
                            else:
                                cleaned_taxa_name = tax_name
                                
                            taxa_rank_list.append(cleaned_taxa_name) 

                    if not taxonomic_rank_exist: 
                        taxa_rank_list.append("NA")

                # Removes space between Genus species to get Genus_species to be used in the header.
                taxa_rank_list.append(scientific_name.replace(' ', '_'))
                taxa_to_fasta = ";".join(taxa_rank_list)

                taxid_entry = f"{species}; {taxid[0]}; {scientific_name}; {taxa_to_fasta}"
                return {'filtered_species_collection': taxid_entry, 'species_found': species}

            consecutive_fail_counter = 0
            break

        except HTTPError as error:
            last_exception = error_occurrence(error, species, retry_delay)
            return {'failed_information_downloads': f"HTTPError: Failed to download info about {species} - {error}"}
        except Exception as e:
            last_exception = error_occurrence(e, species, retry_delay)
            return {'failed_information_downloads': f"Exception: Failed to download info about {species} - {e}"}
    else:
        return {'failed_information_downloads': f"Max retries reached for {species}"}


def process_species_accession_number(specie, search_settings, max_search, batch_size, max_retries, sort_by_sequence_length, consecutive_fail_counter, log_file, failed_accession_downloads, type_of_download, run_name, create_directory, repeated_failures, error_occurrence, over_cap_csv_list):
    """Download and process new accession number for the species from NCBI."""
    species_directory = f"Species/{specie}"
    create_directory(species_directory)
    accession_file_1 = f"{species_directory}/{specie}_{run_name}_accession_numbers.txt"
    
    search = str(specie.replace("_", " ") + search_settings)
    retry_delay = 10
    repeated_failures(consecutive_fail_counter, log_file, failed_accession_downloads, type_of_download)

    for i in range(max_retries):
        try:
            accession_collection = []
            for start in range(0, max_search, batch_size):
                search_handle = Entrez.esearch(
                    db="nucleotide",
                    term=search,
                    retmax=batch_size,
                    retstart=start,
                    sort="Sequence Length" if sort_by_sequence_length else None,
                    idtype="acc"
                )
                search_record = Entrez.read(search_handle)
                search_handle.close()   
                accession_collection.append(search_record["IdList"])

                count = int(search_record["Count"])
                if count > max_search:
                    over_cap_csv_list.append(f"{specie};{str(max_search)};{str(count)}\n")
                    count = max_search

                time.sleep(0.5)

                if count <= 10000:
                    break

            consecutive_fail_counter = 0
            break
                
        except HTTPError as error:
            last_exception = error_occurrence(error, specie, retry_delay)
        except Exception as error:
            last_exception = error_occurrence(error, specie, retry_delay)

    else:
        print(f"Unable to retrieve accession numbers for {specie}.")
        log_entry = f"Error: Failed to download accession numbers for {specie} - {last_exception}"
        failed_accession_downloads.append(log_entry)
        consecutive_fail_counter += 1
    
    with open(accession_file_1, "w") as file:
        for accession_batch in accession_collection:
            for accession_number in accession_batch:
                file.write(accession_number + "\n")

    return specie

def fetch_taxid(species, email, api_key):
    """
    Fetches the unique NCBI TaxID for a given species name.
    
    Implements exponential backoff for handling API rate limits (HTTP 429) 
    and server errors (HTTP 5xx).
    """
    Entrez.email = email
    Entrez.api_key = api_key
    
    max_retries = 3
    retry_delay = 10  # Starts at 10 seconds, doubles on subsequent failures
    
    for attempt in range(max_retries):
        try:
            # Perform the NCBI Taxonomy database search
            handle = Entrez.esearch(db="Taxonomy", term=species.strip())
            record = Entrez.read(handle)
            handle.close()
            
            if record.get("IdList"):
                return record["IdList"][0]
            return None  # No TaxID found for the provided species name
            
        except HTTPError as e:
            # Retry only on rate limits (429) or server-side errors (500-599)
            if e.code == 429 or 500 <= e.code < 600:
                if attempt < max_retries - 1:
                    print(f"\033[33mHTTP Error {e.code} for '{species}'. Retrying in {retry_delay}s...\033[0m")
                    time.sleep(retry_delay)
                    retry_delay *= 2
                    continue
            # Do not retry on client errors like 400 (Bad Request) or 404 (Not Found)
            return None
            
        except Exception as e:
            # Catch other connection or xml parsing errors
            if attempt < max_retries - 1:
                print(f"\033[33mError fetching '{species}': {e}. Retrying in {retry_delay}s...\033[0m")
                time.sleep(retry_delay)
                retry_delay *= 2
                continue
            return None

    return None
def fetch_fasta_in_chunks(accession_list, batch_size=1000, log_file=None):
    """Downloads FASTA sequences in smaller sub-batches to prevent NCBI IncompleteRead timeouts."""
    all_records = []
    
    for i in range(0, len(accession_list), batch_size):
        chunk = accession_list[i : i + batch_size]
        max_retries = 3
        
        for attempt in range(max_retries):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=chunk,
                    rettype="fasta",
                    retmode="text"
                )
                data = handle.read()
                handle.close()
                all_records.append(data)
                break  # Success, move to next chunk
                
            except Exception as e:
                if attempt < max_retries - 1:
                    time.sleep(10 * (attempt + 1))  # Exponential backoff: 10s, 20s...
                else:
                    # If 1000 failed, try emergency micro-batch of 250
                    print(f"Warning: Chunk download failed for batch starting at index {i}. Retrying with micro-batches...")
                    for micro_i in range(0, len(chunk), 250):
                        micro_chunk = chunk[micro_i : micro_i + 250]
                        try:
                            h = Entrez.efetch(db="nucleotide", id=micro_chunk, rettype="fasta", retmode="text")
                            all_records.append(h.read())
                            h.close()
                        except Exception:
                            pass

    return "".join(all_records)
    
def download_species_data(species_list, email, api_key, length_threshold, max_length, retmax, custom_query, thread_number):
    """
    Concurrently downloads sequence data for a list of species from NCBI.
    
    Successfully retrieved fasta sequences are written to 'preformated_sequences.fasta'.
    Returns a list of species that returned no sequences or failed to download.
    """
    list_of_not_found = []
    
    with open("preformated_sequences.fasta", "w") as file:
        with ThreadPoolExecutor(max_workers=thread_number) as executor:
            futures = {}

            # 1. Submit all download tasks to the thread pool
            for species in species_list:
                time.sleep(random.uniform(0.1, 0.4))  # Stagger requests to avoid API spikes
                
                future = executor.submit(
                    fetch_species_data, species, email, api_key, length_threshold,
                    max_length, retmax, custom_query
                )
                futures[future] = species

            # 2. Process results as they finish and update the progress bar
            with tqdm(total=len(species_list), desc="Downloading species data") as pbar:
                for future in as_completed(futures):
                    species_name = futures[future]
                    try:
                        fetched_data = future.result()
                        if fetched_data:
                            file.write(fetched_data)
                        else:
                            list_of_not_found.append(species_name)
                    except Exception as e:
                        # Catch unexpected thread failures and mark species as not found
                        print(f"\n\033[33mError processing {species_name}: {e}\033[0m")
                        list_of_not_found.append(species_name)
                    
                    # Update progress bar only after a task actually completes
                    pbar.update(1)

    return list_of_not_found

def fetch_species_data(species, email, api_key, length_threshold, max_length, retmax, custom_query, retry_delay=10, max_retries=4):
    """
    Searches and fetches nucleotide FASTA sequences for a given species from NCBI.
    
    Uses NCBI's E-utilities (esearch and efetch) with history tracking.
    Implements exponential backoff for error handling to respect API limits.
    """
    Entrez.email = email
    Entrez.api_key = api_key
    species_clean = species.strip()
    
    # Construct the highly specific eSearch term
    search_term = f'{species_clean}[Organism] AND ("{length_threshold}"[SLEN] : "{max_length}"[SLEN]) AND biomol_genomic[PROP] NOT "unverified" {custom_query}'
    
    for attempt in range(max_retries):
        try:
            # 1. Search for records and post to NCBI History server
            search_handle = Entrez.esearch(
                db="nucleotide",
                term=search_term,
                retmode="xml",
                retmax=retmax,
                usehistory="y"
            )
            search_record = Entrez.read(search_handle)
            search_handle.close()

            # If no sequences match the criteria, exit cleanly (not an error)
            if int(search_record.get("Count", 0)) == 0:
                return None
            
            # 2. Fetch the records using the WebEnv history tracking
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

            return fetch_record

        except HTTPError as e:
            # Handle NCBI server overloads or rate limits
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
                retry_delay *= 2  # Exponential backoff
                continue
        except Exception as e:
            # Handle connection drops or incomplete XML parsing
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
                retry_delay *= 2  # Exponential backoff
                continue

    # Only print to terminal if a species completely fails after all retries
    print(f"\n\033[31mMax retries exceeded for {species_clean}. Moving to next species.\033[0m", flush=True)
    return None

def run_mafft_parallel(files_to_align, length_threshold, longest_amplicon_size, aligned_mafft, thread_number):
    """
    Runs MAFFT alignment on multiple batch files in parallel.
    
    Collects sequences that pass the length filters and writes them to a final FASTA file.
    Returns a list of sequences that were rejected based on amplicon size constraints.
    """
    all_filtered_sequences = []
    all_not_accepted_sequences = []
    
    print(f"Running MAFFT alignment on {len(files_to_align)} batches...")
    
    with ThreadPoolExecutor(max_workers=thread_number) as executor:
        # Submit all alignment tasks to the executor
        futures = {
            executor.submit(process_file, file_path, length_threshold, longest_amplicon_size): file_path
            for file_path in files_to_align
        }

        # Process results as they finish and update the progress bar
        with tqdm(total=len(files_to_align), desc="Aligning sequence batches") as pbar:
            for future in as_completed(futures):
                try:
                    filtered_sequences, not_accepted_sequences = future.result()
                    all_filtered_sequences.extend(filtered_sequences)
                    all_not_accepted_sequences.extend(not_accepted_sequences)
                except Exception as e:
                    file_path = futures[future]
                    print(f"\n\033[31mError processing alignment batch {file_path}: {e}\033[0m")
                
                pbar.update(1)

    # Write all the successfully filtered sequences to the output file
    with open(aligned_mafft, "w") as file:
        for description, sequence in all_filtered_sequences:
            file.write(f">{description}\n{sequence}\n")

    return all_not_accepted_sequences

def process_file(aligned_file, length_threshold, longest_amplicon_size):
    """
    Runs MAFFT alignment on a batch file, extracts the amplicon region based on 
    primer positions, and filters sequences based on length constraints.
    
    Returns a tuple containing a list of accepted sequences and a list of rejected headers.
    """
    not_accepted_sequences = []
    filtered_sequences = []
    
    # 1. Run MAFFT alignment on the batch
    mafft_cline = MafftCommandline(input=aligned_file, auto=True)
    stdout, stderr = mafft_cline()

    # Parse the alignment output from MAFFT
    alignment_data = StringIO(stdout)
    
    # Read all records into memory to avoid parsing the string stream twice
    records = list(AlignIO.read(alignment_data, "fasta"))

    start_pos_forward = None
    end_pos_reverse = None
    sequence_lengths = []

    # 2. Determine primer boundaries and calculate raw sequence lengths
    for record in records:
        aligned_seq = str(record.seq)
        if record.id == "Forward_primer":
            start_pos_forward = find_5_end_fwd_position(aligned_seq)
        elif record.id == "Reverse_primer":
            end_pos_reverse = find_5_end_rev_position(aligned_seq)
        else:
            # Calculate the raw, unaligned length (no gaps)
            sequence_length = len(aligned_seq.replace("-", ""))
            sequence_lengths.append(sequence_length)

    # Ensure we have a median length (handle edge cases gracefully)
    median_length = statistics.median(sequence_lengths) if sequence_lengths else 0

    # 3. Extract amplicons and filter based on size constraints
    for record in records:
        if record.id in ["Forward_primer", "Reverse_primer"]:
            # Keep primers completely intact (removing alignment gaps)
            filtered_sequences.append((record.description, str(record.seq).replace("-", "")))
            continue

        # Extract the targeted amplicon region between the two primers
        extracted_region = str(record.seq)[start_pos_forward:end_pos_reverse]
        amplicon_length = len(extracted_region.replace("-", ""))

        # Constraint: must be longer than the threshold and shorter than max median multiplier
        if amplicon_length > length_threshold and amplicon_length <= (longest_amplicon_size * median_length):
            filtered_sequences.append((record.description, extracted_region.replace("-", "")))
        else:
            # Log rejected sequences for user reference
            not_accepted_sequences.append(record.description)

    return filtered_sequences, not_accepted_sequences

def download_fasta_sequences(specie, run_name, taxid_collection, sequence_batch_size, max_retries, log_file, consecutive_fail_counter, failed_information_downloads, type_of_download):
    """Download fasta sequences with adaptive batch size reduction after 1 re-attempt."""
    header_lineage = None
    with open(taxid_collection, "r") as file:
        species_info = file.readlines()
        for line in species_info:
            if specie in line:
                header_lineage = line.strip().split("; ")[3]
                break

    if header_lineage is None:
        print(f"Header lineage for {specie} not found.")
        return

    path_to_new_acc = f"Species/{specie}/{specie}_new_accession_numbers_{run_name}.txt"
    try:
        with open(path_to_new_acc, "r") as file:
            new_accession_numbers = [line.strip() for line in file.readlines()]
    except FileNotFoundError:
        print(f"\033[31mAccession numbers file not found for {specie.replace('_', ' ')}\033[0m")
        return

    retry_delay = 10
    repeated_failures(consecutive_fail_counter, log_file, failed_information_downloads, type_of_download)

    current_batch_size = sequence_batch_size

    for attempt in range(max_retries):
        try:
            with open(f"Species/{specie}/{specie}_new_seqs_{run_name}.fasta", "w") as file:
                for i in range(0, len(new_accession_numbers), current_batch_size):
                    batch = new_accession_numbers[i:i + current_batch_size]
                    with Entrez.efetch(db="nucleotide", idtype="acc", id=batch, rettype="fasta", retmode="text") as handle:
                        for seq_record in SeqIO.parse(handle, "fasta"):
                            seq_record.id = seq_record.id.replace(' ', '').replace('_','')
                            file.write(f">gb|{seq_record.id}|{header_lineage}\n{seq_record.seq}\n") 
                    time.sleep(1) 
            consecutive_fail_counter = 0
            break

        except (HTTPError, Exception) as error:
            # Check if this was the second failure (attempt >= 1) before halving
            if attempt >= 1:
                current_batch_size = max(250, current_batch_size // 2)
                print(f"\033[33mRepeated download issue for {specie.replace('_', ' ')}. Reducing batch size to {current_batch_size}...\033[0m")
            else:
                print(f"\033[33mDownload issue for {specie.replace('_', ' ')}. Retrying once more with same batch size ({current_batch_size})...\033[0m")
                
            last_exception = error_occurrence(error, specie, retry_delay)
            
    else:
        log_entry = f"Error: Failed to download information about {specie} - {last_exception}"
        failed_information_downloads.append(log_entry)
        consecutive_fail_counter += 1

def perform_blast_operations(new_accession_species, run_name, database_path, analysed, min_length, outfmt_string, evalue, reuse_sequences, log_file, specie):
    """Perform BLAST operations for each species in the list."""
    failed_blasts = []
    nothing_new_to_blast = []
    species_directory = f"Species/{specie}"

    if reuse_sequences:
        blast_again = []
        for sequence_file in os.listdir(species_directory):
            if "new_seqs" in sequence_file:
                blast_again.append(os.path.join(species_directory, sequence_file))
        file_to_blast = f"{species_directory}/{specie}_temp_concatenated_blast_file.fasta"
        concatenate_files(blast_again, file_to_blast)

        old_accession_collection = []
        for old_accession_file in os.listdir(species_directory):
            if "new_accession_numbers" in old_accession_file:
                old_accession_collection.append(os.path.join(species_directory, old_accession_file))
        finished_accession_numbers = f"{species_directory}/{specie}_temp_concatenated_accession_file.txt"
        concatenate_files(old_accession_collection, finished_accession_numbers)

    else:
        file_to_blast = f"{species_directory}/{specie}_new_seqs_{run_name}.fasta"
        finished_accession_numbers = f"{species_directory}/{specie}_{run_name}_accession_numbers.txt"

    blast_output = f"{species_directory}/{specie}_{run_name}_fasta_blast_results.txt"
    analysed_accession_numbers = f"{species_directory}/{specie}_{analysed}.txt"

    if not reuse_sequences and os.path.isfile(file_to_blast) and os.path.getsize(file_to_blast) == 0:
        analysis_completed(finished_accession_numbers, analysed_accession_numbers, reuse_sequences)
        os.remove(file_to_blast)
        return

    if reuse_sequences or os.path.exists(file_to_blast): 
        blastn_cline = NcbiblastnCommandline(
            query=file_to_blast,
            db=database_path,
            out=blast_output, 
            max_target_seqs=100,
            evalue=evalue,
            word_size=11,
            reward=2,
            penalty=-3,
            gapopen=5,
            gapextend=2,
            outfmt=outfmt_string
        )
        
        stdout, stderr = blastn_cline()

        if stderr:
            log_entry = f"Error: BLAST failed for {specie} - {stderr}"
            failed_blasts.append(log_entry)

        elif os.path.getsize(blast_output) == 0:
            analysis_completed(finished_accession_numbers, analysed_accession_numbers, reuse_sequences)
            nothing_new_to_blast.append(specie.replace('_', ' '))

        else: 
            blast_df = columnize(blast_output)
            blast_df['length_minus_gapopen'] = blast_df['length'] - blast_df['gapopen']

            blast_result = (
                blast_df
                .sort_values(by='qseqid')
                .loc[lambda x: x['length'] > min_length]
                .assign(qseq=lambda x: x['qseq'].str.replace('-', ''))
                .pipe(lambda x: x.loc[x.groupby('qseqid')['length_minus_gapopen'].idxmax()])
                .assign(identical_seqs=lambda x: x.groupby('qseq')['qseqid'].transform('nunique'))
                .drop_duplicates(subset='qseq', keep='first')
                .sort_values(by='qseq') 
                .reset_index(drop=True)
            )
            blast_output_basename, extension = os.path.splitext(blast_output)
            write_fasta(blast_result, blast_filename=blast_output_basename + "_check.fasta")
            analysis_completed(finished_accession_numbers, analysed_accession_numbers, reuse_sequences)
            
    else:
        print(f"\033[31mNo blasting for {specie.replace('_', ' ')}: {file_to_blast} not found. \033[0m")

    if nothing_new_to_blast:
        with open(log_file, "a") as file:
            file.write(f'However, these species had no BLAST results from their sequences:\n')
            file.write('\n'.join(nothing_new_to_blast) + '\n')
        with open(f"BLAST_results/{run_name}_species_without_BLAST_results.txt", "w") as file:
            file.write('\n'.join(nothing_new_to_blast) + '\n')

    if reuse_sequences:
        os.remove(file_to_blast)
        os.remove(finished_accession_numbers)

def process_fasta_file(file_path, species_data, accessions_to_keep=None, log_file=None):
    """Parses, filters, and merges FASTA files into species_data dictionary."""
    for record in SeqIO.parse(file_path, "fasta"):
        # Extract accession number
        if '|' in record.id:
            accession = record.id.split('|')[1]
        elif '_' in record.id:
            accession = record.id.split('_')[1]
        else:
            accession = record.id

        # Filter unapproved/invalid sequences
        if accessions_to_keep is not None and accession not in accessions_to_keep:
            if log_file:
                with open(log_file, "a") as f: f.write(f"{accession}\n")
            continue
        
        sequence = str(record.seq)
        if not all(base in "ACGT" for base in sequence.upper()):
            if log_file:
                with open(log_file, "a") as f: f.write(f"{accession}\n")
            continue

        # Extract counter and species
        clean_header = record.description.rstrip('*')
        if '|' in clean_header:
            parts = clean_header.split('|')
            counter = int(parts[-1]) if parts[-1].isdigit() else 1
            species = parts[2].split(';')[-1]
            header_base = '|'.join(parts[:-1])
        else:
            parts = clean_header.split('_')
            counter = int(parts[-1]) if parts[-1].isdigit() else 1
            species = f"{parts[-3]}_{parts[-2]}"
            header_base = '_'.join(parts[:-1])

        # Merge into dict
        if species not in species_data:
            species_data[species] = {'sequences': {}}
        
        if sequence in species_data[species]['sequences']:
            species_data[species]['sequences'][sequence]['counter'] += counter
        else:
            species_data[species]['sequences'][sequence] = {
                'counter': counter, 'header_without_counter': header_base
            }

def filter_sequences(input_file, output_file, number_ambigious_nucleotides):
    """Filters sequences by ambiguity threshold and length bounds."""
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            seq = str(record.seq)
            seq_len = len(seq)

            # Check length constraints first
            if not (min_length <= seq_len <= max_length):
                continue

            if all(base in "ACGT" for base in seq):
                out_handle.write(f">{record.id}\n{seq}\n")
            else:
                if sum(1 for b in seq if b not in "ACGT") <= number_ambigious_nucleotides:
                    out_handle.write(f">{record.id}\n{seq}\n")

def add_symbol_to_new_sequence(new_potential_sequences):
    """Adds '*' to headers of new sequences for tree visualization."""
    marked = []
    with open(new_potential_sequences, 'r') as f:
        for line in f:
            if line.startswith('>'):
                marked.append(line.strip() + "*\n")
            else:
                marked.append(line)
    return marked

def align_batches(input_dir, output_dir, thread_number=7):
    """Aligns batches of FASTA files using MAFFT."""
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith(".fasta"):
            in_path = os.path.join(input_dir, file)
            out_path = os.path.join(output_dir, f"aligned_{file}")
            mafft_cline = MafftCommandline(input=in_path, thread=thread_number, auto=True)
            stdout, _ = mafft_cline()
            with open(out_path, "w") as f: f.write(stdout)

def merge_alignments(input_dir, output_file):
    """Merges aligned FASTA batches."""
    with open(output_file, "w") as outfile:
        for file in sorted(os.listdir(input_dir)):
            if file.startswith("aligned_") and file.endswith(".fasta"):
                with open(os.path.join(input_dir, file)) as infile:
                    outfile.write(infile.read())

def assess_alignment(file_path):
    """Prints basic alignment stats."""
    align = AlignIO.read(file_path, "fasta")
    print(f"Sequences: {len(align)} | Length: {align.get_alignment_length()}")
    gaps = sum(1 for i in range(align.get_alignment_length()) if '-' in align[:, i])
    print(f"Total gaps: {gaps}")

def filter_unique_fasta(file_path):
    """Removes non-unique sequences based on headers."""
    records = list(SeqIO.parse(file_path, "fasta"))
    counts = Counter([r.id for r in records])
    unique = [r for r in records if counts[r.id] == 1]
    SeqIO.write(unique, file_path, "fasta")

def write_output(output_path, species_data):
    """Writes the merged dictionary back to a FASTA file."""
    with open(output_path, 'w') as f:
        for species in sorted(species_data.keys()):
            for seq, info in species_data[species]['sequences'].items():
                f.write(f">{info['header_without_counter']}|{info['counter']}\n{seq}\n")

def read_sequences_to_remove(file_path):
    """Lists accessions to exclude based on curation file."""
    to_remove = set()
    if file_path and os.path.exists(file_path):
        with open(file_path, "r") as f:
            for line in f:
                if line.strip() and not line.startswith("Write down"):
                    to_remove.add(line.strip())
    return to_remove

def parse_fasta(file_path):
    """Counts total sequences per species."""
    counters = defaultdict(int)
    for record in SeqIO.parse(file_path, "fasta"):
        parts = record.description.split('|')
        counters[parts[2].split(';')[-1]] += int(parts[-1])
    return counters

def plot_histogram_sequence_per_species(species_counters, run_name):
    """Plots sequence count distribution."""
    vals = sorted(species_counters.values(), reverse=True)
    mean_v, med_v = np.mean(vals), np.median(vals)
    plt.figure(figsize=(10, 6))
    plt.bar(range(1, len(vals) + 1), vals)
    plt.axvline(x=len(vals)/2, color='r', linestyle='--', label=f'Median: {med_v:.1f}')
    plt.title('Sequences per species')
    plt.savefig(f'{run_name}_histogram_sequences_per_species.png')
    plt.close()

def plot_histogram_sequence_lenghts(updated_database, run_name):
    """Plots sequence length distribution."""
    lengths = [len(r.seq) for r in SeqIO.parse(updated_database, "fasta")]
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=50)
    plt.title('Sequence length distribution')
    plt.savefig(f'{run_name}_histogram_sequence_lengths.png')
    plt.close()

def find_duplicates(reference_database):
    """Identifies duplicate sequences in a FASTA file and groups their headers."""
    sequence_dict = defaultdict(list)
    
    with open(reference_database, "r") as file:
        current_header = None
        current_sequence = []
        
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    seq_str = "".join(current_sequence)
                    sequence_dict[seq_str].append(current_header)
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
                
        # Catch the final sequence block
        if current_header and current_sequence:
            seq_str = "".join(current_sequence)
            sequence_dict[seq_str].append(current_header)
            
    return {seq: headers for seq, headers in sequence_dict.items() if len(headers) > 1}

def get_positions(sequence_input):
    """Get a list of positions used for checking mismatches."""
    return list(range(len(sequence_input)))

# Degenerate bases mapping
degenerate_base_matches = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
    'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}, 'I': {'A', 'C', 'G', 'T'}
}

def positions_to_check(mode, positions, aligned_primer, aligned_sequence, record_result):
    """Shows if there is a mismatch at a particular position."""
    number_of_mismatches = 0
    mismatch_details = []

    mismatches_mode = f"{mode}_mismatches"
    count_mode = f"{mode}_number_of_mismatches"
    
    for pos in positions:
        sequence_char = aligned_sequence[pos]
        primer_char = aligned_primer[pos]

        if sequence_char != primer_char:
            # Check regular or degenerate character matches
            if primer_char in degenerate_base_matches and sequence_char in degenerate_base_matches[primer_char]:
                continue  # Valid degenerate match
            else:
                number_of_mismatches += 1
                bad_position = len(positions) - pos 
                mismatch_details.append((bad_position, primer_char, sequence_char))

    record_result[count_mode] = number_of_mismatches
    record_result[mismatches_mode] = mismatch_details if mismatch_details else "NA"

def primer_alignment(mode, primer_sequence, sequence_to_align, record_result):
    """Align the primer with the sequence and check for mismatches."""
    aligner = PairwiseAligner()
    aligner.match = 1
    
    # Check if the primer contains degenerate bases
    standard_bases = set("ACGT")
    is_degenerate = any(base not in standard_bases for base in primer_sequence.upper())
    
    # Apply mismatch penalty dynamically
    aligner.mismatch = 0 if is_degenerate else -1  
    
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    aligner.end_gap_score = 0
    aligner.mode = 'global'
    
    alignments = aligner.align(primer_sequence, sequence_to_align)
    primer_length = len(primer_sequence)
    max_score = aligner.match * primer_length
    score_cutoff = max_score * 0.5
    length_cutoff = primer_length * 1.3

    sequence_key = f"{mode}_sequence"
    
    best_alignment = None 
    best_score = -float('inf')

    for alignment in alignments:
        if alignment.score > score_cutoff and alignment.length <= length_cutoff:
            if alignment.score > best_score:
                best_score = alignment.score
                best_alignment = alignment

    if best_alignment:
        aligned_primer_full = best_alignment[0]
        aligned_sequence_full = best_alignment[1]
        
        # Strip leading/trailing gaps from the primer alignment
        start = 0
        end = len(aligned_primer_full)
        while start < end and aligned_primer_full[start] == '-':
            start += 1
        while end > 0 and aligned_primer_full[end - 1] == '-':
            end -= 1
            
        aligned_primer = aligned_primer_full[start:end]
        aligned_sequence = aligned_sequence_full[start:end]

        record_result[sequence_key] = aligned_sequence
        positions = get_positions(aligned_primer)
        positions_to_check(mode, positions, aligned_primer, aligned_sequence, record_result)
    else:
        record_result[sequence_key] = "NA"
        record_result[f"{mode}_number_of_mismatches"] = "NA"
        record_result[f"{mode}_mismatches"] = "NA"

def get_family_name(leaf, pattern=r'\.\d+_|[\d.]+'):
    """
    Extracts the family name from a leaf node's name in the phylogenetic tree.
    
    Assumes a standardized FASTA header format where the lineage string 
    is parsed and the family identifier sits at index 5.
    """
    name_parts = re.split(pattern, leaf.name)
    parts = [part for part in name_parts if part]
    
    # parts[1] contains the lineage string separated by underscores
    return parts[1].split('_')[5]


def is_monophyletic(tree, species_leaves):
    """
    Determines if a given set of leaves forms a monophyletic group (clade).
    
    A group is monophyletic if their lowest common ancestor only has 
    these specific leaves as descendants.
    """
    # Special case: A species with only one sequence is inherently monophyletic
    if len(species_leaves) == 1 and species_leaves[0].is_leaf():
        return True

    # Find the lowest common ancestor of the specified species leaves
    common_ancestor = tree.get_common_ancestor(species_leaves)
    
    # Extract the names of all leaves descending from this ancestor
    ancestor_leaves = set(common_ancestor.get_leaf_names())
    target_leaves = set(leaf.name for leaf in species_leaves)
    
    # If the ancestor's descendants exactly match our target list, it is monophyletic
    return ancestor_leaves == target_leaves

def check_mismatches(mismatches, position_checks, mismatch_sets):
    """Checks for mismatches between the input primer and sequences."""
    for item in mismatches:
        if isinstance(item, tuple):
            position, primer_char, sequence_char = item
            if position in position_checks:
                if primer_char in degenerate_base_matches and sequence_char in degenerate_base_matches[primer_char]:
                    continue
                elif (primer_char, sequence_char) in mismatch_sets:
                    return True
    return False

def has_gap(mismatches):
    """Checks if there are any gaps in either the primer or in the sequence."""
    for item in mismatches:
        if isinstance(item, tuple):
            _, primer_char, sequence_char = item
            if primer_char == '-' or sequence_char == '-':
                return True
    return False

def nucleotide_proportions_diagram(df, mode):
    """Generates a stacked bar diagram showing nucleotide proportions."""
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    df['Most_Common'] = df.idxmax(axis=1)
    plot_data = df.drop(columns='Most_Common')

    fig, ax = plt.subplots(figsize=(12, 8))
    plot_data.plot(kind='bar', stacked=True, ax=ax, colormap='coolwarm')

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Most_Common'], rotation=0, ha='center')

    plt.xlabel("Position 5' to 3'")
    plt.ylabel('Proportion')
    plt.title('Nucleotide proportions at each position')
    plt.legend(title='Nucleotide/Gaps', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout() 

    file_name = "Forward_primer_nucleotide_proportions.png" if mode == 'forward_sequence' else "Reverse_primer_nucleotide_proportions.png"
    plt.savefig(file_name)
    plt.close(fig) # Free up memory

def proportions_table(df, mode, padded_primer):
    """Generates a CSV table of nucleotide frequencies."""
    file_name = "Forward_primer_nucleotide_frequencies.csv" if mode == "forward_sequence" else "Reverse_primer_nucleotide_frequencies.csv"
    
    with open(file_name, 'w') as f:
        # Write custom headers safely without file-seeking hacks
        headers = ['Consensus_sequence_5_to_3'] + list(df.columns)
        f.write(';'.join(str(h) for h in headers) + '\n')
        
        # Write padded primer row
        f.write(';'.join(padded_primer) + '\n')
        
        # Append DataFrame (skipping its default header)
        df.to_csv(f, sep=";", header=False)

def multiple_sequence_alignment(df, mode): 
    """Generates sequence alignments using the linsi algorithm with MAFFT."""
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as temp_input:
        for index, row in df.iterrows():
            if row[mode] != "NA":
                temp_input.write(f">{row['species_name']}_{row['accession_number']}\n{row[mode]}\n")
        temp_input_path = temp_input.name
    
    file_name = "Forward_sequences_aligned.fasta" if mode == "forward_sequence" else "Reverse_sequences_aligned.fasta"
    
    # Run MAFFT alignment
    mafft_cline = MafftCommandline(input=temp_input_path, auto=True)
    stdout, stderr = mafft_cline()
    
    with open(file_name, "w") as handle:
        handle.write(stdout)
        
    os.remove(temp_input_path)

def analyze_primer_frequencies(df, mode, primer): 
    """Calculates nucleotide frequencies and generates reports. (Formerly omega_function)"""
    multiple_sequence_alignment(df, mode) 
    
    # Filter out NA sequences
    sequences = [seq for seq in df[mode] if seq != "NA"]
    if not sequences:
        return
        
    max_len = max(len(seq) for seq in sequences) 
    
    # Add gaps at the end of sequences that are shorter than max_len
    padded_sequences = [seq.ljust(max_len, '-') for seq in sequences] 
    adjusted_primer = primer.ljust(max_len, '-')
    
    padded_primer = ["Primer_used_5_to_3"] + list(adjusted_primer)
    
    frequencies = [] 
    for i in range(max_len):
        column = [seq[i] for seq in padded_sequences]
        frequencies.append(Counter(column))

    proportions = []
    most_common_seq = []
    
    for freq in frequencies:
        total = sum(freq.values())
        if total > 0:
            proportions.append({char: count / total for char, count in freq.items()})
        else:
            proportions.append({})
        most_common_seq.append(freq.most_common(1)[0][0] if freq else '-')

    proportions_df = pd.DataFrame(proportions).fillna(0)

    nucleotide_order = ['A', 'C', 'G', 'T', '-']
    for col in nucleotide_order:
        if col not in proportions_df.columns:
            proportions_df[col] = 0
            
    proportions_df = proportions_df[nucleotide_order]
    
    df_to_transpose = pd.DataFrame(frequencies).fillna(0).astype(int)
    for col in nucleotide_order:
        if col not in df_to_transpose.columns:
            df_to_transpose[col] = 0
    df_to_transpose = df_to_transpose[nucleotide_order]

    df_transposed = df_to_transpose.T
    new_headers = [freq.most_common(1)[0][0] if freq else '-' for freq in frequencies]
    df_transposed.columns = new_headers

    nucleotide_proportions_diagram(proportions_df, mode)
    proportions_table(df_transposed, mode, padded_primer)


# =============================================================================
# PIPELINE SUBCOMMANDS
# =============================================================================

def run_template(args):
    """Executes the Reference Template logic."""
    log_file, run_name = get_log_file_name()
    program_start = time.time()
    
    log_start(log_file, "Template Creation")
    
    append_and_print_message(log_file, f"\nRunning template creation...\nCommand: {get_command_string()}\n")
    
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
    thread_number = number_threads(args.threads)
    command_string = get_command_string()
    output_prefix = "temp_split_file"
    forward_primer = args.forward.lower()
    forward_primer = forward_primer.replace('i', 'n')
    reverse_primer = reverse_complement(args.reverse).lower()
    reverse_primer = reverse_primer.replace('i', 'n')

    if not args.Complete:
        if not all ([args.input_file, args.forward, args.reverse, args.email, args.api_key]):
            print("An input file with a list of species along with the following are required:\nForward primer (-f, --forward)\nReverse primer (-r, --reverse)\nEmail address (-e, --email)\nAPI key (-a, --api_key)\n\nThe flag to finish the database (-C, --Complete) requires no other inputs.\n")
            sys.exit(1)
        else:
            input_file = os.path.abspath(args.input_file)
            save_config(args.email, args.api_key, args.forward, args.reverse)
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
                            f'Search term: [Organism] AND ("{length_threshold}"[SLEN] : "{max_length}"[SLEN]) AND biomol_genomic[PROP] NOT "unverified" {args.query}\n')
                
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


            finally:
                
                os.chdir("..")

                print("\nA draft of the reference template database has been created.")
                print(f"Make sure to review {working_directory}{to_be_curated} before finalizing it to a reference template database with -C (--Complete).")
                print("reference_template.fasta has been created and is ready for use with Echopipe_database_creation.py.\n")

                next_step_cmd = f"python echopipe.py template -C {args.input_file}"

                log_end(log_file, program_start, f"\033[32m{next_step_cmd}\033[0m")

    else:

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

        if args.Complete:
            next_step_cmd = f"python echopipe.py create {input_file} reference_template_database.fasta"

            log_end(log_file, program_start, f"\033[32m{next_step_cmd}\033[0m")
        
        else:
            next_step_cmd = f"python echopipe.py template -C {input_file}"
            
            # Combine the tutorial text into the log_end message
            msg = (
                "First remove sequences arising from other gene regions; see tutorial: "
                "https://github.com/EivindStensrud/EchoPipe/tree/main\n\n"
                "Next code line:\n"
                f"\033[32m{next_step_cmd}\033[0m"
            )
            log_end(log_file, program_start, msg)

def run_create(args):
    """Executes the Database Creation logic."""
    log_file, run_name = get_log_file_name()
    program_timer = time.time()
    
    log_start(log_file, "Database Creation")
    
    append_and_print_message(log_file, f"\nRunning database creation...\nCommand: {get_command_string()}\n")

    config = load_config()
    if args.email == "email@email.com": args.email = config.get("email", args.email)
    if args.api_key == "api_key": args.api_key = config.get("api_key", args.api_key)

    if not args.repeat:
        if args.email == "email@email.com" or args.api_key == "api_key" or not all([args.input_file, args.input_database]):
            print(f"\nThe command used to run the script was: python {get_command_string()}\n")
            print("ERROR: An input file with a list of species, input database, email address and API key are required.\n(Email address and API key are not needed when re-analyzing sequences with -R, --repeat).")
            sys.exit(1)

    input_species = args.input_file
    input_fasta = args.input_database
    Entrez.email = args.email
    Entrez.api_key = args.api_key
    sort_by_sequence_length = args.sort
    max_search = args.maxcount
    max_length_input = args.maxlength
    min_length = args.ampliconsize
    max_length_search = f' AND ("{str(min_length)}"[SLEN] : "{str(max_length_input)}"[SLEN])'
    mitochondria_on = ' AND mitochondrion[filter]' if args.mitochondria else ""
    ribosomal_on = ' AND 12S' if args.ribosomal else ""
    custom_query = f' {args.query}' if args.query else ""
    sequence_batch_size = args.batch_size
    use_old_taxid = args.taxid
    evalue = f"5e-{args.evalue}"
    reuse_sequences = args.repeat

    logs_dir = "Log_files/"
    thread_number = number_threads(args.threads)

    counter_file = os.path.join(logs_dir, "run_counters.txt")
    DATE = datetime.today().strftime('%Y-%m-%d')
    i = read_counter(counter_file, DATE)
    update_counter(counter_file, DATE, i + 1 )
    
    log_file, run_name = get_log_file_name(logs_dir)
    print(f"Log file for this run: {os.path.relpath(log_file)}")

    database_directory = "Template_databases/" + os.path.splitext(input_fasta)[0] 
    search_settings = '[Organism] AND biomol_genomic[PROP]' + max_length_search + mitochondria_on + ribosomal_on + custom_query 
    
    base_dir = os.path.abspath(os.getcwd())
    taxid_collection = os.path.join(base_dir, "taxid_collection.txt") 
    logs_dir_abs = os.path.join(base_dir, logs_dir)

    taxid_log_file = os.path.join(logs_dir_abs, f"{run_name}_taxid_collection.txt")  
    error_collection = os.path.join(logs_dir_abs, f"{run_name}_species_not_found.txt") 
    duplicate_collection = os.path.join(logs_dir_abs, f"{run_name}_duplicate_species_entries.txt") 
    species_found = os.path.join(logs_dir_abs, f"{run_name}_species_found.txt") 

    taxonomic_ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus'] 
    create_directory(database_directory) 
    timestamp = datetime.today().strftime("%H:%M:%S") 

    with open(log_file, "a") as file:
        file.write(
           f"\nThe command used to run the script was: python {get_command_string()}\n\n"
            f"\nReference database was created with the following settings:\n"
            f"This log file contains information from the run {run_name} at {timestamp}.\n\n"
            "------------------------------------------------------------------------------\n\n"
            "The script was run with the following settings:\n"
            f"Input file: {args.input_file}\n"
            f"Input reference database: {args.input_database}\n"
            f"Email address: {args.email}\n"
            f"API key: {args.api_key}\n"
            f"Date: {DATE}\n"
        )

        if not reuse_sequences:
            file.write(
               f"Sorted by largest sequence lengths: {args.sort}\n"
                f"Maximum accession number per species: {args.maxcount}\n"
                f"Minimum amplicon size: {args.ampliconsize}\n"
                f"Search term: {search_settings}\n\n"
                "--------------------------------------------------------------------------\n\n"
            )
        else:
            file.write(
                f"Existing sequence files are BLASTed against the input database.\n\n"
                "--------------------------------------------------------------------------\n\n"
            )

    database_name = os.path.splitext(input_fasta)[0] 
    database_path = os.path.join(database_directory, database_name) 

    makeblastdb_cline = NcbimakeblastdbCommandline( 
        dbtype="nucl", 
        input_file=input_fasta, 
        out=database_path, 
        title=input_fasta 
    )
    
    stdout, stderr = makeblastdb_cline() 
    if stderr:
        append_and_print_message(log_file, f"Error: {stderr}\n")
    else:
        append_and_print_message(log_file,
            f"Template database creation completed successfully in directory: {database_directory}\n\n"
            "--------------------------------------------------------------------------\n\n")

    consecutive_fail_counter = 0 
    max_retries = 5

    if not use_old_taxid: 
        type_of_download = "Species information download"
        append_and_print_message(log_file, f"{type_of_download} is starting...\n")

        with open(input_species, "r") as file:
            input_species_list = file.readlines()

        input_species_list = sorted([line.strip().split(";")[0] for line in input_species_list if line.strip()])
        input_species_list = sorted(list(set(input_species_list)))
        print(f"Input_species_list: {input_species_list}")

        species_not_found = []
        duplicate_species = []
        filtered_species_collection = []
        failed_information_downloads = []

        loop_timer = time.time()

        with ThreadPoolExecutor(max_workers=thread_number) as executor:
            futures = [
                executor.submit(process_species_taxonomy, species, filtered_species_collection, taxonomic_ranks, max_retries) 
                for species in input_species_list
            ]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing species taxonomy"):
                result = future.result()
                if 'species_not_found' in result:
                    species_not_found.append(result['species_not_found'] + "\n")
                elif 'duplicate_species' in result:
                    duplicate_species.append(result['duplicate_species'] + "\n")
                elif 'filtered_species_collection' in result:
                    filtered_species_collection.append(result['filtered_species_collection'] + "\n")
                elif 'failed_information_downloads' in result:
                    failed_information_downloads.append(result['failed_information_downloads'])

        loop_duration = round(time.time() - loop_timer, 2)
        loop_finished(failed_information_downloads, log_file, type_of_download, loop_duration)

        if species_not_found:
            species_not_found.insert(0, f"The following species had no taxonomic search results with the input species list used on the run {run_name}.\n\n")
            write_lines_from_list(species_not_found, error_collection, log_file, f"Non-valid species name found. See {os.path.relpath(error_collection)}\n")
        if duplicate_species:
            duplicate_species.insert(0, f"The following species were duplicates from the input species list from the run {run_name}.\n\n")
            write_lines_from_list(duplicate_species, duplicate_collection, log_file, f"Duplicate species name found. See {os.path.relpath(duplicate_collection)}\n")
        write_lines_from_list(filtered_species_collection, taxid_collection, log_file, f"The following species had sequences on NCBI. See {os.path.relpath(species_found)}\n")
        print(f"{type_of_download} took {loop_duration} seconds to finish.")

    database_curation = f"Database_curation/{run_name}/"
    create_directory(database_curation) 

    with open(taxid_collection, "r") as file:
        species_list = file.readlines() 
        
    write_lines_from_list(species_list, taxid_log_file) 
    species_list = sorted([line.strip().split("; ")[2].replace(' ', '_') for line in species_list]) 

    if not reuse_sequences:
        type_of_download = "Accession number download"
        append_and_print_message(log_file, f"{type_of_download} is starting...\n")
                                 
        over_cap_csv_list = [] 
        failed_accession_downloads = [] 
        create_directory("Species") 
        consecutive_fail_counter = 0
        failed_information_downloads = [] 

        batch_size = 10000 
        if max_search <= 10000: 
            batch_size = max_search

        loop_timer = time.time() 
        consecutive_fail_counter = 0

        with ThreadPoolExecutor(max_workers=thread_number) as executor:
            futures = [executor.submit(process_species_accession_number, specie, search_settings, max_search,
                               batch_size, max_retries, sort_by_sequence_length, consecutive_fail_counter, log_file,
                               failed_accession_downloads, type_of_download, run_name, create_directory,
                               repeated_failures, error_occurrence, over_cap_csv_list)
                for specie in species_list]

            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing species accession number"):
                result = future.result() 
        
        loop_duration = round(time.time() - loop_timer, 2)
        loop_finished(failed_information_downloads, log_file, type_of_download, loop_duration)

        write_lines_from_list(over_cap_csv_list, f"{run_name}_over_max.csv", log_file,
            f"Note: Not all sequences were downloaded for one or more species since --maxcount ({max_search}) < available accession numbers. See {run_name}_over_max.csv\n"
            f"Suggestion: Run program again with {run_name}_over_max.csv as input after database curation.\n"
            f"Make sure to adjust parameters such as --maxcount and/or --sort by sequence length.\n\n")

        accession_collection = [] 
    else: 
        append_and_print_message(log_file,
            "Existing sequences are to be re-analysed.\n"
            "No accession numbers are downloaded.\n\n"
            "--------------------------------------------------------------------------\n\n\n")

    print(f"{type_of_download} took {loop_duration} seconds to finish.")

    analysed = "analysed.accession_numbers" 
    new_accession_species = [] 

    if not reuse_sequences:
        for species in species_list: 
            accession_file_1 = f"Species/{species}/{species}_{run_name}_accession_numbers.txt" 
            accession_file_2 = f"Species/{species}/{species}_{analysed}.txt" 
            does_file_exist(accession_file_2) 
            
            with open(accession_file_1, "r") as file: 
                accession_file_content_1 = set(file.readlines()) 
            with open(accession_file_2, "r") as file: 
                accession_file_content_2 = set(file.readlines())

            differences = sorted(accession_file_content_1.difference(accession_file_content_2)) 

            if not differences: 
                print(f"No new accession numbers detected for \033[31m{species.replace('_', ' ')}\033[0m")
            else:
                print(f"New accession numbers detected for \033[32m{species.replace('_', ' ')}\033[0m") 
                new_accession_species.append(species) 
                new_accession = f"Species/{species}/{species}_new_accession_numbers_{run_name}.txt" 
                write_lines_from_list(differences, new_accession)

        append_and_print_message(log_file,
            f"{len(new_accession_species)} out of {len(species_list)} species have new accession numbers.\n\n"
            "--------------------------------------------------------------------------\n\n\n")
    else: 
        for species in species_list:
            accession_numbers_file = f"Species/{species}/{species}_{analysed}.txt"
            previous_accession_numbers_path = f"Species/{species}"
            
            check_if_analysed_exist = os.path.isfile(accession_numbers_file)
            
            if os.path.exists(previous_accession_numbers_path):
                previous_accession_numbers_exist = any("new_accession_numbers" in file for file in os.listdir(previous_accession_numbers_path))
            else:
                previous_accession_numbers_exist = False
                
            if not check_if_analysed_exist and not previous_accession_numbers_exist:
                continue 
            else:
                new_accession_species.append(species)

    if not reuse_sequences:
        type_of_download = "Sequence download"
        append_and_print_message(log_file, f"{type_of_download} is starting...\n")

        failed_sequence_download = []
        consecutive_fail_counter = 0
        failed_information_downloads = [] 

        loop_timer = time.time()
        
        with ThreadPoolExecutor(max_workers=thread_number) as executor:
            futures = [executor.submit(download_fasta_sequences, specie, run_name, taxid_collection, sequence_batch_size, max_retries, log_file, consecutive_fail_counter, failed_information_downloads, type_of_download)
                for specie in new_accession_species]

            for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading fasta sequences"):
                result = future.result() 
        
        loop_duration = round(time.time() - loop_timer, 2)
        loop_finished(failed_information_downloads, log_file, type_of_download, loop_duration)

    else: 
        append_and_print_message(log_file,
            "Existing sequences are to be re-analysed.\n"
            "No sequences are downloaded.\n\n"
            "--------------------------------------------------------------------------\n\n\n")

    print(f"{type_of_download} took {loop_duration} seconds to finish.\n"
        "--------------------------------------------------------------------------\n\n")

    loop_timer = time.time()
    type_of_download = "BLAST extraction"
    outfmt_string = "6 qseqid qseq sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    append_and_print_message(log_file, f"{type_of_download} is starting...\n")
    
    create_directory("BLAST_results")

    with ThreadPoolExecutor(max_workers=thread_number) as executor:
        futures = [executor.submit(perform_blast_operations, new_accession_species, run_name, database_path, analysed, min_length,
        outfmt_string, evalue, reuse_sequences, log_file, specie)
            for specie in new_accession_species]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Conducting BLAST"):
            result = future.result() 
        
        loop_duration = round(time.time() - loop_timer, 2)

    loop_finished(failed_information_downloads, log_file, type_of_download, loop_duration)

    print(f"{type_of_download} took {loop_duration} seconds to finish.\n"
        "--------------------------------------------------------------------------\n\n")

    all_species_new_results = [] 

    for species in new_accession_species: 
        new_results_check = f"Species/{species}/{species}_{run_name}_fasta_blast_results_check.fasta" 
        
        if os.path.exists(new_results_check): 
            with open(new_results_check, "r") as file: 
                new_blast_results = file.readlines() 
                all_species_new_results.extend(new_blast_results) 
        else:
            print(f"Nothing new from {run_name} for {species}.") 

    if all_species_new_results: 
        result_file = f"BLAST_results/{run_name}_to_curate.fasta"
        with open(result_file, "w") as file: 
            for line in all_species_new_results:
                file.write(line)
        append_and_print_message(log_file,
            f"\n\nNew sequence(s) found! {result_file} is ready to be curated.\n")

        species_in_final_file = set()  

        with open(result_file, "r") as file:
            final_file_content = file.readlines()

            for line in final_file_content:
                if line.startswith('>'):  
                    header = line.strip()
                    species_name = header.split('|')[2].split(';')[-1].replace('_', ' ')
                    species_in_final_file.add(species_name)

        sorted_final_file = sorted(species_in_final_file)

        with open(f"{database_curation}{run_name}_species_found.txt", "w") as output_file:
            for species_name in sorted_final_file:
                output_file.write(species_name + "\n")

    error_collection_date = f"{database_curation}{run_name}_species_not_found.txt"
    duplicate_collection_date = f"{database_curation}{run_name}_duplicate_species_entries.txt"

    if species_not_found:
        with open(error_collection, "r") as outfile:
            with open(error_collection_date, "w") as infile:
                for line in outfile:
                    infile.writelines(line)
    if duplicate_species:
        with open(duplicate_collection, "r") as outfile:
            with open(duplicate_collection_date, "w") as infile:
                for line in outfile:
                    infile.writelines(line)
    
    next_cmd = f"python echopipe.py curate BLAST_results/{run_name}_to_curate.fasta"
    
    # print recommended next code line
    recommendation_msg = (
        "To help curate the database (you can optionally filter by sequence length using --min_length and --max_length), run the following code:\n"
        f"\033[32m{next_cmd}\033[0m"
    )
    
    log_end(log_file, program_timer, recommendation_msg)


def run_curate(args):
    """Executes the Database Curation logic (Alignments & Trees)."""
    log_file, run_name = get_log_file_name()
    program_timer = time.time()
    
    log_start(log_file, "Database Curation")

    # 2. Paths & Vars
    new_potential_sequences = os.path.abspath(args.input_file)
    most_recent_database = os.path.abspath(args.old_database) if args.old_database else ""
    mafft_online_alignment = os.path.abspath(args.mafft_online) if args.mafft_online else None
    num_ns = args.number_ns
    thread_number = number_threads()
    
    # 3. Directories
    database_curation_dir = os.path.abspath("Database_curation/")
    working_directory = os.path.join(database_curation_dir, run_name)
    os.makedirs(working_directory, exist_ok=True)
    
    # Files
    concatenated_file = os.path.join(working_directory, f"{run_name}_concatenated_file.fasta")
    maffted_file = os.path.join(working_directory, f"{run_name}_aligned.fasta")
    tree_string = os.path.join(working_directory, f"{run_name}_tree_string.newick")
    temp_filtered_file = os.path.join(working_directory, "temp_sequence_file.fasta")

    # 4. Preparation
    filter_unique_fasta(new_potential_sequences)
    base_dir = os.getcwd()
    os.chdir(working_directory)

    # 5. Alignment Logic
    min_len = args.min_length
    max_len = args.max_length
    if mafft_online_alignment and os.path.exists(mafft_online_alignment):
        append_and_print_message(log_file, "Using provided MAFFT online alignment.\n")
        concatenated_file = mafft_online_alignment
        filter_sequences(new_potential_sequences, temp_filtered_file, num_ns, min_len, max_len)
    else:
        filter_sequences(new_potential_sequences, temp_filtered_file, num_ns, min_len, max_len)        
        records_to_write = []
        
        # 1. Parse new potential sequences, add '*', and clean forbidden characters from IDs
        for record in SeqIO.parse(temp_filtered_file, "fasta"):
            # Append '*' to mark new sequences for tree visualization
            record.id = f"{record.id}*"
            record.description = ""  # Clear description to avoid redundant text
            record.id = record.id.replace('|', '_').replace(';', '_').replace(' ', '')
            records_to_write.append(record)
            
        # 2. Parse the old/existing database records if available
        if os.path.exists(most_recent_database):
            for record in SeqIO.parse(most_recent_database, "fasta"):
                record.id = record.id.replace('|', '_').replace(';', '_').replace(' ', '')
                records_to_write.append(record)
                
        # 3. Write out using SeqIO to guarantee perfect newline separation every time
        SeqIO.write(records_to_write, concatenated_file, "fasta")

        n = sum(1 for line in open(concatenated_file) if line.startswith(">"))
        
        # Determine MAFFT algorithm
        mafft_args = {"input": concatenated_file, "thread": thread_number, "reorder": True}
        if n <= 1000: mafft_args.update({"localpair": True, "maxiterate": 100})
        else: 
            mafft_args.update({"auto": True})
        
        mafft_cline = MafftCommandline(**mafft_args)
        stdout, _ = mafft_cline()
        with open(maffted_file, "w") as f:
            f.write(stdout)

    # 6. Tree & Analysis
    dup_seqs = find_duplicates(concatenated_file)
    if dup_seqs:
        with open(f"{run_name}_duplicate_sequences.txt", "w") as f:
            for seq, headers in dup_seqs.items():
                f.write("\n".join(headers) + f"\n{seq}\n\n")

    # FastTree
    is_large = n >= 10000
    fasttree_kwargs = {
        "nt": True,
        "fastest": is_large,
        "input": maffted_file,
        "out": tree_string
    }
    if is_large:
        fasttree_kwargs["cat"] = 20

    fasttree_cline = FastTreeCommandline(**fasttree_kwargs)
    stdout, stderr = fasttree_cline()
    
    if stderr:
        append_and_print_message(log_file, f"FastTree messages: {stderr.strip()}\n")
    # Monophyly Check
    tree = PhyloTree(tree_string)
    pattern = r'\.\d+_|[\d.]+'
    species_leaves = defaultdict(list)
    for leaf in tree: species_leaves[get_species_name(leaf, pattern)].append(leaf)

    mono, para = [], []
    table_data = {'Species': [], 'Family': [], 'Is Monophyletic': []}
    
    for sp, leaves in species_leaves.items():
        is_mono = is_monophyletic(tree, leaves)
        if is_mono: mono.append(sp)
        else: 
            para.append(sp)
        for fam in set(get_family_name(l, pattern) for l in leaves):
            table_data.update({'Species': table_data['Species']+[sp], 'Family': table_data['Family']+[fam], 'Is Monophyletic': table_data['Is Monophyletic']+[is_mono]})

    pd.DataFrame(table_data).to_csv(f"{run_name}_dataframe.csv", index=False)
    
    # Save text outputs
    with open(f"{run_name}_monophyletic_group.txt", "w") as f:
        f.write("Species found as monophyletic:\n\n" + "\n".join(sorted(mono)) + "\n")
    with open(f"{run_name}_paraphyletic_group.txt", "w") as f:
        f.write("Species found as paraphyletic:\n\n" + "\n".join(sorted(para)) + "\n")

    # Clean up
    if os.path.exists(temp_filtered_file): os.remove(temp_filtered_file)
    os.chdir(base_dir)

    # 7. Final Recommendation
    next_step = f"\033[32mpython echopipe.py complete -b BLAST_results/{run_name}_to_curate.fasta -c Database_curation/{run_name}/{run_name}_aligned.fasta -u Database_name_{run_name}.fasta\033[0m"
    log_end(log_file, program_timer, next_step)

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


def run_complete(args):
    """Executes the Database Completion logic (Merging, Filtering, Trees & Stats)."""
    # 1. Setup paths and directories
    log_file, run_name = get_log_file_name()
    program_timer = time.time()
    
    log_start(log_file, "Database Completion")
    
    new_db_content = os.path.abspath(args.blast_file) if args.blast_file else ""
    aligned_curated_file = os.path.abspath(args.curated_file) if args.curated_file else ""
    updated_database = os.path.abspath(args.updated_database)
    old_database = os.path.abspath(args.old_database) if args.old_database else ""
    thread_number = number_threads()

    command_string = get_command_string()
    curated_directory = os.path.abspath(f"Database_curation/{run_name}/Curated_content")
    os.makedirs(curated_directory, exist_ok=True)
    
    removed_accession_log = os.path.abspath(f"Database_curation/{run_name}/{run_name}_removed_accessions.txt")

    append_and_print_message(log_file, f"\nRunning database completion...\nCommand: {command_string}\n")

    # 2. Process and merge sequences
    species_data = defaultdict(lambda: {'sequences': {}})
    accessions_to_delete = read_sequences_to_remove(aligned_curated_file)
    accessions_to_keep = read_sequences_to_keep(aligned_curated_file) if aligned_curated_file else None

    if old_database:
        print(f"Loading old database: {os.path.basename(old_database)}")
        process_fasta_file(old_database, species_data)

    if new_db_content:
        print(f"Loading and filtering new sequences: {os.path.basename(new_db_content)}")
        process_fasta_file(new_db_content, species_data, accessions_to_keep, removed_accession_log)
        
    if updated_database:
        print(f"Writing updated database to: {os.path.basename(updated_database)}")
        write_output(updated_database, species_data)

    # 3. Post-Processing & Statistics (Switching to curated directory)
    base_directory = os.getcwd()
    os.chdir(curated_directory)
    print(f"\nGenerating informative files in directory: {curated_directory}\n")

    duplicate_sequences = find_duplicates(updated_database)
    if duplicate_sequences:
        duplicate_sequences_file = f"{run_name}_duplicate_sequences.txt"
        with open(duplicate_sequences_file, "w") as file:
            for sequence, headers in duplicate_sequences.items():
                for header in headers:
                    file.write(header + "\n")
                file.write(sequence + "\n\n")
        print(f"Duplicate sequences that persisted are noted in {duplicate_sequences_file}.\n")

    # Plotting
    species_counters = parse_fasta(updated_database)
    plot_histogram_sequence_per_species(species_counters, run_name)
    plot_histogram_sequence_lenghts(updated_database, run_name)

    # 4. Alignment Preparation
    temp_file = os.path.abspath("temporary_file.fasta")
    maffted_database = os.path.abspath(f"aligned_{os.path.basename(updated_database)}")
    tree_string = os.path.abspath(f"{run_name}_curated_tree_string.newick")

    with open(updated_database, "r") as outfile:
        n = sum(1 for line in outfile if line.startswith(">"))
        print(f"Total number of sequences for alignment: {n}")
        
        outfile.seek(0)
        with open(temp_file, "w") as infile:
            for line in outfile:
                infile.write(line.replace("|", "_").replace(";", "_"))

    # 5. MAFFT Alignment
    sequences_per_file_linsi = 1000
    if n <= sequences_per_file_linsi:
        append_and_print_message(log_file, f"Aligning with MAFFT linsi (accurate) using {n} sequences...")
        mafft_cline = MafftCommandline(input=temp_file, localpair=True, maxiterate=100, thread=thread_number, reorder=True)
    elif n <= 10000:
        append_and_print_message(log_file, f"Aligning with MAFFT auto using {n} sequences with {thread_number} threads...")
        mafft_cline = MafftCommandline(input=temp_file, auto=True, thread=thread_number, reorder=True)
    else:
        append_and_print_message(log_file, f"Aligning with MAFFT auto using {n} sequences...")
        mafft_cline = MafftCommandline(input=temp_file, auto=True, reorder=True)

    try:
        stdout, stderr = mafft_cline()
        with open(maffted_database, "w") as out_handle:
            out_handle.write(stdout)
        print("Alignment complete.\n")
    except Exception as e:
        print(f"MAFFT error: {e}")
        print(f"Try with larger RAM, or upload {os.path.basename(temp_file)} to MAFFT Online: https://mafft.cbrc.jp/alignment/server/large.html")
        print("Recommended MAFFT online settings:\n- UPPERCASE/lowercase: Same as input\n- Strategy: PartTree\n")

    mafft_duplicates = find_duplicates(maffted_database)
    if mafft_duplicates:
        dup_mafft_file = f"{run_name}_aligned_duplicate_sequences.txt"
        with open(dup_mafft_file, "w") as file:
            for sequence, headers in mafft_duplicates.items():
                for header in headers:
                    file.write(header + "\n")
                file.write(sequence + "\n\n")
        print(f"Aligned duplicate sequences are noted in {dup_mafft_file}.")
    else:
        print("No duplicate sequences found in alignment.\n")

    # 6. FastTree Generation
    if n >= 10000:
        append_and_print_message(log_file, f"Running FastTree (fastest settings) on {n} sequences...")
        fasttree_cline = FastTreeCommandline(nt=True, fastest=True, cat=20, input=maffted_database, out=tree_string)
    else:
        append_and_print_message(log_file, f"Running FastTree (standard settings) on {n} sequences...")
        fasttree_cline = FastTreeCommandline(nt=True, input=maffted_database, out=tree_string)

    if os.path.exists(maffted_database):
        try:
            stdout, stderr = fasttree_cline()
            if stderr:
                print(f"FastTree running messages: {stderr.strip()}")
        except Exception as e:
            print(f"Error when running FastTree: {e}")

    # 7. PhyloTree / ete3 Evaluation
    tree = PhyloTree(tree_string)
    species_leaves_dict = defaultdict(list)
    pattern = r'\.\d+_|[\d.]+'

    for leaf in tree:
        species_name = get_species_name(leaf, pattern)
        species_leaves_dict[species_name].append(leaf)

    monophyletic_species = []
    paraphyletic_species = []
    
    for species_name, species_leaves in species_leaves_dict.items():
        if is_monophyletic(tree, species_leaves):
            monophyletic_species.append(species_name)
        else:
            paraphyletic_species.append(species_name)

    monophyletic_species = sorted(monophyletic_species)
    with open(f"{run_name}_post_curation_monophyletic_group.txt", "w") as file:
        file.write("These species came out as monophyletic in the gene tree after curation.\n\n")
        for mono in monophyletic_species:
            file.write(f"{mono.replace('_', ' ')}\n")

    paraphyletic_species = sorted(paraphyletic_species)
    with open(f"{run_name}_post_curation_paraphyletic_group.txt", "w") as file:
        file.write("These species came out as paraphyletic in the gene tree after curation.\n\n")
        for para in paraphyletic_species:
            file.write(f"{para.replace('_', ' ')}\n")

    # Cleanup temporary files and return to base directory
    if os.path.exists(temp_file):
        os.remove(temp_file)
    os.chdir(base_directory)

    # 8. Archive old database
    if old_database:
        path_to_old_databases = "Old_databases/"
        os.makedirs(path_to_old_databases, exist_ok=True)
        try:
            os.rename(old_database, os.path.join(path_to_old_databases, os.path.basename(old_database)))
        except OSError as e:
            print(f"Notice: Could not move old database. {e}")

    # 9. Final Recommendations
    # 9. Final Recommendations
    program_duration = round(time.time() - program_timer, 2)
    relative_path_database = os.path.relpath(updated_database, base_directory)
    end_timestamp = datetime.today().strftime("%H:%M:%S")
    eval_monophyly_path = f"Database_curation/{run_name}/Curated_content/{run_name}_post_curation_monophyletic_group.txt"

    # Define the command string for the log_end helper
    next_cmd = (f"python echopipe.py evaluate {relative_path_database} {eval_monophyly_path} "
                f"[-f <forward_primer> -r <reverse_primer>]")

    append_and_print_message(log_file,
        f"\nIt took {program_duration} seconds from start to finish.\n"
        f"The program finished at {end_timestamp}.\n"
        f"The reference database {relative_path_database} has successfully been created!\n"
        "##############################################################################\n\n")

    next_cmd = f"python echopipe.py evaluate {relative_path_database} {eval_monophyly_path}"
    
    recommendation_msg = (
        "Conduct an additional evaluation of the database. See the tutorial for more instructions:\n"
        "https://github.com/EivindStensrud/EchoPipe/blob/main/Tutorial.md\n\n"
        "Evaluate the database (your original primers will be loaded automatically!):\n"
        f"\033[32m{next_cmd}\033[0m\n\n"
        "Want to test slightly different primer variations? Override the saved ones manually:\n"
        f"\033[32m{next_cmd} -f <variant_fwd> -r <variant_rev>\033[0m"
    )
    
    log_end(log_file, program_timer, recommendation_msg)

def run_evaluate(args):

    """Executes the Additional Evaluation logic (Dataframes & Primers)."""
    # 1. Setup paths and logging
    log_file, run_name = get_log_file_name()
    program_timer = time.time()

    log_start(log_file, "Database Evaluation")
    
    reference_database = os.path.abspath(args.reference_database)
    monophyletic_group = os.path.abspath(args.monophyletic_group)
    database_name = os.path.splitext(os.path.basename(args.reference_database))[0]
    evaluation_directory = os.path.join("Evaluation", database_name)

    config = load_config()
    forward_primer = args.forward_primer if args.forward_primer else config.get("forward", "")
    reverse_primer = args.reverse_primer if args.reverse_primer else config.get("reverse", "")

    os.makedirs(evaluation_directory, exist_ok=True)
    
    # 2. Logic Body
    with open(monophyletic_group, "r") as file:
        monophyletic_list = [line.strip() for line in file.readlines()]

    duplicate_sequences = find_duplicates(reference_database)

    # Dataframe structure definitions
    df_one_dict = {key: [] for key in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 
                                      'monophyletic', 'species_sequences', 'total_count', 'duplicate_sequences', 
                                      'average_gc_content', 'lowest_gc_content', 'highest_gc_content']}
    
    df_two_dict = {key: [] for key in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 
                                      'accession_number', 'sequence', 'count', 'length', 'gc_content']}

    lineage_keys = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    species_gc_content, family_gc_content = defaultdict(list), defaultdict(list)

    # Parse FASTA
    for record in SeqIO.parse(reference_database, "fasta"):
        lineage = record.id.split("|")[2].split(";")
        species, family = lineage[6], lineage[4]
        counter = int(record.id.split("|")[3])
        gc_val = gc_fraction(record.seq) * 100

        species_gc_content[species].append(gc_val)
        family_gc_content[family].append(gc_val)

        if species not in df_one_dict['species']:
            dup_count = sum(1 for entries in duplicate_sequences.values() for entry in entries if species in entry)
            for key, val in zip(lineage_keys, lineage):
                df_one_dict[key].append(val)
            df_one_dict['monophyletic'].append("Yes" if species.replace('_', ' ') in monophyletic_list else "No")
            df_one_dict['species_sequences'].append(1)
            df_one_dict['total_count'].append(counter)
            df_one_dict['duplicate_sequences'].append(dup_count)
        else:
            idx = df_one_dict['species'].index(species)
            df_one_dict['total_count'][idx] += counter
            df_one_dict['species_sequences'][idx] += 1

        for key, val in zip(lineage_keys, lineage):
            df_two_dict[key].append(val)
        df_two_dict.update({
            'accession_number': df_two_dict.get('accession_number', []) + [record.id.split("|")[1]],
            'sequence': df_two_dict.get('sequence', []) + [str(record.seq)],
            'count': df_two_dict.get('count', []) + [counter],
            'length': df_two_dict.get('length', []) + [len(record.seq)],
            'gc_content': df_two_dict.get('gc_content', []) + [round(gc_val, 2)]
        })

    # Calculate the final GC content stats for each species
    for sp in df_one_dict['species']:
        gc_list = species_gc_content[sp]
        df_one_dict['average_gc_content'].append(round(sum(gc_list) / len(gc_list), 2))
        df_one_dict['lowest_gc_content'].append(round(min(gc_list), 2))
        df_one_dict['highest_gc_content'].append(round(max(gc_list), 2))

    # Save CSVs
    pd.DataFrame(df_one_dict).to_csv(os.path.join(evaluation_directory, f"{database_name}_species_summary.csv"), sep=";", index=False)
    pd.DataFrame(df_two_dict).to_csv(os.path.join(evaluation_directory, f"{database_name}_sequence_summary.csv"), sep=";", index=False)

    # Primer Alignment Logic (Keep your existing primer_alignment loop here)
    if forward_primer or reverse_primer: # This blocked is skipped entirely if neither -f nor -r is used.
        append_and_print_message(log_file, "Primer compatibility check has been performed. Please review the relevant files.\n")

        if forward_primer:
            forward_max = int(len(forward_primer) * 1.2)
        if reverse_primer:
            reverse_max = int(len(reverse_primer) * 1.2)

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

        # Define mismatch sets for the different rules (based on Stadhouders el al., 2010).
        purine_purine_mismatches = {('A', 'T'), ('A', 'G'), ('G', 'A'), ('G', 'C'), ('C', 'G')}
        pyrimidine_pyrmidine_mismatches = {('T', 'A'), ('T', 'C'), ('C', 'T')}
        purine_pyrimidine_mismatches = {('C', 'A'), ('A', 'C'), ('G', 'T'), ('T', 'G')}

        filtered_records_perfect = []
        filtered_records_ok = []
        filtered_records_bad = []
        filtered_records_no_data = []
        filtered_records_all_data = []

        for record in primer_alignment_results:
            avoid_forward = False
            dodge_forward = False
            avoid_reverse = False
            dodge_reverse = False

            if forward_primer:
                forward_mismatch = record['forward_mismatches']
                forward_mismatches = record['forward_number_of_mismatches']
                forward_sequence = record['forward_sequence']
                avoid_forward = has_gap(forward_mismatch)

            if reverse_primer:
                reverse_mismatch = record['reverse_mismatches']
                reverse_mismatches = record['reverse_number_of_mismatches']
                reverse_sequence = record['reverse_sequence']
                avoid_reverse = has_gap(reverse_mismatch)

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
                if avoid_forward:
                    forward_sequence_status = "Bad"
                elif dodge_forward:
                    forward_sequence_status = "Ok"
                elif forward_sequence == "NA":
                    forward_sequence_status = "NA"
                else:
                    forward_sequence_status = "Perfect"

                record['forward_sequence_status'] = forward_sequence_status


            if reverse_primer:
                if avoid_reverse:
                    reverse_sequence_status = "Bad"
                elif dodge_reverse:
                    reverse_sequence_status = "Ok"
                elif reverse_sequence == "NA":
                    reverse_sequence_status = "NA"
                else:
                    reverse_sequence_status = "Perfect"

                record['reverse_sequence_status'] = reverse_sequence_status
                
            if forward_primer and reverse_primer: # This is triggered when both primers have been used as inputs.
                if forward_sequence_status == "Perfect" and reverse_sequence_status == "Perfect":
                    filtered_records_perfect.append(record)
                elif avoid_forward or avoid_reverse:
                    filtered_records_bad.append(record)
                elif dodge_forward or dodge_reverse:
                    filtered_records_ok.append(record)
                elif forward_sequence == "NA" and reverse_sequence == "NA":
                    filtered_records_no_data.append(record)

            elif forward_primer or reverse_primer: # If only one primer is used, this happens instead. Ensures functionality in both scenarios.
                if forward_sequence_status == "Perfect" and reverse_sequence_status == "Perfect":
                        filtered_records_perfect.append(record)
                elif avoid_forward or avoid_reverse:
                    filtered_records_bad.append(record)
                elif dodge_forward or dodge_reverse:
                    filtered_records_ok.append(record)
                elif forward_sequence == "NA" and reverse_sequence == "NA":
                    filtered_records_no_data.append(record)
                
        df_perfect = pd.DataFrame(filtered_records_perfect)
        df_ok = pd.DataFrame(filtered_records_ok)
        df_bad = pd.DataFrame(filtered_records_bad)
        df_no_data = pd.DataFrame(filtered_records_no_data)
        df_combined = pd.concat([df_perfect, df_ok, df_bad, df_no_data], ignore_index=True)
        
        df_perfect.to_csv(f"{database_name}_primer_result_perfect_entries.csv", sep=";", index=False)
        df_ok.to_csv(f"{database_name}_primer_result_ok_entries.csv", sep=";", index=False)
        df_bad.to_csv(f"{database_name}_primer_result_bad_entries.csv", sep=";", index=False)
        df_no_data.to_csv(f"{database_name}_primer_result_no_data_entries.csv", sep=";", index=False)
        df_combined.to_csv(f"{database_name}_primer_result_all_entries_info.csv", sep=";", index=False)


    if forward_primer:
        analyze_primer_frequencies(df_combined, "forward_sequence", forward_primer)
    if reverse_primer:
        analyze_primer_frequencies(df_combined, "reverse_sequence", reverse_primer)


    os.chdir("../..")

    # 3. Finalize and Print Status
    current_working_directory = os.getcwd()
    relative_path_database = os.path.relpath(evaluation_directory, current_working_directory)
    
    # Terminal prints
    print(f"\nThe reference database has successfully been evaluated!")
    print(f"\nDataframes have been created and saved to: {evaluation_directory}")
    
    if args.forward_primer or args.reverse_primer:
        print("Primer compatibility check has been performed. Please review the relevant files.")

    relative_ref_db = os.path.relpath(reference_database)

    # Suggestion for reformatting
    print("\n------------------------------------------------------------------------------")
    print("Suggestion:")
    print("If you need to prepare this database for downstream taxonomic classifiers (like SINTAX, QIIME 2, or DADA2),")
    print("you can optionally reformat the headers using the new reformat command:")
    print(f"\033[32mpython echopipe.py reformat {relative_ref_db} <format>\033[0m")
    print("Available formats: sintax, rdp, dadt, dads, idt, qiime")
    print("------------------------------------------------------------------------------\n")

    # 4. Final Logging
    next_step_msg = (f"Evaluation complete. Dataframes are located in: {evaluation_directory}\n"
                     f"Optional reformatting suggestion provided for taxonomic classifiers.")
    
    log_end(log_file, program_timer, next_step_msg)

def run_reformat(args):
    """Executes the Database Reformatting logic for taxonomy classifiers."""
    log_file, run_name = get_log_file_name()
    program_timer = time.time()
    
    log_start(log_file, "Database Reformatting")
    
    reference_database = os.path.abspath(args.reference_database)
    fmt = args.format
    
    database_name = os.path.splitext(reference_database)[0]
    reformated_database = f"{database_name}_{fmt}.fasta"

    def write_new_fasta(header, seq, file_handle):
        file_handle.write(header + "\n")
        file_handle.write(seq + "\n")

    append_and_print_message(log_file, f"\nRunning database reformatting to '{fmt}' format...\nCommand: {get_command_string()}\n")

    with open(reformated_database, "w") as file:
        for record in SeqIO.parse(reference_database, "fasta"):
            split_header = record.id.split("|")
            accession_number = split_header[1]
            lineage = split_header[2].split(";")
            sequence = str(record.seq)
            
            if fmt == "sintax":            
                sintax_header = f">{accession_number};tax=d:{lineage[0]},p:{lineage[1]},c:{lineage[2]},o:{lineage[3]},f:{lineage[4]},g:{lineage[5]},s:{lineage[6]}"
                write_new_fasta(sintax_header, sequence, file)

            elif fmt == "rdp":
                rdp_header = f">{accession_number}\troot;{lineage[0]};{lineage[1]};{lineage[2]};{lineage[3]};{lineage[4]};{lineage[5]};{lineage[6]}"
                write_new_fasta(rdp_header, sequence, file)

            elif fmt == "dadt":
                dad_header = f">{lineage[0]};{lineage[1]};{lineage[2]};{lineage[3]};{lineage[4]};{lineage[5]}"
                write_new_fasta(dad_header, sequence, file)

            elif fmt == "dads":
                dads_header = f">{accession_number} {lineage[5]} {lineage[6]}"
                write_new_fasta(dads_header, sequence, file)

            elif fmt == "idt":
                idt_header = f">Root;{lineage[0]};{lineage[1]};{lineage[2]};{lineage[3]};{lineage[4]};{lineage[5]};{lineage[6]}"
                write_new_fasta(idt_header, sequence, file)

            elif fmt == "qiime":
                qiif_header = f">{accession_number}"
                write_new_fasta(qiif_header, sequence, file)
                
    if fmt == "qiime": # QIIME 2 also requires a companion taxonomy text file.
        qiif_txt_file = f"{database_name}_{fmt}.txt"
        with open(qiif_txt_file, "w") as txt_handle:
            for record in SeqIO.parse(reference_database, "fasta"):
                split_header = record.id.split("|")
                accession_number = split_header[1]
                lineage = split_header[2].split(";")
                species = "_".join(lineage[6].split("_")[1:])
                qiif_txt_header = f"{accession_number}\tk__{lineage[0]}; p__{lineage[1]}; c__{lineage[2]}; o__{lineage[3]}; f__{lineage[4]}; g__{lineage[5]}; s__{species}"
                txt_handle.write(qiif_txt_header + "\n")
        print(f"QIIME companion taxonomy file saved to: {qiif_txt_file}")

    print(f"\nReformatting completed. Saved as: {reformated_database}\n")
    
    next_step_msg = f"Database successfully reformatted to {fmt}. Output file: {reformated_database}"
    log_end(log_file, program_timer, next_step_msg)

# =============================================================================
# COMMAND LINE INTERFACE (CLI) ROUTER
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        prog="echopipe",
        description="EchoPipe: A complete pipeline for reference database creation and curation.",
        epilog="Example usage: python echopipe.py create --help")
    
    subparsers = parser.add_subparsers(title="commands", dest="command", required=True)

    # 1. Template Parser
    parser_template = subparsers.add_parser("template", help="Generate a template reference database.")
    gen_args = parser_template.add_argument_group("Uncurated template generation arguments")

    gen_args.add_argument('input_file', type=str, nargs='?', help="Txt or CSV file species names or a fasta file.")
    gen_args.add_argument('-f', '--forward', type=str, default="", help="The forward primer used to find region of interest, (5'-3').")
    gen_args.add_argument('-r', '--reverse', type=str, default="", help="The reverse primer used to find region of interest, (5'-3').")
    gen_args.add_argument('-e', '--email', type=str, default="email@email.com", help="Your email if NCllBI needs to contact you.")
    gen_args.add_argument('-a', '--api_key', type=str, default="api_key", help="The user's NCBI API key.")
    gen_args.add_argument('-q', '--query', help='Custom query additions.')
    gen_args.add_argument('-t', '--threshold', type=int, default=150, help="The minimum length of a sequence.")
    gen_args.add_argument('-l', '--length', type=int, default=22000, help="The longest allowed sequence length.")
    gen_args.add_argument('-m', '--max', type=int, default=1, help="Number of sequences downloaded per species. Default = 1.")
    gen_args.add_argument('-p', '--provided_sequences', action='store_true', default='', help="Use a fasta file as reference template.")
    gen_args.add_argument('-z', '--longest_amplicon_size', type=float, default=2, help="Multiplier for median length.")
    gen_args.add_argument("-n", "--random_subset", type=int, help="Number of random species to use.")
    gen_args.add_argument("-sf", "--subset_file", type=str, help="Path to a file containing a specific subset of species.")
    gen_args.add_argument('-T', '--threads', type=int, default=None, help="Number of parallel threads to use (default: auto-detected, max 7).")
    
    comp = parser_template.add_argument_group("Argument used to finish the curated reference template database.")
    comp.add_argument('-C', '--Complete', action="store_true", help="Completes the reference template database.")
    comp.add_argument('input_file_species', type=str, nargs='?', help="Txt or CSV file species names.")
    
    # 2. Create Parser
    parser_create = subparsers.add_parser(
        "create", 
        help="Mine reference sequences using BLAST.",
        description="This command mines NCBI for reference sequences and creates a BLAST-ready database.") 

    create_args = parser_create.add_argument_group("Database creation arguments")
    create_args.add_argument('input_file', type=str, help="A txt file or CSV with a list of species names.")
    create_args.add_argument('input_database', type=str, help="Path to the input reference database fasta file.")
    create_args.add_argument('-e', '--email', type=str, default="email@email.com", help="User's email address.")
    create_args.add_argument('-a', '--api_key', type=str, default="api_key", help="User's NCBI API key.")
    create_args.add_argument('-s', '--sort', action="store_true", help="Sort by length (Not recommended).")
    create_args.add_argument('-c', '--maxcount', type=int, default=10000, help="Maximum accession numbers per species.")
    create_args.add_argument('-l', '--maxlength', type=int, default=22000, help="Longest allowed sequence length.")
    create_args.add_argument('-z', '--ampliconsize', type=int, default=50, help="Minimum size an amplicon may be.")
    create_args.add_argument('-m', '--mitochondria', action='store_true', help="Search targets mitochondrial sequences.")
    create_args.add_argument('-r', '--ribosomal', action='store_true', help="Search for mitochondrial 12S ribosomal DNA.")
    create_args.add_argument('-q', '--query', help='Custom NCBI search term.')
    create_args.add_argument('-b', '--batch_size', type=int, default=5000, help="Batch size for downloading sequences.")
    create_args.add_argument('-t', '--taxid', action='store_true', help="Use last saved taxid list.")
    create_args.add_argument('-E', '--evalue', type=int, default=20, help="E-value for BLAST. Default = 20. Increase for longer markers; keep lower for shorter markers to avoid introducing non-target gene regions.")
    create_args.add_argument('-R', '--repeat', action='store_true', help="Repeat curation on previously downloaded sequences.")
    create_args.add_argument('-T', '--threads', type=int, default=None, help="Number of parallel threads to use (default: auto-detected, max 7).")

    # 3. Curate Parser
    parser_curate = subparsers.add_parser("curate", 
        help="Align and generate trees for manual curation.",
        description="This command help with the curation of the database.")

    curate_args = parser_curate.add_argument_group("Curation arguments")
    curate_args.add_argument('input_file', type=str, help="Database to revise.")
    curate_args.add_argument('-o', '--old_database', type=str, default="", help="The previous version of database.")
    curate_args.add_argument('-N', '--number_ns', type=int, default=0, help="Number of N's and ambiguous nucleotides allowed.")
    curate_args.add_argument('-M', '--mafft_online', type=str, default="", help="Path to MAFFT online alignment file.")
    curate_args.add_argument('--min_length', type=int, default=150, help="Minimum sequence length to keep.")
    curate_args.add_argument('--max_length', type=int, default=float('inf'), help="Maximum sequence length to keep.")

    # 4. Complete Parser 
    parser_complete = subparsers.add_parser("complete",
        help="Filter, merge, and finalize the database.",
        description="This command filter, merges and finalizes the reference database.")

    parser_complete.add_argument('-b', '--blast_file', type=str, help="New database BLAST file.")
    parser_complete.add_argument('-c', '--curated_file', type=str, help="Curated aligned FASTA.")
    parser_complete.add_argument('-o', '--old_database', type=str, default="", help="Existing master database.")
    parser_complete.add_argument('-u', '--updated_database', type=str, default="database.fasta", help="Output database.")

    # 5. Evaluate Parser
    parser_evaluate = subparsers.add_parser("evaluate", 
        help="Evaluate the database (GC content, primers).",
        description="This command evaluates the reference database.")

    eval_args = parser_evaluate.add_argument_group("Evaluation arguments")
    eval_args.add_argument('reference_database', type=str, help="Path to the reference database.")
    eval_args.add_argument('monophyletic_group', type=str, help="Path to the monophyletic groups text file.")
    eval_args.add_argument('-f', '--forward_primer', type=str, default=None, help="The forward primer sequence to check (5'-3').")
    eval_args.add_argument('-r', '--reverse_primer', type=str, default=None, help="The reverse primer sequence to check (5'-3').")

    # 6. Reformat Parser
    parser_reformat = subparsers.add_parser(
        "reformat", 
        help="Reformat database headers for popular taxonomic classifiers.",
        description=(
            "Reformats the input database headers to a pre-defined format:\n"
            "sintax = SINTAX\n"
            "rdp = Ribosomal Database Project (RDP)\n"
            "dadt = DADA2 assignTaxonomy\n"
            "dads = DADA2 assignSpecies\n"
            "idt = IDTAXA\n"
            "qiime = QIIME 2\n"))

    parser_reformat.add_argument('reference_database', type=str, help="The reference database for which the header format will be changed.")
    parser_reformat.add_argument('format', type=str, choices=['sintax', 'rdp', 'dadt', 'dads', 'idt', 'qiime'], help="Write in one of the available formats.")

    args = parser.parse_args()

    commands = {
        "create": run_create,
        "template": run_template,
        "curate": run_curate,
        "complete": run_complete,
        "evaluate": run_evaluate,
        "reformat": run_reformat
    }
    commands[args.command](args)

if __name__ == "__main__":
    main()