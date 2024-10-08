#!/usr/bin/env python3

# ---------------------------------------------------------------------------
# Copyright (c) 2024, Daniel Borg, Eivind Stensrud, Alexander Eiler
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ---------------------------------------------------------------------------


"""
This pipeline mines potential reference sequences for eDNA metabarcoding studies using sequence similarity.
The required input files are a species list (both .txt and .csv works) and a few reference sequences of the genomic area of interest of related species (.fasta).
Additional inputs are the user's email address and their API key to NCBI, which are required for optimal performance.
"""
import warnings
from Bio import BiopythonDeprecationWarning
# Suppress the specific Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from Bio import Entrez # Used for working with and obtaining data from NCBI.
from Bio import SeqIO # Used for parsing sequence information.
from Bio.Blast.Applications import NcbimakeblastdbCommandline # Used to make our local database for the purpose of BLASTing.
from Bio.Blast.Applications import NcbiblastnCommandline # So that we can BLAST locally on the computer. NOTE: MUST DOWNLOAD BLAST!!!
# https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html check this out for downloading BLAST.
import os # In order to allow the creation of new directories/folders.
import os.path # Retrieves a file's name and/or path.
import time # In order to add delays, so that servers are not put under much strain (and if one would like to see how long it takes for the code to execute).
import pandas as pd # Used to make a dataframe. Columns. Dealing with BLAST output.
import sys # Allows for the program to be shut down. Executed if NCBI is unresponsive.
from tqdm import tqdm # For progress bar. Entirely optional. Note it may need to be installed!
from urllib.error import HTTPError # Is useful for handling potential HTTP errors.
from datetime import datetime # So that we can get the current day's date.
import argparse # Allows the program to be run from the command line.


# Input parameters
#
#
parser = argparse.ArgumentParser(
    prog="EchoPipe - Database creation",
    description="This script creates an uncurated reference database.",
    epilog="Version 1.0")

group_1 = parser.add_argument_group("Arguments used to generate the uncurated reference database.\nThe required arguments are input_file, input_database, -e (--email) and -a (--api_key)")
group_1.add_argument('input_file', type=str,
    help="A txt file or CSV with a list of species names.")
group_1.add_argument('input_database', type=str,
    help="Path to the input reference database fasta file.")
group_1.add_argument('-e', '--email', type=str, default="eivisten@uio.no",
    help="User's email address associated with the NCBI API key.")
group_1.add_argument('-a', '--api_key', type=str, default="5538bd02e0d704e3416263ad4d51c74b7608",
    help="User's NCBI API key.")
group_1.add_argument('-s', '--sort', action="store_true",
    help="Sort by length, if implemented, targets longer sequences (Not recommended).")
group_1.add_argument('-c', '--maxcount', type=int, default=10000,
    help="Maximum number of accession numbers downloaded per species, default = 10 000.")
group_1.add_argument('-l', '--maxlength', type=int, default=22000,
    help="Longest allowed sequence length for an accession number to be analysed, default = 22 000.")
group_1.add_argument('-z', '--ampliconsize', type=int, default=50,
    help="Minimum size an amplicon may be in order to accepted. Consider adjusting this based on the marker region. Default = 50.")
group_1.add_argument('-m', '--mitochondria', action='store_true',
    help="Search targets mitochondrial sequences.")
group_1.add_argument('-r', '--ribosomal', action='store_true',
    help="Search result will only return accession numbers that have been annotated with mitochondrial 12S ribosomal DNA.")
group_1.add_argument('-q', '--query', action='store_true',
    help='Search result will include the user input search term(s). Example, limit the search to 12S region: -q "AND 12S". To exclude a term write "NOT 12S".')
group_1.add_argument('-b', '--batch_size', type=int, default=5000,
    help="Amount of sequences that can be downloaded simultaneously. Default = 5000. Only need to consider using this with values lower than the default and only use it in conjunction with --maxlength and expecting to download very large sequences (size of chromosomes).")
group_1.add_argument('-t', '--taxid', action='store_true', default="",
    help="Use last saved taxid list. Can be used to change NCBI taxonomy (Not recommended)")
group_1.add_argument('-E', '--evalue', type=int, default=5,
    help="The E-value used for BLAST. The higher the value entered, the more stringent the BLAST becomes. Default = 5.")
group_1.add_argument('-R', '--repeat', action='store_true', default="",
    help="No new sequences are downloaded. Instead, all previously downloaded sequences are BLASTed against the input database. Note: this options resets the file containing analysed accession numbers and the output will be curated from the start.")

args = parser.parse_args()

if not args.repeat:
    if not all ([args.input_file, args.input_database, args.email, args.api_key]):
        parser.error("An input file with a list of species, input database, email address and API key are required.\n(Email address and API key are not needed when re-analyzing sequences with -R, --repeat).")

# Arguments from argparse are stored as variables used downstream.
input_species = args.input_file
input_fasta = args.input_database
Entrez.email = args.email
Entrez.api_key = args.api_key
sort_by_sequence_length = args.sort
max_search = args.maxcount
max_length_input = args.maxlength
min_length = args.ampliconsize
max_length = f' AND ("20"[SLEN] : "{str(max_length_input)}"[SLEN])'
mitochondria_on = ' AND mitochondrion[filter]' if args.mitochondria else ""
ribosomal_on = ' AND 12S' if args.ribosomal else ""
custom_query = f' {args.query}' if args.query else ""
sequence_batch_size = args.batch_size
use_old_taxid = args.taxid
evalue = f"5e-{args.evalue}"
reuse_sequences = args.repeat


DATE = datetime.today().strftime('%Y-%m-%d')
i = 1
run_name = f"{DATE}_{i}"
search_settings = '[Organism] AND biomol_genomic[PROP]' + max_length + mitochondria_on + ribosomal_on + custom_query # Baseline and complimentary search terms. Can, for example, filter away all non-annotated mitochondrial sequences.

program_timer = time.time()

###
### Functions below
###

# Adds the new accession numbers from run_name to species.analysed_accession_numbers.txt.
def analysis_completed(finished_accession_numbers, analysed_accession_numbers):
    if reuse_sequences: # This is checked first in order to get the accession file cleared.
        with open(analysed_accession_numbers, "w") as file:
            file.truncate(0) # The content of the file is reduced to 0 bytes, i.e., it is made empty.

    with open(finished_accession_numbers, "r") as file: # Opens up and reads the accession numbers file...
        finished_accession_content = set(file.readlines()) # ... and stores them as a set in finished_accession_content.
    with open(analysed_accession_numbers, "r") as file: # Does the same thing as above, but with the analysed accession numbers.
        previous_accession_content = set(file.readlines())
    new_additions = finished_accession_content.difference(previous_accession_content) # Accession numbers that are not in previous_accession_content are stored in "new_additions".
    
    with open(analysed_accession_numbers, "a") as file: # Appends the content from the variable (i.e., the finished accession numbers) onto species.analysed_accession_numbers.txt.
        for line in new_additions:
            file.write(line)
    with open(analysed_accession_numbers, "r") as file: # Opens, reads and sorts the accession numbers and saves them to a new variable.
        sorted_accession_numbers = file.readlines()
        sorted_accession_numbers.sort()
    with open(analysed_accession_numbers, "w") as file: # Writes the now sorted accession numbers to species_analysed_accession_numbers.txt.
        file.writelines(sorted_accession_numbers)

# Used to concatenate files when re-analysing previous sequence files.
def concatenate_files(input_list, concatenated_file):
    with open(concatenated_file, "w") as outfile:
        for previous_file in input_list:
            with open(previous_file, "r") as infile:
                for line in infile:
                    outfile.write(line)

# Used to make a data frame from blast_output.
def columnize(blast_output):
    data_frame = pd.read_table(blast_output, header=None)
    data_frame.columns = ["qseqid", "qseq", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return data_frame

# Used to create a fasta file named after blast_filename with the desired header for each respective sequence.
def write_fasta(data, blast_filename):
    with open(blast_filename, "w") as file_conn:
        for _, row in data.iterrows():
            header = f'>{row["qseqid"]}|{row["identical_seqs"]}\n'
            sequence = str(row["qseq"])
            file_conn.write(header)
            file_conn.write(sequence + '\n')

# Used to create directories if they do not already exist.
def create_directory(name_of_directory):
    if not os.path.exists(name_of_directory): # This part checks if it exists...
        os.makedirs(name_of_directory) # ... and if it does not, the directory with the input name is created.

# This checks if a file exists, and creates it if it does not. Useful on certain locations in the code where a file is not created on the spot, but is still expected to read one.
def does_file_exist(file_to_check_path):
    check_file = os.path.isfile(file_to_check_path) # If the defined path does not exist, then check_file = False.
    if check_file is False: #
        with open(file_to_check_path, "w") as file: # This, as well with the next line below creates the file...
            pass # ... and the file is empty. Note: `pass` could be replaced with file.write("") for the same result! Just a matter of user readability.

# Update the last line (in the log file) in an ongoing loop.
def update_last_line(file_path, progress_info):
    with open(file_path, 'r') as file:
        lines = file.readlines() # Attempts to read the lines of the file...
    if lines and not idx==1: # ... and if they file contains lines, then the "if lines" is True. However, if idx==1 we want to skip this part, as we want to append the first line and not overwrite anything.
        lines[-1] = progress_info + '\n' # This overwrites the current last line with progress_info and a newline (\n).
    else: # This occurs if either the file has no lines or if idx==1, which it should be on the first iteration of the loop.
        lines.append(progress_info + '\n') # This line of data is appended, rather than overwriting another line.
    with open(file_path, 'w') as file:
        file.writelines(lines) # The line gets actually written on to the file.


# A function that prints a message (msg) and appends it on to a file.
def append_and_print_message(log_file, msg):
    with open(log_file, "a") as file:
        file.write(msg)
    print(msg)


# Shuts down the program if an error persists. Could be problems with NCBI servers and not worth wasting time trying.
def repeated_failures(consecutive_fail_counter, log_file, check_for_fail, type_of_download):
    if consecutive_fail_counter == 2:
        message = "\nTwo consecutive download errors encountered. Pausing for 10 minutes before reattempting.\n"
        append_and_print_message(log_file, message) # Runs the function that prints the message and writes it to the log_file.
        time.sleep(600) # 10 minute pause before trying once more.
    elif consecutive_fail_counter == 3:
        message = "\nErrors persist. Shutting down the program. Try again later."
        if check_for_fail:
            with open(log_file, "a") as file:
                file.write(f'\n{type_of_download} is finished. Please note:\n')
                file.write('\n'.join(check_for_fail) + '\n')
        append_and_print_message(log_file, message)
        sys.exit() # Kills the program.

# Appends a message to the log_file stating a loop (type_of_download) is finished and how long it took.
def loop_finished(check_for_fail, log_file, type_of_download):
    with open(log_file, "a") as file:
        if check_for_fail:
            file.write(f'\n{type_of_download} is finished. Please note:\n')
            file.write('\n'.join(check_for_fail) + '\n')
        else:
            file.write(f'\n{type_of_download} is finished! No errors detected.\n'
                f'It took {loop_duration} seconds to finish.\n\n'
                '--------------------------------------------------------------------------\n\n\n')

# Writes lines from a list to a file.
def write_lines_from_list(list_with_lines, lines_to_file_path, log_file=None, msg=None): # log_file=None and msg=None makes these optional arguments.
    #if list_with_lines: # First checks if the list has any content. Otherwise an empty file would be created.
        with open(lines_to_file_path, "w") as file: # Opens the file in write-mode.
            for line in list_with_lines:
                file.writelines(line) # Writes all of the lines from the list on to the file.
        if msg is not None: # If msg is set to be a string, it is used and calls the function append_and_print_message.
            append_and_print_message(log_file, msg) # Uses the argument to msg and appends it to whatever log_file is defined as when the function write_lines_from_list is called.

# Used when errors may occur. Prints what the error is called, for what species and also announces how long of a pause is initiated.
def error_occurrence(error):
    global retry_delay
    print(f"An error occurred: {e}") # Let's us know an error has occurred. Not as specific as the HTTP one above, but this covers for less common errors as well.
    print(f"Error occurred for {species}")
    print(f"Initiating a {retry_delay} seconds pause before retrying.")

    time.sleep(retry_delay) # Pause if there is an error. Repeats the amount of time max_retries is set to. See below for duration.
    retry_delay *=2 # Exponentially increases retry_delay. If max_retries=3 and retry_delay is initially set to 15, then it increases to 30 and finally 60 for the last attempt.

    return error # This return statement allows "error" to be saved from the function. saved_variable = error_occurrence(e) means that "e" is saved in saved_variable.

###
###
###


database_directory = "Template_databases/" + os.path.splitext(input_fasta)[0] # Where the database should be stored. With the placeholder used above, the location will be Template_databases/Swedish_amphibians_ref_db

logs = "Log_files/" # Directory where files containing information from the search is stored. Information such as if the search had duplicate species (and which), what parametres were set (e.g. maximum amount of accession numbers).
log_file = f"{logs}{run_name}_log.txt" # Log file where both completed and ongoing progress is stored.

taxid_collection = "taxid_collection.txt" # Species' unique taxids are stored here and used for retrieving data.
taxid_log_file = f"{logs}{run_name}_{taxid_collection}" # A copy of the taxid_collection.txt file that is either re-used or created. This file, however, is placed in Log_files/.
error_collection = f"{logs}species_not_found.txt" # If a species is not found from input_species, it will be listed here.
duplicate_collection = f"{logs}duplicate_species_entries.txt" # If there are duplicate entries of the same species in input_species, it will be listed here (this is based on taxid, so species with more than one name may be detected).

taxonomic_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus'] # Could be modified by the user. Used to make a FASTA header downstreams. Note: species are included downstreams.

###
###
###


create_directory(database_directory) # Creates a directory for the database to be used.

update_run_name = False # Makes sure update_run_name is flagged as False.
# Not the most elegant while loop. Got the job done, though.
while os.path.exists(log_file): # If log_file exists, then it means the program already has been ran with the current run_name. A new variable is assigned to run_name down below based on how many runs has been made on the initial DATE.
    if update_run_name:
        run_name = f"{DATE}_{i}" # run_name is updated if update_run_name is True. Updated run_name could look like 2023_12_24_2, where the _2 simply indicates it is the second run of the day.
    i += 1 # Since i is set to 1 before entering the loop, it will always turn into 2 immediately if the loop begins.
    update_run_name = True # Flags for run_name to be updated with an extension based on the value of i.
    log_file = f"{logs}{run_name}_log.txt"


create_directory(logs) # Creates a directory where non-species specific log files are stored.
timestamp = datetime.today().strftime("%H:%M:%S") # Creates a variable called timestamp with the format HH:MM:SS.
with open(log_file, "w") as file:
    file.write(
        f"This log file contains information from the run {run_name} at {timestamp}.\n\n"
        "--------------------------------------------------------------------------\n\n\n"
        "The script was run with the following settings:\n"
        f"Input file: {args.input_file}\n"
        f"Input reference database: {args.input_database}\n"
        f"Email address: {args.email}\n"
        f"API key: {args.api_key}\n"
        f"Date: {DATE}\n")
    if not reuse_sequences:
        file.write(
            f"Sorted by largest sequence lengths: {args.sort}\n"
            f"Maximum accession number per species: {args.maxcount}\n"
            f"Minimum amplicon size: {args.ampliconsize}\n"
            f"Search term: {search_settings}\n\n"
            "--------------------------------------------------------------------------\n\n")
    else: # I.e., if reuse_sequences is True.
        file.write(
            f"Existing sequence files are BLASTed against the input database.\n\n"
            "--------------------------------------------------------------------------\n\n"
        )

database_name = os.path.splitext(input_fasta)[0] # The name of the database. Can be something else. The code removes the extension from a file and makes it the second object [1]. Example: [0] = Swedish_amphibians_ref_db
database_path = os.path.join(database_directory, database_name) # Basically where the database will be. database_path will be used when performing BLASTs.

makeblastdb_cline = NcbimakeblastdbCommandline( # This will create the database based on the input done by the user (input_fasta).
    dbtype="nucl", # Nucleotide, as opposed to protein ("prot").
    input_file=input_fasta, # The user's input.
    out=database_path, # The output. How we make use of the database.
    title=input_fasta, # Optional.
)
stdout, stderr = makeblastdb_cline() # This executes the makeblastdb command.
if stderr:
    append_and_print_message(log_file,
        f"Error: {stderr}\n")
    #print("Error:", stderr) # This was the old code. Kept here for reference in case the changed code fails.
else:
    append_and_print_message(log_file,
        f"Template database creation completed successfully in directory: {database_directory}\n\n"
        "--------------------------------------------------------------------------\n\n")
    #print(f"Database creation completed successfully in directory: {database_directory}") # This was the old code. Kept here for reference in case the changed code fails.


###
###
###

consecutive_fail_counter = 0 # If it reaches its threshold, a 10 minutes pause is initiated. If subsequent failure occurrs (e.g., NCBI not responding), then the program shuts down.
max_retries = 4 # The amount of time we wish to retry to download the data from NCBI.


if not use_old_taxid: # Skips this step if the user decided to go with the previous taxid collection - saves time proportional to the number of input species.
    type_of_download = "Species information download" # This block downloads species information, such as taxid.
    append_and_print_message(log_file, f"{type_of_download} is starting...\n")

    with open(input_species, "r") as file: # Takes the species the user has listed in the file...
        input_species_list = file.readlines() # ... reads the list and stores it in input_species_list.
    input_species_list = sorted([line.strip().split(";")[0] for line in input_species_list if line.strip()]) # input_species_list is sorted and unwanted white space are removed.

    species_not_found = [] # Empty lists are created that are to be filled with the correct corresponding species to later on be written to the appropriate file.
    duplicate_species = []
    filtered_species_collection = []
    failed_information_downloads = [] # If there are any species that fails to have their information downloaded, they are stored here and will be written to the log_file.


    loop_timer = time.time() # Used to calculate the duration of the loop. At this location it calculates the time it takes to download all accession numbers.

    ### I think the new stuff have been covered here...
    ###

    for idx, species in tqdm(enumerate(input_species_list, start=1), total=len(input_species_list), desc="Going through the species list"): # Loops through each species that originates from the users input.
        retry_delay = 10 # Delay in seconds between each new attempt to download data after an error.

        repeated_failures(consecutive_fail_counter, log_file, failed_information_downloads, type_of_download)
        
        for _ in range(max_retries):
            try:
                taxid_handle = Entrez.esearch( # Performs an esearch where we are interested in the taxids from the species in the list.
                    db="Taxonomy", 
                    term=species, # Our search term.
                    retmode="xml",
                    usehistory="y"
                )
                taxid_record = Entrez.read(taxid_handle) # Reads the data from the esearch and stores it in taxid_record.
                taxid_handle.close() # Close the handle after we are done with it.
                taxid = taxid_record["IdList"] # "IdList" has the taxid here, since the search was done in the taxonomy database. With this, the current species' taxid is store in the variable called taxid.

                time.sleep(0.5) # A brief pause between the esearch and the efetch.

                if taxid_record["Count"] == str(0) or "ErrorList" in taxid_record: # When this is True it means that the search (i.e. from a misspelled name or a non-existing one) came out as negative and no data was retrieved.
                    species_not_found.append(species + "\n") # The species' name is added to the species_not_found list.
                    break
                elif any(taxid[0] == line.split(';')[1].strip() for line in filtered_species_collection): # Checks to see if the taxid already exists in the list filtered_species_collection.
                    duplicate_species.append(species + "\n") # If it does exist, it means it is a duplicate entry and is thus added to the duplicate_species list.
                    break
                else: # If the taxid for the current species exists and is not a duplicate entry, the loop proceeds.
                    webenv = taxid_record["WebEnv"] # WebEnv is used for efetch in order to be nice towards NCBI by helping out with the upcoming efetch.
                    query_key = taxid_record["QueryKey"] # Same purpose as WebEnv.

                    fetch_handle = Entrez.efetch( # Performs an efetch on the Taxonomy database in order to get taxonimic information about the species.
                        db="Taxonomy",
                        retmode="xml",
                        webenv=webenv,
                        query_key=query_key
                    )
                    fetch_taxonomy_data = Entrez.read(fetch_handle) # Taxonomic data is stored in fetch_taxonomy_date.
                    fetch_handle.close() # Close the handle after we are done with it.

                    scientific_name = fetch_taxonomy_data[0]["ScientificName"] # The appropriate scientific name is used (regardless of the search input; old names may still be used to find results).
                    LineageEx = fetch_taxonomy_data[0]["LineageEx"] # Stores the lineage data in LineageEx.
                    taxa_rank_list = [] # Creates a list for taxonomic ranks to be added onto for each species. Used in a FASTA header downstreams.
                    for taxonomic_rank in taxonomic_ranks: 
                        taxonomic_rank_exist = False # Set to default. Used to see if a particular taxonomic rank we are interested is exists/is assigned to the species.

                        for item in LineageEx: # Loops through all the items in LineageEx, wherein "Rank" corresponds to the taxonomic ranks we are looking at. "ScientificName" represents the name of such ranks. E.g., for humans: genus = Homo.
                            if item["Rank"] == taxonomic_rank: # if item["Rank"] is the same as taxonomic_rank...
                                taxonomic_rank_exist = True # ... then it is a match; it exists...
                                taxa_rank_list.append(item["ScientificName"]) # ... and item["ScientificName"] (e.g. genus = Homo) is added to taxa_rank_list.

                        if not taxonomic_rank_exist: # If the taxonomic rank we are looking at exists/is assigned to the species...
                            taxa_rank_list.append("NA") # ... then the phrase NA is appended to the position that such a rank normally would occupy in the header.
                    taxa_rank_list.append(scientific_name.replace(' ', '_')) # Removes space between Genus species to get Genus_species to be used in the header.
                    taxa_to_fasta = ";".join(taxa_rank_list) # Joins each entry inside of taxa_rank_list together by ;. Example: Eukaryota;Chordata;Amphibia;Caudata;Salamandridae;Triturus;Triturus_cristatus

                    taxid_entry = f"{species}; " + f"{taxid[0]}; " + f"{scientific_name}; " + taxa_to_fasta # taxid_entry gets the name of the species that came from user input, taxid, scientific name (which can very well be the same as the input) and the species lineage we use for FASTA downstreams.
                    filtered_species_collection.append(taxid_entry + "\n") # taxid_entry is appended to the list filtered_species_collection.
                time.sleep(0.5)

                consecutive_fail_counter = 0 # Makes sure the timer is set to 0 after successful iterations.
                
                break # Is essential. Without this 'break', it will loop each species per max_retries. If a loop is successful, then it moves on to the next species.
            except HTTPError as e: # Specific for HTTP errors.
                last_exception = error_occurrence(e)

            except Exception as e: # Non-specified errors.
                last_exception = error_occurrence(e)
       
        else:
            print(f"{species} failed to download.")
            log_entry = f"Error: Failed to download information about {species} - {last_exception}"
            failed_information_downloads.append(log_entry)

            consecutive_fail_counter += 1
            
        log_progress = f"Progress: {idx} out of {len(input_species_list)} species from the input list has been processed. Most recent species: {species}."
        update_last_line(log_file, log_progress)

    loop_duration = round(time.time() - loop_timer, 2)

    loop_finished(failed_information_downloads, log_file, type_of_download)

    print(f"{type_of_download} took {loop_duration} seconds to finish.")


    if species_not_found:
        species_not_found.insert(0, f"The following species had no taxonomic search results with the input species list used on the run {run_name}.\n\n")
    write_lines_from_list(species_not_found, error_collection, log_file, # Calls the function write_lines_from_list with additional parameters that ensures a message is printed and appended to the log.
        f"Non-valid species name found. See {error_collection}\n") # This is the message that gets written to log_file.

    if duplicate_species:
        duplicate_species.insert(0, f"The following species were duplicates from the input species list from the run {run_name}.\n\n")
    write_lines_from_list(duplicate_species, duplicate_collection, log_file,
        f"Duplicate species name found. See {duplicate_collection}\n") # This is the message that gets written to log_file.

    write_lines_from_list(filtered_species_collection, taxid_collection) # Calls the function write_lines_from_list without additional parameters, thus it does not append to the log.

else:
    append_and_print_message(log_file,
        "Species list from a previous run is used.\n"
        "No new species information is downloaded.\n\n"
        "--------------------------------------------------------------------------\n\n\n")

###
###
###

with open(taxid_collection, "r") as file: # Opens taxid_collection in read mode.
    species_list = file.readlines() # Stores the content as species_list.
write_lines_from_list(species_list, taxid_log_file) # Creates an identical file from the content in taxid_collection, however, at a different location for logging purposes.
species_list = sorted([line.strip().split("; ")[2].replace(' ', '_') for line in species_list]) # species_list is slightly redefined. Index nr 2 [2] after splitting the "; ", which is the scientific name retrieved from NCBI, is formated as Homo_sapiens.

if not reuse_sequences:
    type_of_download = "Accession number download"
    append_and_print_message(log_file, f"{type_of_download} is starting...\n")
                             
    #max_search = 10000 # The maximum amount of accession numbers that will be downloaded for a species. Can be changed.


    over_cap_csv_list = [] # An empty list where species that exceeded the max_search count are added to, so users can be informed that the search does not account for all accession numbers.
    failed_accession_downloads = [] # If there are any species that fails to have their accession numbers downloaded, they are stored here and will be written to the log_file.

    create_directory("Species") # Creates a directory called Species (where all subdirectories for each species will be located) if the directory does not already exist.

    consecutive_fail_counter = 0
    failed_information_downloads = [] # If there are any species that fails to have their information downloaded, they are stored here and will be written to the log_file.

    batch_size = 10000 # 10000 is the maximum amount that can be retrieved at once.
    if max_search <= 10000: # If max_search is made to be lower than 10000, then batch_size must be adjusted as well.
        batch_size = max_search


    loop_timer = time.time() # Used to calculate the duration of the loop. At this location it calculates the time it takes to download all accession numbers.


    for idx, species in tqdm(enumerate(species_list, start=1), total=len(species_list), desc="Downloading accession numbers"):
        species_directory = f"Species/{species}" # Species' directory location from the working directory, e.g., Species/Alosa_alabamae
        create_directory(species_directory) # Creates a directory for each species (if such a directory does not already exist).
        
        search = str(species.replace("_", " ")+search_settings) # Defining the search term we want to use with esearch. E.g. Homo sapiens[Organism]. We also have implemented a sequence length (SLEN) limit to avoid crashes by ignoring massive sequences.
        
        retry_delay = 10

        repeated_failures(consecutive_fail_counter, log_file, failed_accession_downloads, type_of_download) # Checks if the program is struggling to perform successful downloads from NCBI. Here it downloads accession numbers.
        
        for _ in range(max_retries):
            try:
                accession_collection = [] # Creates an empty list where the accession numbers from each batch is to be stored before being written into a file.

                for start in range(0, max_search, batch_size):

                    search_handle = Entrez.esearch( # An esearch is performed...
                        db="nucleotide", #... on the nucleotide database.
                        term=search,
                        retmax=batch_size, # This number can be bypassed with retstart shennanigans.
                        retstart=start, # In order to bypass the 10000 limit of retmax.
                        sort="Sequence Length" if sort_by_sequence_length else None, # Sorts by sequence length instead of relevance/default. Bigger sequences comes first. Included as an option.
                        idtype="acc"
                    )
                    search_record = Entrez.read(search_handle) # Search data is stored in search_record.
                    search_handle.close()    
                    accession_collection.append(search_record["IdList"])

                    count = int(search_record["Count"]) # Makes an integer of the amount of accession numbers a species have.
                    if count > max_search: # This ensures that our desired max_search is upheld.
                        over_cap_csv_list.append(f"{species};{str(max_search)};{str(count)}\n")
                        count = max_search

                    time.sleep(0.5)

                    if count <= 10000: # If the existing amount of accession numbers is less than 10000 (the highest allowed for batch_size) then will not attempt to download more numbers when max_search is higher than 10000.
                        break

                consecutive_fail_counter = 0 # Makes sure the timer is set to 0 after successful iterations.
                
                break
                    
            except HTTPError as e:
                last_exception = error_occurrence(e)
                
            except Exception as e: # Non-specified errors.
                last_exception = error_occurrence(e)
        
        else:
            print(f"Unable to retrieve accession numbers for {species}.")
            log_entry = f"Error: Failed to download accession numbers for {species} - {last_exception}"
            failed_accession_downloads.append(log_entry)

            consecutive_fail_counter += 1
        
        with open(f"./{species_directory}/{species}_{run_name}_accession_numbers.txt", "w") as file: # File where the accession numbers from run_name are stored.
            for accession_batch in accession_collection:
                for accession_number in accession_batch:
                    file.write(accession_number + "\n")

        log_progress = f"Progress: {idx} out of {len(species_list)} species has been processed. Most recent species: {species.replace('_', ' ')}."
        update_last_line(log_file, log_progress)
        
        time.sleep(1) # Gives it a quick rest before the loop performs an esearch again (or code simply continues).

    loop_duration = round(time.time() - loop_timer, 2)

    loop_finished(failed_accession_downloads, log_file, type_of_download)

    write_lines_from_list(over_cap_csv_list, f"{run_name}_over_max.csv", log_file,
        f"Note: Not all sequences were downloaded for one or more species since --maxcount ({max_search}) < available accession numbers. See {run_name}_over_max.csv\n"
        f"Suggestion: Run program again with {run_name}_over_max.csv as input after database curation.\n"
        f"Make sure to adjust parameters such as --maxcount and/or --sort by sequence length.\n\n")

    accession_collection = [] # Makes the list empty again in order to release some of the memory. Could be smart in case the list holds a high amount of accession numbers.

    print(f"{type_of_download} took {loop_duration} seconds to finish.")

else: # This is triggered if the user runs the program with the flag -R, --repeat.
    append_and_print_message(log_file,
        "Existing sequences are to be re-analysed.\n"
        "No accession numbers are downloaded.\n\n"
        "--------------------------------------------------------------------------\n\n\n")

###
###
###

analysed = "analysed.accession_numbers" # Makes code a bit shorter whenever this text is to be used.
new_accession_species = [] # Creates an empty list to hold species where there has been additions of new sequences.

if not reuse_sequences:
    for species in species_list: # Loops through each species in order to compare the accession numbers from run_name and of those that have been previously analysed.
        accession_file_1 = f"Species/{species}/{species}_{run_name}_accession_numbers.txt" # The accession numbers from run_name.
        accession_file_2 = f"Species/{species}/{species}_{analysed}.txt" # The previously analysed accession numbers.
        does_file_exist(accession_file_2) # Creates the file (if it does not exist) we need to read and compare accession_file_1 to.
        
        with open(accession_file_1, "r") as file: # Opens up and reads the accession number from run_name...
            accession_file_content_1 = set(file.readlines()) # ... and stores them as a set in accession_file_content_1.
        with open(accession_file_2, "r") as file: # Does the same thing as above, but with the analysed accession numbers.
            accession_file_content_2 = set(file.readlines())

        differences = sorted(accession_file_content_1.difference(accession_file_content_2)) # Accession numbers that are not in accession_file_content_2 are stored in "differences" and are sorted.

        if not differences: # If nothing is stored in differences, then it means the recently downloaded accession numbers has been previouly analysed.
            print(f"No new accession numbers detected for \033[32m{species.replace('_', ' ')}\033[0m.")
        else:
            print(f"New accession numbers detected for \033[32m{species.replace('_', ' ')}\033[0m.") # If the content does not match, it means new accession numbers have been found.
            new_accession_species.append(species) # Adds the species to the new list that will be iterated through downstream, as there is no point in looping through species with no new sequences.
            new_accession = f"Species/{species}/{species}_new_accession_numbers_{run_name}.txt" # Defines the path to the file we want to store new accession numbers from.
            write_lines_from_list(differences, new_accession)


    append_and_print_message(log_file,
        f"{len(new_accession_species)} out of {len(species_list)} species have new accession numbers.\n\n"
        "--------------------------------------------------------------------------\n\n\n")

else: # This is triggered if the user runs the program with the flag -R, --repeat. This code checks that there has been previous runs performed for the species.
    for species in species_list:
        accession_numbers_file = f"Species/{species}/{species}_{analysed}.txt"
        previous_accession_numbers_path = f"Species/{species}"
        check_if_analysed_exist = os.path.isfile(accession_numbers_file) # If the defined path does not exist, then check_if_analysed_exist = False.
        previous_accession_numbers_exist = any("new_accession_numbers" in file for file in os.listdir(previous_accession_numbers_path))
        if not check_if_analysed_exist and not previous_accession_numbers_exist: # This parts checks if the current species has 1) previously analysed accession numbers 2) at least a file with accession numbers to be re-analysed.
            continue # If both of the above are false, then the loop for the current species stops and moves on to the next.
        else:
            new_accession_species.append(species) # Adds the species to the new list that will be iterated through downstream.

        
###
###
###

if not reuse_sequences:
    type_of_download = "Sequence download"
    append_and_print_message(log_file, f"{type_of_download} is starting...\n")

    failed_sequence_download = []

    consecutive_fail_counter = 0
    failed_information_downloads = [] # If there are any species that fails to have their information downloaded, they are stored here and will be written to the log_file.


    loop_timer = time.time()

    for idx, species in enumerate(new_accession_species, start=1): # Looping through each species in order to find new FASTA sequences from the run_name's accession numbers.
        
        with open(taxid_collection, "r") as file: # Opens taxid_collection...
            species_info = file.readlines() # ... and stores the content inside species_info.
            for line in species_info: # Loops through each line
                if species in line: # Checks to see if species (e.g. Triturus_cristatus) exists inside of each line as it loops through.
                    header_lineage = line.strip().split("; ")[3] # Header for the FASTA file. The line is split at "; " and [3] (fourth object, if you will) is stored in header_lineage.
                    break # Stops the current loop - no point in looping again once the species has been found.

        path_to_new_acc = f"Species/{species}/{species}_new_accession_numbers_{run_name}.txt" # Path to the file with the new accession numbers, whose sequences are to be retrieved.

        with open(path_to_new_acc, "r") as file:
            new_accession_numbers = file.readlines() # The species' (new) accession numbers are stored in new_accession_numbers.
        new_accession_numbers = [line.strip() for line in new_accession_numbers] # Removes "\n" from each line.

        retry_delay = 10

        repeated_failures(consecutive_fail_counter, log_file, failed_information_downloads, type_of_download)
        
        for _ in range(max_retries):
            try:
                with open(f"Species/{species}/{species}_new_seqs_{run_name}.fasta", "w") as file: # Creates and opens a file where the sequences from the (new) accession numbers are to be written into.
                    if sequence_batch_size > 5000: # If the user has input a higher value than 5000 it is set back to the default, as this value should only be lowered when dealing with very large sequences.
                        sequence_batch_size = 5000 # The amount of accession numbers, and thus sequences, we download per request. Can be 10 000 at most, but keeping it lower puts less stress on the machine.
                    for i in tqdm(range(0, len(new_accession_numbers), sequence_batch_size), desc=f"Downloading new sequence(s) for {species.replace('_', ' ')}"): # tqdm adds a progres bar. The loop itself is for iterating through the accession numbers in batches.
                        batch = new_accession_numbers[i:i + sequence_batch_size] # Updates the variable of batch for each loop in order to go through the next set of accession numbers.
                        with Entrez.efetch( # Performs an efetch to get the sequences.
                            db="nucleotide",
                            idtype="acc",
                            id=batch, # The accession numbers in batch are the search term - many at the same time!
                            rettype="fasta",
                            retmode="text"
                        ) as handle:
                            for seq_record in SeqIO.parse(handle, "fasta"): # The handle, i.e. the data from the search, is parsed through with SeqIO as seq_record.
                                file.write(str(f">gb|{seq_record.id}|{header_lineage}\n" + seq_record.seq + "\n")) # A header, newline, sequence, newline is written onto the file.
                        time.sleep(1) # Brief pause between batch search.
                print(f"New sequence(s) for \033[32m{species.replace('_', ' ')}\033[0m has been downloaded.") # Prints a message that new sequence(s) have been downloaded. Optional.

                consecutive_fail_counter = 0 # Makes sure the timer is set to 0 after successful iterations.
                
                break
            except HTTPError as e: # Specific for HTTP errors.
                last_exception = error_occurrence(e)
                
            except Exception as e: # Non-specified errors.
                last_exception = error_occurrence(e)
        else:
            print(f"Unable to download sequence(s) for {species}.")
            log_entry = f"Error: Failed to download information about {species} - {last_exception}"
            failed_information_downloads.append(log_entry)

            consecutive_fail_counter += 1

        log_progress = f"Progress {idx} out of {len(new_accession_species)} species. Most recent species: {species.replace('_', ' ')}."
        update_last_line(log_file, log_progress)

    loop_duration = round(time.time() - loop_timer, 2)

    loop_finished(failed_information_downloads, log_file, type_of_download)

    print(f"{type_of_download} took {loop_duration} seconds to finish.")

else: # This is triggered if the user runs the program with the flag -R, --repeat.
    append_and_print_message(log_file,
        "Existing sequences are to be re-analysed.\n"
        "No sequences are downloaded.\n\n"
        "--------------------------------------------------------------------------\n\n\n")


###
###
###


type_of_download = "BLAST"
append_and_print_message(log_file, f"{type_of_download} is starting...\n")

failed_blasts = []
nothing_new_to_blast = []


loop_timer = time.time()

outfmt_string = "6 qseqid qseq sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" # The format we want is nr 6, which is a tabular format. qseqid, qseq etc are the labels for the columns.

for idx, species in enumerate(new_accession_species, start=1):
    species_directory = f"Species/{species}"
    if reuse_sequences: # If previous sequence files are to be run, file_to_blast is made into a list in order to be able to hold multiple file names.
        blast_again = []
        for sequence_file in os.listdir(species_directory):
            if "new_seqs" in sequence_file:
                blast_again.append(os.path.join(species_directory, sequence_file))
        file_to_blast = f"{species_directory}/{species}_temp_concatenated_blast_file.fasta"
        concatenate_files(blast_again, file_to_blast)

        old_accession_collection = []
        for old_accession_file in os.listdir(species_directory):
            if "new_accession_numbers" in old_accession_file:
                old_accession_collection.append(os.path.join(species_directory, old_accession_file))
        finished_accession_numbers = f"{species_directory}/{species}_temp_concatenated_accession_file.txt"
        concatenate_files(old_accession_collection, finished_accession_numbers)

    else: # This happens during a normal run of the program; the newly downloaded sequence file is BLASTed.
        file_to_blast = f"{species_directory}/{species}_new_seqs_{run_name}.fasta" # We want to BLAST the current run_name's sequence(s).
        finished_accession_numbers = f"{species_directory}/{species}_{run_name}_accession_numbers.txt"
    blast_output = f"{species_directory}/{species}_{run_name}_fasta_blast_results.txt" # The name of the BLAST output.
    analysed_accession_numbers = f"{species_directory}/{species}_{analysed}.txt"

    # In order to check if the new_seqs file is empty. If it is empty, the file is removed and the current species' loop is stopped with "continue".
    if not reuse_sequences and os.path.isfile(file_to_blast) and os.path.getsize(file_to_blast) == 0: # Check if the file exists and is empty. Checking reuse_sequnces first will ensure the if statement is skipped if it is True.
        # Accession numbers from run_name are added to species_analysed_accession_numbers.txt in order to ensure the accession number is not used again.
        analysis_completed(finished_accession_numbers, analysed_accession_numbers)

        # Remove the file, since it is empty and unwanted.
        os.remove(file_to_blast)

        # Remove the name from the list.
        new_accession_species.remove(species)

        continue # Moves on to the next entry in new_accession_species.
    
    if reuse_sequences or os.path.exists(file_to_blast): # Proceeds with BLAST if a file based on the current run_name exists.
        blastn_cline = NcbiblastnCommandline(
            query = file_to_blast, # Query is thus the (new) sequences.
            db = database_path, # Specifying our database.
            out = blast_output, # The file that is modified in the next step. E.g., Lissotriton_vulgaris_2023-09-13_fasta_blast_results.txt
            max_target_seqs = 100, # The limit to our output. When set to 10, we get the 10 best hits.
            evalue = evalue,
            word_size = 11,
            reward = 2,
            penalty = -3,
            gapopen = 5,
            gapextend = 2,
            outfmt = outfmt_string, # The format we want to use.
        )
        stdout, stderr = blastn_cline()

        if stderr: # If there is no error, then stderr is empty and no error message is printed.
            print("Error:", stderr)
            log_entry = f"Error: BLAST failed for {species} - {stderr}"
            failed_blasts.append(log_entry)
            
        elif os.path.getsize(blast_output) == 0: # Checks if blast_output is empty or not (i.e. if its size is 0).            
            analysis_completed(finished_accession_numbers, analysed_accession_numbers) # Accession numbers are added to {species}_analysed_accession_numbers.txt once BLAST is done.
            nothing_new_to_blast.append(f"{species.replace('_', ' ')}")
            print(f"BLAST completed for {species.replace('_', ' ')}. However, without any new results.") # If blast_output is empty, then the loop ends here. ### Old code. Kept as reference.
        
        else:
            blast_df = columnize(blast_output) # Creates a dataframe based on blast_output.
            blast_df['length_minus_gapopen'] = blast_df['length'] - blast_df['gapopen'] # Adds another column to the dataframe (df), which in this case simply is length minus gapopen.
            
            blast_result = ( # Modifications to blast_df that are performed and the result is blast_result.
                blast_df
                .sort_values(by='qseqid') # Sorts the content by the column 'qseqid'.
                .loc[lambda x: x['length'] > min_length] # Removes and row where the 'length' is shorter than the user defined min_length.
                .assign(qseq=lambda x: x['qseq'].str.replace('-', '')) # Simply removes all the hyphens (-), if there are any, from the column 'qseq'.
                .pipe(lambda x: x.loc[x.groupby('qseqid')['length_minus_gapopen'].idxmax()]) # Groups up based on 'qseqid' and spares only the row with the highest length_minus_gapopen from each group.
                .assign(identical_seqs=lambda x: x.groupby('qseq')['qseqid'].transform('nunique')) # Adds a column called 'identical_seqs'. Groups up based on 'qseq', and each row within the group with a unique 'qseqid' is counted once and added to the total counter in identical_seqs.
                .drop_duplicates(subset='qseq', keep='first') # Removes rows with duplicate entries of the same 'qseq'.
                .sort_values(by='qseq') # Sorts by 'qseq'. Maybe not necessary?
                .reset_index(drop=True) # Removes added grouping structures and restores the default index back to 0.
            )
            blast_output_basename, extension = os.path.splitext(blast_output) # This removes the .txt extension from blast_output in the new variable called blast_output_basename.
            write_fasta(blast_result, blast_filename=blast_output_basename + "_check.fasta")
            analysis_completed(finished_accession_numbers, analysed_accession_numbers) # Accession numbers are added to {species}_analysed_accession_numbers.txt once BLAST is done (just like in the elif-statement).
            print(f"BLAST completed for {species.replace('_', ' ')}.") 
    else:
        print(f"No blasting here. {file_to_blast} not found.") # Occurs if file_to_blast does not exist (which is the case if no new sequence(s) was found). Might replace this.

    if reuse_sequences: # After BLAST is finished and if sequences are being reused, the temporary files created are then deleted in this step.
        os.remove(file_to_blast)
        os.remove(finished_accession_numbers)

    log_progress = f"Progress {idx} out of {len(new_accession_species)} species. Most recent species: {species.replace('_', ' ')}."
    update_last_line(log_file, log_progress)

loop_duration = round(time.time() - loop_timer, 2)

create_directory("BLAST_results")

loop_finished(failed_blasts, log_file, type_of_download)
if nothing_new_to_blast:
    with open(log_file, "a") as file: # Similar to the loop_finished function, but with specific modification for this one instance.
        file.write(f'However, these species had no BLAST results from their sequences:\n')
        file.write('\n'.join(nothing_new_to_blast) + '\n') # Writes which species had no BLAST results into the log file...
    with open(f"BLAST_results/{run_name}_species_without_BLAST_results.txt", "w") as file:
        file.write('\n'.join(nothing_new_to_blast) + '\n') # ... and into a separate file for a more straight forward review.

print(f"{type_of_download} took {loop_duration} seconds to finish.")


###
###
###


all_species_new_results = [] # Creates an empty list to be filled with new sequences from run_name (if there are any).

for species in new_accession_species: # Loops through the species that has been BLASTed.
    new_results_check = f"Species/{species}/{species}_{run_name}_fasta_blast_results_check.fasta" # Creates a variable based on a file name we wish to know the existance of.
    
    if os.path.exists(new_results_check): # Checks if the file exists. If yes, then the subsequent code is executed.
        with open(new_results_check, "r") as file: # The file is read...
            new_blast_results = file.readlines() # ... and its content is stored in new_blast_results.
            all_species_new_results.extend(new_blast_results) # Extends the list called all_species_new_results with the content from the read file.
    else:
        print(f"Nothing new from {run_name} for {species}.") # Prints a message if there is nothing new from run_name. Optional print statement.

if all_species_new_results: # This 'if' statement makes the subsequent code to not be executed if the list is empty.
    result_file = f"BLAST_results/{run_name}_to_curate.fasta"
    with open(result_file, "w") as file: # Creates and writes content into a file with all the new sequences retrieved from run_name.
        for line in all_species_new_results:
            file.write(line)
    append_and_print_message(log_file,
        f"\n\nNew sequence(s) found! {result_file} is ready to be curated.\n")


    species_in_final_file = set()  # Using a set to ensure uniqueness.

    with open(result_file, "r") as file:
        final_file_content = file.readlines()

        for line in final_file_content:
            if line.startswith('>'):  # Check if it's a header line.
                header = line.strip()
                species_name = header.split('|')[2].split(';')[-1].replace('_', ' ')

                species_in_final_file.add(species_name)

    sorted_final_file = sorted(species_in_final_file)

    # Write unique species names to a file.
    with open(f"BLAST_results/{run_name}_species_with_results.txt", "w") as output_file:
        for species_name in sorted_final_file:
            output_file.write(species_name + "\n")


###
###
###


end_timestamp = datetime.today().strftime("%H:%M:%S") # Creates a variable called end_timestamp with the format HH:MM:SS.
program_duration = round(time.time() - program_timer, 2)
append_and_print_message(log_file,
    f"\n\n\nIt took {program_duration} seconds from start to finish.\n"
    f"The program finished at {end_timestamp}.\n")

print("\nRecommendation:\n")
print("To help curate the database:\n")
print(f"Next code line: python Echopipe_database_curation.py BLAST_results/{run_name}_to_curate.fasta.fasta\n")


database_curation = f"Database_curation/{run_name}/"
create_directory(database_curation) # Creates a directory for files related to database curation if it does not already exist.
error_collection_date = f"{database_curation}{run_name}_species_not_found.txt"
duplicate_collection_date = f"{database_curation}{run_name}_duplicate_species_entries.txt"

# Makes a copy of which species were not found in NCBI's taxonomic database and which ones it recognized as duplicate entries. These copies are put in Database_curation/{run_name}/
with open(error_collection, "r") as outfile:
    with open(error_collection_date, "w") as infile:
        for line in outfile:
            infile.writelines(line)
with open(duplicate_collection, "r") as outfile:
    with open(duplicate_collection_date, "w") as infile:
        for line in outfile:
            infile.writelines(line)
            