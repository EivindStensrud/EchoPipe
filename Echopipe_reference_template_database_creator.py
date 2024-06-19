#!/usr/bin/env python3

# ---------------------------------------------------------------------------
# This script is for generating a template reference database.
#
# The user may either choose to download sequences based on a species list of their own, or they may provide a fasta file with sequences.
# ---------------------------------------------------------------------------


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

###
###
###

parser = argparse.ArgumentParser(
    prog="EchoPipe - Reference template database creator",
    description="Use this script to create a reference template for downstream analysis in the EchoPipe workflow.",
    epilog="Version 1.0")

group_1 = parser.add_argument_group("Arguments used to generate an uncurated reference template.\nThe required arguments are input_file, -f (--forward), -r (--reverse), -e (--email) and -a (--api_key)")
group_1.add_argument('input_file', type=str, nargs='?',
    help="A txt file or CSV with a list of species names or a fasta file that is to be converted into the reference template database (-p is then required).")
group_1.add_argument('-f', '--forward',
    help="The forward primer used to find region of interest, (5'-3' orientation).")
group_1.add_argument('-r', '--reverse',
    help="The reverse primer used to find region of interest, (5'-3' orientation).")
group_1.add_argument('-e', '--email', type=str, default=None,
    help="Your email if NCBI needs to contact you.")
group_1.add_argument('-a', '--api_key', type=str, default=None,
    help="The user's NCBI API key, allows for faster downloads.")
group_1.add_argument('-q', '--query',
    help='The search result will include the user input search term(s). Example, limit the search to 12s region: -q "AND 12s" followed by the search term. To exclude a term write  "NOT 12s".')
group_1.add_argument('-t', '--threshold', type=int, default=150,
    help="The minimum length of a sequence, including the primer regions. Any sequence shorter than this is discarded. Default cutoff is set to 150 bases.")
group_1.add_argument('-l', '--length', type=int, default=22000,
    help="The longest allowed sequence length for template creation. WARNING: The longer the sequence the more computational power is required to align the sequences.")
group_1.add_argument('-m', '--max', type=int, default=1,
    help="The number of sequences that are downloaded per species. Incresing the number may increase the coverage, while increasing the computational power. The total number of downloaded sequences is recommended to not exceed 500. Default = 1.")
group_1.add_argument('-p', '--provided_sequences', action='store_true', default="",
    help="Use if a fasta file is provided to be used as a reference template.")

group_2 = parser.add_argument_group("Argument used to finish the curated reference template database.")
group_2.add_argument('-C', '--Complete', action="store_true",
    help="Completes the reference template database.")

args = parser.parse_args()

working_directory = "Reference_template_creation/"
to_be_curated = "aligned_sequences_to_curate.fasta"

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
                    f"Minimum sequence length: {args.threshold}\n")
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
                            except Exception as e:
                                print(f"An error occurred:\n{e}")

                if list_of_not_found:
                    with open("no_search_results.txt", "w") as file:
                        for line in list_of_not_found:
                            file.writelines(line)

            ###
            ###
            ###

            # Function to find the position of the first base in the forward primer.
            def find_5_end_fwd_position(aligned_seq):
                position_5_end_fwd = 1
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

            batch_size = 40

            records = list(SeqIO.parse(mafft_input_file, "fasta"))
            num_batches = (len(records) + batch_size - 1) // batch_size

            for i in range(num_batches):
                start_idx = i * batch_size
                end_idx = min((i + 1) * batch_size, len(records))
                output_file = f"{output_prefix}_{i + 1}.fasta"
                with open(output_file, "w") as file:
                    SeqIO.write(records[start_idx:end_idx], file, "fasta")
                    file.write(f">Forward_primer\n{forward_primer}\n")
                    file.write(f">Reverse_primer\n{reverse_primer}\n")


            input_prefix = "temp_split_file_"
            maffted_file = "aligned_sequences.fasta"

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
                                if seq_record.id == "Forward_primer" or seq_record.id == "Reverse_primer":
                                    trimmed_and_filtered_sequences.append(f">{str(seq_record.description)}\n")
                                    trimmed_and_filtered_sequences.append(str(seq_record.seq.replace("-", "")) + "\n")
                                elif len(seq_record.seq[start_pos_forward:end_pos_reverse].replace("-", "")) > length_threshold:
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

            with open(to_be_curated, "w") as file:
                mafft_cline = MafftCommandline(
                    input=maffted_file,
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

    print("reference_template_database.fasta has been created and ready for use with ECHoPipe_db_creation.py")

