#!/usr/bin/env python3

### This script will reformat the headers of the database to specific formats.

from Bio import SeqIO
import os
import argparse
from argparse import RawDescriptionHelpFormatter

###
###
###

parser = argparse.ArgumentParser(
	prog="EchoPipe - Reformat",
	description= (
        "Reformats the input database to a pre-defined format:\n"
        "sintax = SINTAX\n"
        "rdp = Ribosomal Database Project (RDP)\n"
        "dadt = DADA2 assignTaxonomy\n"
        "dads = DADA2 assignSpecies\n"
        "idt = IDTAXA\n"
        "qiime = QIIME 2\n"
    ),
    formatter_class=RawDescriptionHelpFormatter,
	epilog="Version 1.0")

parser.add_argument('reference_database', type=str,
	help="The reference database for which the header format will be changed.")
parser.add_argument('format', type=str, choices=['sintax', 'rdp', 'dadt', 'dads', 'idt', 'qiime'],
	help="Write in one of the available formats as seen above.")

args = parser.parse_args()

reference_database = args.reference_database
format = args.format

###
###
###

def write_new_fasta(header):
    file.write(header + "\n")
    file.write(sequence + "\n")

###
###
###

database_name = os.path.splitext(reference_database)[0] # Removes the file extension from the file name.
reformated_database = f"{database_name}_{format}.fasta"

with open(reformated_database, "w") as file:
    for record in SeqIO.parse(reference_database, "fasta"):
        split_header = record.id.split("|")
        accession_number = split_header[1]
        lineage = split_header[2].split(";")
        sequence = str(record.seq)
        
        if format == "sintax":            
            sintax_header = f">{accession_number};tax=d:{lineage[0]},p:{lineage[1]},c:{lineage[2]},o:{lineage[3]},f:{lineage[4]},g:{lineage[5]},s:{lineage[6]}"
            write_new_fasta(sintax_header)

        elif format == "rdp":
            rdp_header = f">{accession_number}\troot;{lineage[0]};{lineage[1]};{lineage[2]};{lineage[3]};{lineage[4]};{lineage[5]};{lineage[6]}"
            write_new_fasta(rdp_header)

        elif format == "dadt":
            dad_header = f">{lineage[0]};{lineage[1]};{lineage[2]};{lineage[3]};{lineage[4]};{lineage[5]}"
            write_new_fasta(dad_header)

        elif format == "dads":
            dads_header = f">{accession_number} {lineage[5]} {lineage[6]}"
            write_new_fasta(dads_header)

        elif format == "idt":
            idt_header = f">Root;{lineage[0]};{lineage[1]};{lineage[2]};{lineage[3]};{lineage[4]};{lineage[5]};{lineage[6]}"
            write_new_fasta(idt_header)

        elif format == "qiime":
            qiif_header = f">{accession_number}"
            write_new_fasta(qiif_header)
                
if format == "qiime": # This one is a bit special as it also creates a text file.
    qiif_txt_file = f"{database_name}_{format}.txt"
    with open(qiif_txt_file, "w") as file:
        for record in SeqIO.parse(reference_database, "fasta"):
            split_header = record.id.split("|")
            accession_number = split_header[1]
            lineage = split_header[2].split(";")
            species = "_".join(lineage[6].split("_")[1:])
            qiif_txt_header = f"{accession_number}\tk__{lineage[0]}; p__{lineage[1]}; c__{lineage[2]}; o__{lineage[3]}; f__{lineage[4]}; g__{lineage[5]}; s__{species}"
            file.write(qiif_txt_header + "\n")

print("\nReformatting completed.\n")