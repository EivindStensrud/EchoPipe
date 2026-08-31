#!/usr/bin/env Rscript

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2); packageVersion("ggplot2")
   }
if (!require("xlsx")) {
   install.packages("xlsx", dependencies = TRUE)
   library(xlsx)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("R.utils")) {
   install.packages("R.utils", dependencies = TRUE)
   library(R.utils)
   }
if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape2)
   }
if (!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
  }
if (!require("psych")) {
   install.packages("psych", dependencies = TRUE)
   library(psych)
   }
 if (!require("optparse")) {
   install.packages("optparse", dependencies = TRUE)
   library(optparse)
   }

option_list = list(
   make_option(c("-i", "--Input_path"), type="character", default=NULL, 
              help="path to fasta.otu.txt", metavar="character"),
   make_option(c("-o", "--Output_path"), type="character", default=NULL, 
              help="path to output", metavar="character")
  
)
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$Input_path)){
  print_help(opt_parser)
  stop("Input argument(s) missing.n", call=FALSE)
}

#Reformats the annotation to work with LCA annotation.

RootPath <- file.path(opt$Input_path) # Change to the path, where the result from DADA2 pipeline or similar pipelines is stored.
FigsPath <- file.path(opt$Output_path) # Path where output will be stored.


print("RootPath")
print(RootPath)

blast_result_file = (list.files(path = RootPath, pattern ="output_blast_results$",full.names=TRUE))

blast_results_db = read.table(file=blast_result_file) # Results from BLAST analysis above.
colnames(blast_results_db) <- c( "Query ID", "Subject", "Identity percentage", 
"Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end","Coverage",
"S.start", "S.end", "Evalue", "Bitscore" ) # EchoPipe build-up, header goes into Subject


tax_db = data.frame(stringr::str_replace_all(blast_results_db$Subject, ";", " / ")) # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis
colnames(tax_db) = c("Subject")
taxonomy = sapply(strsplit(tax_db$Subject, split="\\|"), function(x) {
   # The taxonomy part is the third element in the split list
   taxonomy = x[3]
  
   # Remove potential trailing characters (like numbers and the last pipe) after the taxonomy
   taxonomy_clean = sub("\\|.*$", "", taxonomy)
  
   # Return the clean taxonomy
   return(taxonomy_clean)
})

blast_results_db$"Taxonomy" = taxonomy # creates a column with taxonomy
blast_results_db$Subject = (stringr::str_replace_all(blast_results_db$Subject, ";", "/")) # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis
blast_results_db$Subject = (stringr::str_replace_all(blast_results_db$Subject, "[[|]]", "/")) # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis


tax_acc = data.frame(unique(blast_results_db$Subject))
colnames(tax_acc)= c("Subject")
tax_accno = str_extract(tax_acc$Subject, "(?<=gb\\/)[^/]+")


tax_accno = data.frame(cbind(tax_acc, tax_accno))
colnames(tax_accno) = c("Subject", "Subject accession") # Stores accession number and subject ID in same column.

blast_results_db = left_join(blast_results_db, tax_accno, by = "Subject") # Joins blast_results by accession and subject ID.
blast_results_db$"Subject Taxonomy ID" = blast_results_db$"Subject" # Data formating
blast_results_db$"Source" = "EchoPipe"

blast_results_db_final = data.frame(cbind(blast_results_db$"Query ID",
                                          blast_results_db$"Subject",
                                          blast_results_db$"Subject accession",
                                          blast_results_db$"Subject Taxonomy ID",
                                          blast_results_db$"Identity percentage",
                                          blast_results_db$"Coverage",
                                          blast_results_db$"Evalue",
                                          blast_results_db$"Bitscore",
                                          blast_results_db$"Source",
                                          blast_results_db$"Taxonomy"
                                          ))
colnames(blast_results_db_final) = c("#Query ID",
                                     "#Subject",
                                     "#Subject accession",
                                     "#Subject Taxonomy ID",
                                     "#Identity percentage",
                                     "#Coverage",
                                     "#Evalue",
                                     "#Bitscore",
                                     "#Source",
                                     "#Taxonomy"
                                          )



write.table(blast_results_db_final,file.path(FigsPath, file="blast_result_file.lca.tabular"), row.names = FALSE, quote = FALSE, sep = "\t") # Stores the BLAST results in a table accepted by the LCA analysis.


