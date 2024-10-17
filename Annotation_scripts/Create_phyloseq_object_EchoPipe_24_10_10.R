#!/usr/bin/env Rscript

if (!require("dada2")) {
   BiocManager::install("dada2")
   library(dada2); packageVersion("dada2")
   }

if (!require("phyloseq")) {
   BiocManager::install("phyloseq", dependencies = TRUE)
   library(phyloseq); packageVersion("phyloseq")
   }

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
if (!require("limma")) {
   BiocManager::install("limma")
   library(limma)
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
   make_option(c("-i", "--Input_path"), type="character", default="", 
      help="path to fasta.otu.txt", metavar="character"),
   make_option(c("-o", "--Output_path"), type="character", default="", 
      help="path to output dir", metavar="character"),
   make_option(c("-b", "--Blast_id_cutoff"), type="double", default=99, 
      help="Blast species annotation identity cutoff", metavar="double"),
   make_option(c("-c", "--Blast_cov_cutoff"), type="double", default=99, 
      help="Blast species annotation identity cutoff", metavar="double")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$Input_path)){
  print_help(opt_parser)
  stop("Input argument(s) missing.n", call=FALSE)
}

RootPath <- file.path(opt$Input_path) # Change to the path, where the result from DADA2 pipeline or similar pipelines is stored.
FigsPath <- file.path(opt$Output_path) # Path where output will be stored. Not in use atm


Blast_id_cutoff=(opt$Blast_id_cutoff)
Blast_cov_cutoff=(opt$Blast_cov_cutoff)


seqtab_file = (list.files(path = RootPath, pattern = "seqtab_nochim.rds$" ,full.names=T, recursive = F)) # use seqtab_nochim.rds

seqtab = readRDS(file=seqtab_file) # Input is the ASV table.

blast_file = (list.files(path = RootPath, pattern ="output_blast_results$",full.names=T, recursive = F))


blast_results <- read.table(blast_file, fill=TRUE, sep ="\t") # Read Blast result from above.

colnames(blast_results) <- c( "QueryID",  "SubjectID", "Perc.Ident",
"Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end", "Perc.Cov",
"S.start", "S.end", "E", "Bits" ) #Order Blast results.


#keeps hit with best bit score and better alignment score than Blast_id_cutoff and Blast_cov_cutoff
blast_results_clean <- blast_results %>% filter(as.numeric(Perc.Ident) >= as.numeric(Blast_id_cutoff)) %>% filter(as.numeric(Perc.Cov) >= as.numeric(Blast_cov_cutoff))

#keeps only the best bit-score per query sequence, aka best hit based on bit-score, could potentially be changed to which.min(as.numeric(d$E))
blast_results_clean = do.call( rbind,
        lapply( split(blast_results_clean, blast_results_clean[,c("QueryID") ]), function(d){
                                         d[which.max(as.numeric(d$Bits)), ] } )
        )

#keep taxonomy from blast annotation
tax_blast=data.frame(blast_results_clean)

#formats subject id to taxa.
tax_blast[c("Origin", "Accesion_number","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "BLAST_Species", "Counter")]<- str_split_fixed(tax_blast$SubjectID, "[|;]",10 )
names(tax_blast)[names(tax_blast) == "QueryID"] <- "Query" # makes it compatible with LCA downstream.
names(tax_blast)[names(tax_blast) == "Perc.Ident"] <- "Blast_Perc.Ident" # makes it compatible with LCA downstream.
names(tax_blast)[names(tax_blast) == "Perc.Cov"] <- "Blast_Perc.Cov" # makes it compatible with LCA downstream.

tax_blast=data.frame(tax_blast)

lca_file = (list.files(path = RootPath, pattern =".lca.tabular_04_98_out$",full.names=T, recursive = F))

lca_table = read.table(lca_file, sep="\t",fill=TRUE) # Read LCA output file
colnames(lca_table) = c("Query","LCA_rank","LCA_taxon","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Method", "Identity", "Coverage") # Combines both LCA and top blast assignment. 

#insert, combine species name in top hit result.
# Function to process each group of sequences

combine_identical_tophits = function(group) {
  # Combine species names for all rows in the group
  combined_query = paste(unique(group$LCA_taxon), collapse = "_")
  group$LCA_taxon = combined_query
  
  # Retain the first occurrence of other columns
  first_row = group[1, ]
  first_row$LCA_taxon = combined_query
  
  return(first_row)
}

# Process the data frame, merge and keep tophits.
lca_results_clean = do.call(rbind, lapply(split(lca_table, lca_table$Query), combine_identical_tophits))


tax_lca = data.frame(lca_results_clean) #makes as matrix

#Works correctly.
tax_merged <- merge(tax_lca, tax_blast[, c("Query", "BLAST_Species", "Blast_Perc.Ident", "Blast_Perc.Cov" )], by = "Query", all.x = TRUE)


tax = as.matrix(tax_merged)
rownames(tax) <- tax[,1] # makes first col as rownames
tax = tax[,-1] #removes col which became rownames
tax = as.matrix(tax) # makes as matrix





################## parse seqtab ####################
#only keeps blast hits.
# keeps all Querys from input blast file, identical to number of queries from LCA.

seqtab_clean = seqtab[, colnames(seqtab) %in% blast_results$QueryID] 

seqtab_clean_names = sub('.*\\/', '', rownames(seqtab_clean))
seqtab_clean_names = sub('\\..*', '', (seqtab_clean_names)) #removes .fq.gz from file name


rownames(seqtab_clean) = seqtab_clean_names


################# sequence analysis evaluation ######
#track <- read.table(file.path(RootPath, "track_sequence_stats.csv"),
#                                   sep = ",", header = TRUE)

#summary(track)

#track$sample_id = str_extract(track$X, "[^_]+")

#seqtab_samplesum = rowSums(seqtab_clean)


################# Sequence data analysis and visualization ###############

#using only seqtab, not seqtab_clean keeps all sequences, good to have to check what organisms are not annotated.

am_otu_n <- otu_table(t(seqtab_clean), taxa_are_rows = TRUE)
#dim(am_otu_n)
am_tax <- tax_table(as.matrix(tax))
am_physeq_clean_n <- phyloseq(am_otu_n, am_tax) # Combines ASV/OTU and Taxonomic table.

phyloseq_name="clean_phyloseq_mifish_lca.rds"
saveRDS(am_physeq_clean_n, file.path(RootPath, file=phyloseq_name)) ## Phyloseq object to use for further analysis on local computer.

print("Use the xx.phyloseq_mifish_lca.rds for downstream analysis")


#using only seqtab, not seqtab_clean keeps all sequences, good to have to check what organisms are not annotated.

#' @title Summarize taxon composition
#'
#' @description This function takes a phyloseq object and returns phyloseq data (OTUs are now taxonomic ranks)
#' @param phylo_seq_object A phyloseq object with an OTU table and phylogenetic tree.
#' @param taxonomic_rank The taxonomic rank at which the data should be summarized
#' @details
#' Nice for making taxon summaries
#' @return A data.frame with the taxon name the mean, standard deviation as well as min and max values per samples.
#' @keywords phyloseq
#' @export
#' @examples
#' phyloseq_summarize_taxa()

phyloseq_summarize_taxa <- function(phylo_seq_object, taxonomic_rank = rank_names(phylo_seq_object)[1], errorIfNULL = TRUE) {
  taxa <- as(phyloseq::tax_table(phylo_seq_object, errorIfNULL)[, taxonomic_rank], 'character')
  sum_tax_table <- summarize_taxa(as(phyloseq::otu_table(phylo_seq_object), 'matrix'), taxa)
  phyloseq::phyloseq(phyloseq::otu_table(sum_tax_table, taxa_are_rows = TRUE),
   phyloseq::sample_data(phylo_seq_object, FALSE))
}

#' @title Summarize taxon composition
#'
#' @description This function takes a phyloseq object and returns phyloseq data (OTUs are now taxonomic ranks)
#' @param phylo_seq_object A phyloseq object with an OTU table and phylogenetic tree.
#' @param taxonomic_rank The taxonomic rank at which the data should be summarized
#' @details
#' Nice for making taxon summaries
#' @return A data.frame with the taxon name the mean, standard deviation as well as min and max values per samples.
#' @keywords phyloseq
#' @export
#' @examples
#' summarize_taxa()

summarize_taxa <- function(counts, taxonomy) {
  if(is.matrix(taxonomy)) {
    #message('multiple taxonomies')
    plyr::alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    #message('multiple counts')
    require('plyr')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    #message('summarize')
    tapply(counts, taxonomy, sum)
  }
}
################# Sequence data analysis and visualization ###############

lca_taxon_table_n <- otu_table(phyloseq_summarize_taxa(am_physeq_clean_n, 'LCA_taxon'))
#lca_taxon_table_n[lca_taxon_table_n < 11] = 0
#lca_taxon_table_n[lca_taxon_table_n > 10] = 1

write.csv(t(lca_taxon_table_n), file.path(FigsPath, "lca_taxon_table_n.csv"))

species_table_n <- otu_table(phyloseq_summarize_taxa(am_physeq_clean_n, 'BLAST_Species'))
#species_table_n[species_table_n < 11] = 0
#species_table_n[species_table_n > 10] = 1

write.csv(t(species_table_n), file.path(FigsPath, "species_table_n.csv"))

