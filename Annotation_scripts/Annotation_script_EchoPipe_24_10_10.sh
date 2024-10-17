#!/bin/bash
# Job name: EchoPipe annotation, 241010
#SBATCH --job-name=EchoPipe
#SBATCH --account=nn9744k
#SBATCH --time=2:30:00
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

### 
### seq, annotation MiFish, LCA + blast
# sequencing run
# 2024_10_10 sequencing run
# Fish Ethiopia

Help()
{
   #Display Help
   echo "Unique Molecular Identifier based error correction for eDNA metabarcoding."
   echo " "
   echo "The script expects files to be in a fastq.gz format."
   echo "The script expects read1 and read2 files being seperated as R1 and R2 files respectively."
   echo " "
   echo "Arguments: [-h|-i|-o|-s|-r|-n]"
   echo " "
   echo "Options:"
   echo "-h    Display this help."
   echo "-i    Path to directory where rawdata file is stored (fastq.gz), char."
   echo "-o    Path to directory to store output data, char"
   echo "-s    Path to script"
   echo "-r    Path to the folder of reference database"
   echo "-n    Name of reference database in the path"
}


##Setup job environment
module purge   # clear any inherited modules
#set -o errexit # exit on errors


while getopts "hi:o:s:r:n:b:c:z:l:v:" option; do
   case $option in
      h) # display Help
         Help
         exit 0
         ;;
      i)
         Input_path="$OPTARG"
         ;;
      o)
         Output_path="$OPTARG"
         ;;
      s)
         Script_path="$OPTARG"
         ;;
      r)
         Reference_database_path="$OPTARG"
         ;;
      n) 
         Reference_database_name="$OPTARG"
         ;;
      b)
         Blast_id_cutoff="$OPTARG"
         ;;
      c)
         Blast_cov_cutoff="$OPTARG"
         ;;
      z)
         LCA_bitscore="$OPTARG"
         ;;
      l)
         LCA_id_cutoff="$OPTARG"
         ;;
      v)
         LCA_cov_cutoff="$OPTARG"
         ;;


   esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # set curret directory as working directory if non is selected.
#SCANDIFISH_DIR=/cluster/projects/nn9745k/03_databases/fish/ScandiFish_12s_v1.4/ # set reference database directory of ScandiFish if non is selected. # on saga

if [ -z "$Input_path" ]; then Input_path=; fi #where the rep-seqs.fna output file is
if [ -z "$Output_path" ]; then Output_path=; fi #Where to store the output, standard the same as $Input_path
if [ -z "$Script_path" ]; then Script_path=; fi #if the script directory is not selected, script path is as standard for saga
if [ -z "$Reference_database_path" ]; then Reference_database_path=; fi #if the reference database is not selected, ScandiFish v1.4 will be used
if [ -z "$Reference_database_name" ]; then Reference_database_name=""; fi #if the reference database is not selected, ScandiFish v1.4 will be used
if [ -z "$Blast_id_cutoff" ]; then Blast_id_cutoff=99; fi #Blast species annotation identity cutoff
if [ -z "$Blast_cov_cutoff" ]; then Blast_cov_cutoff=99; fi #Blast species annotation coverage cutoff
if [ -z "$LCA_bitscore" ]; then LCA_bitscore=8; fi #LCA Bitscore top percentage threshold
if [ -z "$LCA_id_cutoff" ]; then LCA_id_cutoff=95; fi #LCA Bitscore top percentage threshold
if [ -z "$LCA_cov_cutoff" ]; then LCA_cov_cutoff=95; fi #LCA Bitscore top percentage threshold


#download reference database
cd $Reference_database_path

wget -nc https://raw.githubusercontent.com/EivindStensrud/ScandiFish/main/ScandiFish_12s_v1.4/ScandiFish_12s_v1.4_nf.fasta -P $Reference_database_path #Downloads database
wget -nc https://raw.githubusercontent.com/EivindStensrud/ScandiFish/main/ScandiFish_12s_v1.4/ScandiFish_12s_v1.4_nf.fasta.md5 -P $Reference_database_path # Downloads md5sum

md5sum -c $Reference_database_path/ScandiFish_12s_v1.4_nf.fasta.md5 > $Output_path/checkmd5.md5 # Creates files which can be checked if md5sum is correct. Needs to updated for new database version
# Downloads LCA script from GitHub.

printf "Initiation the annotation using BLAST with $Reference_database_path/$Reference_database_name \n"
echo " "
module purge

module load BLAST+/2.13.0-gompi-2022a

## Assigning taxonomy by BLASTN
#
##makeblastdb -in ScandiFish_12s_v1.4_nf.fasta -out ScandiFish_12s_v1.4_nf_db -dbtype nucl -title "12S ScandiFish_12s_v1.4 Scandinavian fish database" #Input is the altered format of the ScandiFish database. # only needs to be ran once.
##blast

blastn -max_target_seqs 100 -evalue 1 -query $Input_path/rep-seqs.fna -out $Input_path/output_blast_results -db $Reference_database_path/$Reference_database_name -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qcovs sstart send evalue bitscore" -num_threads 2 # BLASTN compares the sequences and keeps up to 100 reference sequences. # Input is DADA2 output from function uniquesToFasta(seqtab.nochim).

#blastn -max_target_seqs 100 -evalue 1 -query /cluster/projects/nn9745k/02_results/47_Nick/241009/Figs/rep-seqs.fna -out /cluster/projects/nn9745k/02_results/47_Nick/241009/Figs/output_blast_results -db /cluster/projects/nn9745k/03_databases/fish/Ethiopia_12s_241010/Ethipoia_Scandi_241015/Ethipoia_Scandi_241015_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qcovs sstart send evalue bitscore" -num_threads 2 # BLASTN compares the sequences and keeps up to 100 reference sequences. # Input is DADA2 output from function uniquesToFasta(seqtab.nochim).



echo " "
printf "Initiation the Rscript, Blast_annotation.R with \n"
echo " "
module purge

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

Rscript $Script_path/Blast_annotation_fish.R -i $Input_path -o $Output_path # starts r script to create sequence table in dada2 and phyloseq compatible format.

printf "Finished with annotation using the Rscript Blast_annotation.R \n"


wget -nc https://raw.githubusercontent.com/naturalis/galaxy-tool-lca/master/lca.py -P $Script_path #downloads lca.py if it do not exist

## Run LCA analysis in bash # run example 3 from link https://github.com/naturalis/galaxy-tool-lca
## Make sure to empty table before re-running the code, otherwise will append data, not overwrite.

rm $Output_path/blast_result_file.lca.tabular_04_98_out # Removes table if existing LCA analysis exist.

## LCA analysis
## Best hit goes to species, taxonomic sorting is conducted if top hit <99%


python $Script_path/lca.py -i $Output_path/blast_result_file.lca.tabular -o $Output_path/blast_result_file.lca.tabular_04_98_out -b $LCA_bitscore -id $LCA_id_cutoff -cov $LCA_cov_cutoff -t best_hits_range -tid $Blast_id_cutoff -tcov $Blast_cov_cutoff -flh unknown

module purge

module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2


echo " "
echo "Finished with LCA annotation"
echo " "


#RAWDIR=$Input_path
#RESDIR=$Input_path/Figs

Rscript $Script_path/Create_phyloseq_object.R -i $Input_path -o $Output_path -b $Blast_id_cutoff -c $Blast_cov_cutoff
Rscript $Script_path/Create_phyloseq_object.R -i $RAWDIR -o $RAWDIR -b $Blast_id_cutoff -c $Blast_cov_cutoff

echo " "
echo " "
echo "Finished to make Phyloseq object"
echo " "
echo " "
printf "The script was finished %(%Y-%m-%d %H:%M:%S)T\n"
