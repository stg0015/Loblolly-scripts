#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-user=stg86742@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=LoblollyBlast
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=30
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/stg86742/loblolly/annotation3/loblollyblast.out
#SBATCH --error=/scratch/stg86742/loblolly/annotation3/loblollyblast.err

module load BLAST/2.2.26-Linux_x86_64
cd /scratch/stg86742/loblolly/annotation3/GCA_000404065.3_Ptaeda2.0

#prepared genome to search as Ptaeda2.01.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < GCA_000404065.3_Ptaeda2.0_genomic.fna > Ptaeda2.01.fa

tail -n +2 Ptaeda2.01.fa > file2.fa


#use this to pull out the scaffolds that contain the methylated cytosines
awk -F '"' '{print $2}' your_input_file > listofsites


#now obtain the scaffold sequences from the above  list
while read i;
do
echo "searching for "$i""
LC_ALL=C grep -A 1 "$i" file2.fa > "$i"searchsite
done<listofsites


#manually pull out the 400 bp regions that contain the cytosine of interest from the exploratory and clock analyses using: cat "$i"searchsite | cut -c (sites here) > "$i"400bpwindow

#need to make a list of all the sites: ex: ls *400bpwindow > 400bpsites

 
#now for the final while loop that blast searches all the 400bpwindows



while read i;
do
blastall -p blastn -d /db/ncbiblast/nrte/latest/nt -a 4 -i "$i" -o "$i"blast.out

#the below are different lists that were used to find all the CpG and CHG sites from the exploratory and clock analyses
#done<remainingcpglist
#done<finalCHGsites
#done<finallistboth
#done<final400cpglist
#done<finalspglist
done<400bpwindows

