#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-user=stg86742@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=LoblollyBed
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=30
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/stg86742/loblolly/annotation2/loblollyCpGBED.out
#SBATCH --error=/scratch/stg86742/loblolly/annotation2/loblollyCpGBED.err


cd /scratch/stg86742/loblolly/annotation2

module load SAMtools/1.10-iccifort-2019.5.281
module load BEDTools/2.29.2-GCC-8.3.0

# cpgplot generates a gff file with all the locations of the cpg islands
# we want to make this a bed file by cutting the "chromosome", start, and end position
grep -v '#' loblollycpg.gff | cut -f1,4,5 > Loblolly_cpg_islands.bed




#use this to obtain info for lengths of each scaffold
#(Note: I had to download the GCA_* file below from the ncbi genome page where the Loblolly genome is, then upload it
grep -v '#' GCA_000404065.3_Ptaeda2.0_assembly_report.txt | cut -f5,9 > Loblolly.ref.fa


#next step is to use bedtools slop to extend the coordinates of the CpG islands by 2000 in each direction
#bedtools requires a file which contains the lengths of each chromosome so that it doesn't create coordinates outside the possible range


#the first 2 columns of the ref.fa file contain the scaffold name and length
#send the lengths only to a new file
cut -f-2 Loblolly.ref.fa > Loblolly.ref.fa.g

#now generate the shores
slopBed -i Loblolly_cpg_islands.bed -g Loblolly.ref.fa.g -b 2000 | mergeBed \
| subtractBed -a - -b Loblolly_cpg_islands.bed > Loblolly_cpg_shores.bed

#to get shelves I need to get the intermediate file of the islands and shores together
slopBed -i Loblolly_cpg_islands.bed -g Loblolly.ref.fa.g -b 2000 | mergeBed > Loblolly_cpg_island_shore.bed
#then I need to extend that out and remove the island-shore sections
slopBed -i Loblolly_cpg_island_shore.bed -g Loblolly.ref.fa.g -b 2000 | mergeBed \
| subtractBed -a - -b Loblolly_cpg_island_shore.bed > Lobolly_cpg_shelves.bed


#need to import bed files into R and convert to GRanges object
#need to find overlap between covered CpGs and each of these features
