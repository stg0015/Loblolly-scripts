#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --mail-user=stg86742@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=LoblollyCpG
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=30
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/stg86742/loblolly/annotation/loblollyCpG.out
#SBATCH --error=/scratch/stg86742/loblolly/annotation/loblollyCpG.err


cd /scratch/stg86742/loblolly/annotation2

module load EMBOSS/6.6.0-GCC-8.3.0-Java-11

# cpg plot is on sapelo2, it is in EMBOSS, so I assume I need to load EMBOSS first, which is done by using the above code

#location
 #/apps/eb/EMBOSS/6.6.0-GCC-8.3.0-Java-11/bin


#how to run:
#genome file: GCA_000404065.3_Ptaeda2.0_genomic.fna


#need sliding window of 100 bp (minimum size of CGI = 200 bp, Minimum average observed to expected ratio = 0.6, Minimum average percentage #of G plus C = 50.0)


cpgplot GCA_000404065.3_Ptaeda2.0_genomic.fna -window 100 -minlen 200 -minoe 0.6 -minpc 50 -graph none GCA_000404065.3_Ptaeda2.0_genomic.fna.cpgplot -outfeat loblollycpg.gff 


