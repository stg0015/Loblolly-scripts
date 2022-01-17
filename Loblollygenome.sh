#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --mail-user=stg86742@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=Loblollygenome
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/stg86742/loblolly/annotation/loblollygenome.out
#SBATCH --error=/scratch/stg86742/loblolly/annotation/loblollygenome.err


cd /scratch/stg86742/loblolly/annotation

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/065/GCA_000404065.3_Ptaeda2.0/GCA_000404065.3_Ptaeda2.0_genomic.fna.gz


gunzip GCA_000404065.3_Ptaeda2.0_genomic.fna.gz

