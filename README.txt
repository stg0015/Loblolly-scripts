Starting with Merged and Sorted Bam files, Run the following codes in the following order to complete the Loblolly analyses:

1. Run the following Rscripts on merged and sorted Bam files: methylCpG2.R, methylCHGall.R, methylCHH.R. This generates the meth and mat files needed for the exploratory and clock analyses. 

2. Run Loblollygenome.sh. This obtains the P. taeda genome for further codes. 

3. Run LoblollyCpG.sh and then Loblollybedscript.sh. These 2 scripts run cpgplot on the P. taeda genome and then give you the island, shore, and shelf bed files needed to run Genomic ranges in R. Note that you need the GCA_000404065.3_Ptaeda2.0assembly_report.txt file to use the Loblolly CPG code. 

4. Run the following R scripts to replicate the exploratoty and clock analyses: Loblolly Exploratory analysis.R, clock cpg manusript.R, clock CHG manuscript.R, and clock CpG and CHG manuscript.R 

5. Run Loblollyblastmanuscript.sh to search for and blast the top correlated sites from the exploratory analyses and the clock sites in the P. taeda genome. Note that manual cutting of the scaffolds is needed to obtain the 400 bp windows for blast to use when searching for related genomic sequences. 
