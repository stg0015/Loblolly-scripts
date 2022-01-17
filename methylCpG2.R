#part 1 of RRBS analysis

library(methylKit)


#reading in files
RRBSfile.list=list("/scratch/stg86742/loblolly/ALL_1_merged.bam", "/scratch/stg86742/loblolly/ALL_2_merged.bam", "/scratch/stg86742/loblolly/ALL_3_merged.bam", "/scratch/stg86742/loblolly/ALL_4_merged.bam","/scratch/stg86742/loblolly/ALL_5_merged.bam", "/scratch/stg86742/loblolly/ALL_6_merged.bam",
                   "/scratch/stg86742/loblolly/ALL_7_merged.bam", "/scratch/stg86742/loblolly/ALL_8_merged.bam", "/scratch/stg86742/loblolly/ALL_9_merged.bam", "/scratch/stg86742/loblolly/ALL_10_merged.bam", "/scratch/stg86742/loblolly/ALL_11_merged.bam", "/scratch/stg86742/loblolly/ALL_12_merged.bam",
                   "/scratch/stg86742/loblolly/ALL_13_merged.bam", "/scratch/stg86742/loblolly/ALL_14_merged.bam", "/scratch/stg86742/loblolly/ALL_15_merged.bam", "/scratch/stg86742/loblolly/ALL_16_merged.bam", "/scratch/stg86742/loblolly/ALL_17_merged.bam", "/scratch/stg86742/loblolly/ALL_18_merged.bam",
                   "/scratch/stg86742/loblolly/ALL_19_merged.bam", "/scratch/stg86742/loblolly/ALL_20_merged.bam", "/scratch/stg86742/loblolly/ALL_21_merged.bam", "/scratch/stg86742/loblolly/ALL_23_merged.bam", "/scratch/stg86742/loblolly/ALL_26_merged.bam", "/scratch/stg86742/loblolly/ALL_28_merged.bam")


#initially trying for all types of cytosine methylation (CpG, CHG, CHH)
#running this code for CpG
RRBSobjs=processBismarkAln(location=RRBSfile.list
                           ,sample.id=list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "26", "28"),
                           assembly="ASM223467v1",
                           save.folder="callfoldercpg",save.context="CpG",read.context="CpG",
                           nolap=FALSE,mincov=1,minqual=20,phred64=FALSE,
                           treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

#filtering based on coverage
#filtering for 10, 8, and 5
#lo.count=8 indicates that we filter any sites/bases with less than 8x coverage
#hi.perc-99.9 tells it to discard bases with more than 99.9% of coverage in samples (completely methylated)
filtered.myobj10=filterByCoverage(RRBSobjs,lo.count=10,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

filtered.myobj8=filterByCoverage(RRBSobjs,lo.count=8,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)


filtered.myobj5=filterByCoverage(RRBSobjs,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)


#merging all files to perform comparative analysis of bases methylated across all samples
# Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide.
#This provides better coverage, but only advised when looking at CpG methylation
#(for CpH methylation this will cause wrong results in subsequent analyses). I may try this jjst to see what happens

#running for both destrand = True and False; 
#also running for min.per.group = 22, and 19 (this is roughly 90 and 80 percent covered across all samples)

#10x
meth1=unite(filtered.myobj10, destrand=TRUE) #gives 100 percent covered across all samples

meth2=unite(filtered.myobj10, destrand=TRUE, min.per.group=22L) #gives 90 percent covered

meth3=unite(filtered.myobj10, destrand=TRUE, min.per.group=19L) #gives 80 percent covered

#8x
meth4=unite(filtered.myobj8, destrand=TRUE) #gives 100 percent covered across all samples

meth5=unite(filtered.myobj8, destrand=TRUE, min.per.group=22L) #gives 90 percent covered

meth6=unite(filtered.myobj8, destrand=TRUE, min.per.group=19L) #gives 80 percent covered

#5x
meth7=unite(filtered.myobj5, destrand=TRUE) #gives 100 percent covered across all samples

meth8=unite(filtered.myobj5, destrand=TRUE, min.per.group=22L) #gives 90 percent covered

meth9=unite(filtered.myobj5, destrand=TRUE, min.per.group=19L) #gives 80 percent covered



#now for destrand equal false
#10x
meth11=unite(filtered.myobj10, destrand=FALSE)

meth12=unite(filtered.myobj10, destrand=FALSE, min.per.group=22L)

meth13=unite(filtered.myobj10, destrand=FALSE, min.per.group=19L)

#8x
meth14=unite(filtered.myobj8, destrand=FALSE)

meth15=unite(filtered.myobj8, destrand=FALSE, min.per.group=22L)

meth16=unite(filtered.myobj8, destrand=FALSE, min.per.group=19L)

#5x
meth17=unite(filtered.myobj5, destrand=FALSE)

meth18=unite(filtered.myobj5, destrand=FALSE, min.per.group=22L)

meth19=unite(filtered.myobj5, destrand=FALSE, min.per.group=19L)






#correlations of CpG sites among samplesvv(initially run without this, but maybe later try to get plots)
#getCorrelation(meth,plot=TRUE)


#clustering the samples based on how similar their methylation patterns are:
#clusterSamples(meth, dist="correlation", method="ward.D", plot=TRUE)

#PCA analysis
#PCASamples(meth, screeplot=TRUE)


#PCASamples(meth)


##################################################################################################################

##making matrix with methylation percentages

#need to have 18 of each (meth 1-9 = destrand = True, meth 11-19 destrand = F)

#meth1
mat1=percMethylation(meth1)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat1, file= "matCpG10destrue100.csv")
write.csv(meth1, file="methCpG10destrue100.csv")

#meth2
mat2=percMethylation(meth2)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat2, file= "matCpG10destrue90.csv")
write.csv(meth2, file="methCpG10destrue90.csv")


#meth3
mat3=percMethylation(meth3)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat3, file= "matCpG10destrue80.csv")
write.csv(meth3, file="methCpG10destrue80.csv")


#meth4
mat4=percMethylation(meth4)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat4, file= "matCpG8destrue100.csv")
write.csv(meth4, file="methCpG8destrue100.csv")


#meth5
mat5=percMethylation(meth5)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat5, file= "matCpG8destrue90.csv")
write.csv(meth5, file="methCpG8destrue90.csv")


#meth6
mat6=percMethylation(meth6)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat6, file= "matCpG8destrue80.csv")
write.csv(meth6, file="methCpG8destrue80.csv")


#meth7
mat7=percMethylation(meth7)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat7, file= "matCpG5destrue100.csv")
write.csv(meth7, file="methCpG5destrue100.csv")


#meth8
mat8=percMethylation(meth8)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat8, file= "matCpG5destrue90.csv")
write.csv(meth8, file="methCpG5destrue90.csv")



#meth9
mat9=percMethylation(meth9)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat9, file= "matCpG5destrue80.csv")
write.csv(meth9, file="methCpG5destrue80.csv")


#meth11
mat11=percMethylation(meth11)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat11, file= "matCpG10desfalse100.csv")
write.csv(meth11, file="methCpG10desfalse100.csv")


#meth12
mat12=percMethylation(meth12)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat12, file= "matCpG10desfalse90.csv")
write.csv(meth12, file="methCpG10desfalse90.csv")


#meth13
mat13=percMethylation(meth13)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat13, file= "matCpG10desfalse80.csv")
write.csv(meth13, file="methCpG10desfalse80.csv")


#meth14
mat14=percMethylation(meth14)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat14, file= "matCpG8desfalse100.csv")
write.csv(meth14, file="methCpG8desfalse100.csv")


#meth15
mat15=percMethylation(meth15)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat15, file= "matCpG8desfalse90.csv")
write.csv(meth15, file="methCpG8desfalse90.csv")


#meth16
mat16=percMethylation(meth16)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat16, file= "matCpG8desfalse80.csv")
write.csv(meth16, file="methCpG8desfalse80.csv")


#meth17
mat17=percMethylation(meth17)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat17, file= "matCpG5desfalse100.csv")
write.csv(meth17, file="methCpG5desfalse100.csv")

#meth18
mat18=percMethylation(meth18)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat18, file= "matCpG5desfalse90.csv")
write.csv(meth18, file="methCpG5desfalse90.csv")


#meth19
mat19=percMethylation(meth19)
#head(mat)

#dim(matdestrue)
#dim(meth)

write.csv(mat19, file= "matCpG5desfalse80.csv")
write.csv(meth19, file="methCpG5desfalse80.csv")

