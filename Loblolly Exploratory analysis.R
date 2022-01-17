#setting working directory
setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/ALL/")
library(psych)
library(caret)
#check to see how to load in libraries on cluster
library(methylKit)
library(impute)
#starting with 5x 80% files. Looking at everything initially, then running lower quality sites, etc.
#need to initially put in header for all of these, otherwise they wont work. Need to look into making
#code that will do this for me in the future to save on time. 

#CpG
#start with 5x, sites across 80 percent of samples covered
mymat1<-read.csv(file="matCpG5desfalse80.csv", header = T, row.names="ID")
head(mymat1)

#transforming data so you can call variables and it is in the right order
lob1<- as.data.frame(t(mymat1))
#head(lob1)
dim(lob1)
# 24 30351


corrtest1 <- corr.test(lob1[-1], lob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval) ##9 that are significantly correlated
lowestpval<-sort(pval, decreasing = F)
lowestpval
#

decreasingcor<-sort(corrtest1$r, decreasing = F)
decreasingcor
#

increasingcor<-sort(corrtest1$r, decreasing = T)
increasingcor

a<-corrtest1$r
b<-corrtest1$p
correlated<-c(a,b)
View(correlated)

#generating histogram of sites correlated with Age
tiff(file="CpG5x80 age corr wlines des false.tiff", res = 800, width = 120, height = 130, units = "mm")
ZZ<-hist(corrtest1$r, main = "CpGx80percent",xaxt = "n", xlab = "Cytosine Age correlation", col = c("blue"), ylim = range(0:15000))
abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()

#filtering each individual set
#CPG
dim(lob1)
# 24 30351 sites prior to filtering


#getting all correlated sites
a<-subset(corrtest1$r, corrtest1$r > 0.5)
b<-subset(corrtest1$r, corrtest1$r < -0.5)
c<-subset(corrtest1$p.adj, corrtest1$r > 0.5)
d<-subset(corrtest1$p.adj, corrtest1$r < -0.5)


corrsites<-rbind(a,b)
corrsites
length(corrsites)
corrpval<-rbind(c,d)


correlated<-data.frame(corrsites, corrpval)

correlated


write.csv(correlated, "CpGCorrelatedsitesnofilt.csv")


#Filtering zero variant sites

#By default, nearZeroVar will return the positions of the variables that are flagged to be problematic

#removing near zero variant sites
nzv <- nearZeroVar(lob1)
filteredlob1 <- lob1[, -nzv]
dim(filteredlob1)
#now 21566 sites (so removed 8784 sites, or ~29% of the sites)
nzv

write.csv(nzv, "InvariantCpGsites.csv")

#running new correlation with age

corrtest1 <- corr.test(filteredlob1[-1], filteredlob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#10 sites significantly correlated w/ age  (after filtering near zero variance sites)
#previously, CPG des false 5x 80 perc had 9 significant sites, so filtering added 1

decreasingcor<-sort(corrtest1$r, decreasing = F)
decreasingcor
#17 strongly neg correlated (between -0.7 and -0.8) (same)
#70 between -0.6 and -0.7 (same)
#239 between -0.5 and -0.6  So no change in numbers before vs after filtering, just affected P values

increasingcor<-sort(corrtest1$r, decreasing = T)
View(increasingcor)
#2 between 0.8 and 0.9 (same)
#15 between 0.7 and 0.8 (same)
#27 between 0.6 and 0.7 (same)
#163 between 0.5 and 0.6 (same)


a<-corrtest1$r
b<-corrtest1$p.adj
correlated<-data.frame(a,b)
correlated
View(correlated)


#getting total -
a<-subset(corrtest1$r,corrtest1$r>0) 
length(a)
#getting total +
b<-subset(corrtest1$r,corrtest1$r<0) 
length(b)
length(corrtest1$r)
b<-subset(corrtest1$r,corrtest1$r=0) 
length(b)
rval <- corrtest1$r
#View(pval)
lowestrval<-sort(rval, decreasing = F)
lowestrval
write.csv(lowestrval, file = "cpgfilteredrvalues.csv")

#generating histogram of sites correlated with Age
tiff(file="CpG5x80 age corr wlines des false filtered.tiff", res = 800, width = 100, height = 120, units = "mm")
ZZ<-hist(corrtest1$r, main = " ",xaxt = "n", xlab = "Cytosine Age correlation", col = c("blue"), ylim = range(0:15000))
abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()


#now looking at selecting only correlated sites with age

#############
#CpG

#head(corrtest1) 
head(corrtest1$p)


#getting all correlated sites
a<-subset(corrtest1$r, corrtest1$r > 0.5)
b<-subset(corrtest1$r, corrtest1$r < -0.5)
c<-subset(corrtest1$p.adj, corrtest1$r > 0.5)
d<-subset(corrtest1$p.adj, corrtest1$r < -0.5)


corrsites<-rbind(a,b)
corrsites

corrpval<-rbind(c,d)


correlated<-data.frame(corrsites, corrpval)

correlated


write.csv(correlated, "CpGCorrelatedsites.csv")

#manually pulled out the top 35 correlated CpG sites and made a new csv file of these
CpGtop35<-read.csv(file = "CpGCorrelatedsitestop35.csv")

list<-as.matrix(CpGtop35[c(1)])
list
#pulling out the info and coordinates for these top 35 sites
setwd("~/Parrott lab/Research/loblolly pine experiment/CPG/DestrandF/meth files")

CpGmeth<-read.csv(file = "methCpG5desfalse80.csv")

CpGmeth[c(list),][,c(1:3)]

dim(methfile[c(list),][,c(1:3)])

top35cpGtable<-CpGmeth[c(list),][,c(1:3)]
top35cpGtable
write.csv(top35cpGtable, file = "Top35CpG.csv")



#now need to figure out how to get sites (first column) and pull these out from the filteredlob1 dataframe
subset(pval,pval<0.05)[ ,1]
length(subset(pval,pval<0.05)[ ,1])#10 sites so I did this correctly


sig<-subset(pval,pval<0.05)[ ,1]
#sigfilteredlob1 <- filteredlob1[, -sig]
#dim(sigfilteredlob1)


#so the list of age sig variables is a vector
is.vector(sig)

sig2<-as.matrix(sig)

#now it is a matrix, and rownames will give the sites that were found to be significantly correlated with age
rownames(sig2)
length(rownames(sig2)) #returns 10, so good

sig2

#where sites in the all
#CpG : 1-30350, CHG: 30351 - 65428. CHH: 65429 - 146918

#the ones it returns are: 
#10 CpG sites: "site2337"  "site8868"  "site10683" "site13231" "site14054" "site14822" "site17458" "site17859" "site23567" "site24067"

#so same ones large model used


#getting R^2 values of significantly correlated sites
head(corrtest1$r)

subset(corrtest1$r, corrtest1$p<0.05)
subset(pval,pval<0.05)

#making new dataframe with significantly correlated sites with age
r<-subset(corrtest1$r, corrtest1$p<0.05)
p<-subset(pval,pval<0.05)
sigsites<-data.frame(r,p)
sigsites

write.csv(sigsites, file = "CpG5x80percinvfilt.csv")

#Genomic ranges

#trying the annotations using genomic ranges

#reading in for CpG

setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/Annotation file/")

#if (!require("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")


library(GenomicRanges)



# read CpGi data set
#filePath=system.file("Loblolly_cpg_islands.bed")
#filePath


cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)



#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/CPG/DestrandF/meth files/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCpG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)
setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/ALL/")


#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks
dim(filepathpeaks)


#removing the invariant sites to obtain proper background
inv<-read.csv(file = "InvariantCpGsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

#dim(filepathpeaks[c(inv),])

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)

#noe filepathpeals has 21566 sites, the same as the background I was running following inv filtering


filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#2906 of our invariant filtered background sites fall into islands

#2906 / 21566 = 13.47% of our background is in Cpg islands

#the below will tell you if sites are falling in ranges. Many will have zeros
counts1=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts1)

#the find overlaps function tells which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.
findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/CPG/DestrandF/meth files")

CpGmeth<-read.csv(file = "methCpG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")


#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
corrsites<-read.csv(file="CpGCorrelatedsites.csv")
corrsites
dim(corrsites)
View(corrsites)

Corrsites<-corrsites[,1]
Corrsites
View(Corrsites)
View(filepathpeaks)

#use Agesites to pull out sites from filepathpeaks 
#is.data.frame(filepathpeaks)
filepathpeaks<-filepathpeaks[-1, ]
filepathpeaks[246,]
filepathpeaks[c(Corrsites),]
dim(filepathpeaks)
filepathpeaks[c(Corrsites),]
newcorrsites<-filepathpeaks[c(Corrsites),]

dim(newcorrsites)
View(newcorrsites)


sitenames<-newcorrsites[,2]
sitenames
length(sitenames)

startsites<-as.numeric(newcorrsites[,3])
endsites<-as.numeric(newcorrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)


overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps
View(overlaps)

cpgisites.gr[32]
newcorrsites[32,]


newcorrsites[32,]
corrsites[32, ]#first was posotive, makes sense because islands undergo hypermethylation with age

overlaps[]
overlapping<-overlaps@from
overlapping

corrsites[c(overlapping),]

#dim(corrsites[c(overlapping),])
#125

Islandsites<-corrsites[c(overlapping),]

length(Islandsites$corrsites[Islandsites$corrsites>0]) # 14
length(Islandsites$corrsites[Islandsites$corrsites<0]) # 111



#islands of interest:


cpgi.gr[32,]# this corresponds to range () on scaffold 
cpgi.gr[,]# this corresponds to range () on scaffold 
cpgi.gr[,]# this corresponds to range () on scaffold 

cpgi.gr[c(1:125),]


#now trying to get the sites significantly correlated with age: (10 of these), following inv. filtering:
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CpG5x80percinvfilt.csv")
agesites

Agesites<-agesites[,1]
Agesites

#use Agesites to pull out sites from filepathpeaks 
dim(filepathpeaks)
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#2 of our significant sites fell in islands

#the below will tell you if sites are falling in ranges. Many will have zeros
counts3=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts3)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from
overlapping

Agesites
IslandSites<-agesites[c(overlapping),]
IslandSites

Islandsites#age associated
IslandSites#significant
Islandsites
write.csv(Islandsites, file ="CpG Age associated island sites.csv")
write.csv(IslandSites, file ="CpG Age significant island sites.csv")
IslandSites


#Now moving to shores:

setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CPG/DestrandF/meth files")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCpG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks
# read the peaks from meth file

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "InvariantCpGsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

#dim(filepathpeaks[c(inv),])

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)

filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)


#this is the one we are interested in
findOverlaps(cpgisites.gr,cpgi.gr)

cpgi.gr[199,] 

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/CPG/DestrandF/meth files")

CpGmeth<-read.csv(file = "methCpG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1, ]
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CpGCorrelatedsites.csv")
corrsites

Corrsites<-corrsites[,1]
Corrsites

#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[2360,]
Corrsites<-filepathpeaks[c(Corrsites),]

Corrsites
dim(Corrsites)

newagesites[,2]

sitenames<-Corrsites[,2]
sitenames

startsites<-as.numeric(Corrsites[,3])
endsites<-as.numeric(Corrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#91 sites that were correlated fell into chelves. So, 91 / 533 = 17.07% falling into shores


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps

overlapping<-overlaps@from

dim(corrsites)
Shoresites<-corrsites[c(overlapping),]
Shoresites
length(Shoresites$corrsites[Shoresites$corrsites>0]) # 35
length(Shoresites$corrsites[Shoresites$corrsites<0]) # 56

write.csv(Shoresites, file="CpGShore age associated sites.csv")

#now trying to get the sites correlated with age: (70 of these), following inv. filtering:
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CpG5x80percinvfilt.csv")
agesites

Agesites<-agesites[,1]
Agesites

#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[2360,]
newagesites<-filepathpeaks[c(Agesites),]

newagesites

newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)


overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

ShoreSites<-agesites[c(overlapping),]
ShoreSites

write.csv(ShoreSites, file = "CpGShore significant age sites.csv")

cpgisites.gr

#now for shelves:
setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")


cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CPG/DestrandF/meth files")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCpG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)


filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

# read the peaks from meth file

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "InvariantCpGsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

#dim(filepathpeaks[c(inv),])

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)


filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)

findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/CPG/DestrandF/meth files")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCpG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)


#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CpGCorrelatedsites.csv")
agesites

Agesites<-agesites[,1]
Agesites

#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[2360,]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#16 correlated sites fell into shelves

#so 16 / 533 = 3% of correlated sites fall into shores 


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(agesites)

Shelfsites<-agesites[c(overlapping),]
Shelfsites
length(Shelfsites$corrsites[Shelfsites$corrsites>0]) # 9
length(Shelfsites$corrsites[Shelfsites$corrsites<0]) # 7
write.csv(Shelfsites, file = "CpGShelf age associated sites.csv")


#now trying to get the sites correlated with age: (10 of these), following inv. filtering:
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CpG5x80percinvfilt.csv")
agesites

Agesites<-agesites[,1]
Agesites

#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[2360,]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)


overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps

#SO, 0 of my 10 age-correlated sites fell into shelves!

###Graphing results

#Now that I have the results:

#CpG islands: 651,745 islands that exist across the loblolly genome
#background: 5027 / 30350 = 16.56% of the CpG sites we have are falling into islands 
#age associated: 3/10 = 30%


#CpG Shores: 1,144,719 shores that exist across the loblolly genome (in the scaffolds
#background: 4427 / 30350 = 14.58% of the sites we have fell into shores
#age-associated: 2/10 = 20%


#CPG shelves: 700,544 shelves that exist across the loblolly genome
#background: 1816 / 30350 =  0.0598 or 5.98% 
#age-associated: 1 / 10 = 10% 


#Open sea (by process of elimination)
#Genome:?
#Background: (Islands + shores = shelves = 11,270. 30350 - 11270 = 19080.) So, 19080 / 30350 = 0.6287 or 62.87%
#so 62.87% fall into open seas
#Age-associated: 4 / 10 = 40%. So 40% of the age-associated sites fall in open seas. 


#getting pos and neg from age associated for open sea
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

ageislands<-read.csv(file="CpG Age associated island sites")
ageislands
ageshores<-read.csv(file = "CpGShore age associated sites.csv")
ageshelves<-read.csv(file="CpGShelf age associated sites.csv")
ageshores
ageshelves

Correlatedsites<-read.csv(file = "CpGCorrelatedsites.csv")
dim(Correlatedsites)
head(Correlatedsites)
head(ageislands)

ageislands<-ageislands[c(2,3,4)]
ageshores<-ageshores[c(2,3,4)]
ageshelves<-ageshelves[c(2,3,4)]

head(ageislands)

island<-ageislands[, 1]
island

head(Correlatedsites)
dim(Correlatedsites)


island
subset(Correlatedsites, (Site==246))
list<-setdiff(Correlatedsites$Site, ageislands$Site)   # elements in a$x NOT in b$y
list

length(list)

dim(Correlatedsites[Correlatedsites$Site %in% list, ])

Correlatedminusisland<-Correlatedsites[Correlatedsites$Site %in% list, ]
dim(Correlatedminusisland)


#removing shores
list<-setdiff(Correlatedminusisland$Site, ageshores$Site)   # elements in a$x NOT in b$y
list

length(list)# going from 408 to 317, means 91 sites coming out

dim(Correlatedminusisland[Correlatedminusisland$Site %in% list, ])

Correlatedminusislandandshore<-Correlatedminusisland[Correlatedminusisland$Site %in% list, ]
dim(Correlatedminusislandandshore)

#removing shelves
list<-setdiff(Correlatedminusislandandshore$Site, ageshelves$Site)   # elements in a$x NOT in b$y
list

length(list)# going from 317 to 301, means  16 sites coming out

dim(Correlatedminusislandandshore[Correlatedminusislandandshore$Site %in% list, ])

Correlatedminusislandandshoreandshelf<-Correlatedminusislandandshore[Correlatedminusislandandshore$Site %in% list, ]
dim(Correlatedminusislandandshoreandshelf)

#now that all other sites removed, look at pos and neg values for open sea regions!

length(Correlatedminusislandandshoreandshelf$corrsites[Correlatedminusislandandshoreandshelf$corrsites>0]) # 149
#149 open sea sites that are positive with age 
length(Correlatedminusislandandshoreandshelf$corrsites[Correlatedminusislandandshoreandshelf$corrsites<0]) # 152
#152 CpG sites falling in open seas that are negative with age

write.csv(Correlatedminusislandandshoreandshelf, file = "CpG age associated open sea sites.csv")

#significant sites
ageislands<-read.csv(file="CpG Age significant island sites")
ageislands
ageshores<-read.csv(file = "CpGShore significant age sites.csv")
#ageshelves<-read.csv(file="CpGshelf age significant sites.csv")
ageshores
ageshelves

Correlatedsites<-read.csv(file = "CpG5x80percinvfilt.csv")
Correlatedsites
Correlatedsites
head(ageislands)

ageislands<-ageislands[c(2,3,4)]
ageshores<-ageshores[c(2,3,4)]

head(ageislands)

island<-ageislands[, 1]
island

head(Correlatedsites)
dim(Correlatedsites)

Correlatedsites$X

island
subset(Correlatedsites, (Site==246))
list<-setdiff(Correlatedsites$X, ageislands$X)   # elements in a$x NOT in b$y
list
Correlatedsites$X
length(list)

dim(Correlatedsites[Correlatedsites$X %in% list, ])

Correlatedminusisland<-Correlatedsites[Correlatedsites$X %in% list, ]
dim(Correlatedminusisland)

#removing shores
list<-setdiff(Correlatedminusisland$X, ageshores$X)   # elements in a$x NOT in b$y
list

length(list)# going from 408 to 317, means 91 sites coming out

dim(Correlatedminusisland[Correlatedminusisland$X %in% list, ])

Correlatedminusislandandshore<-Correlatedminusisland[Correlatedminusisland$X %in% list, ]
dim(Correlatedminusislandandshore)
Correlatedminusislandandshore

write.csv(Correlatedminusislandandshore, file = "CpG Sig open sea sites.csv")

#now that all other sites removed, look at pos and neg values for open sea regions!

#Graphing CpG genomic distribitions 

library(ggplot2)
#vs ggplot
Region <- c(rep("Islands" , 2) , rep("Shores" , 2) , rep("Shelves" , 2) , rep("Open Seas" , 2) )
Sites <- rep(c("Background" , "Age-Associated") , 4)
value <- c(13.47, 23.45, 13.86,17.07, 6.32, 3, 66.34, 56.47)
data <- data.frame(Region,Sites,value)
data


#setting the order for how I want it to show up in the plot
data$Region <- factor(data$Region, levels = c("Islands", "Shores", "Shelves", "Open Seas"))
data$Sites <- factor(data$Sites, levels = c("Background", "Age-Associated"))
data$Sites
data$Region

# Grouped
annotation<-ggplot(data, aes(fill=Sites, y=value, x=Region)) + 
  geom_bar(position="dodge", size = 0.5, color = "black", stat="identity")+ scale_fill_manual(values = c("black", "blue"))

annotation

#cytosine<-ggplot(data, aes(fill=context, y=value, x=status)) + 
#  geom_bar(position="stack", size = 0.5, color = "black", stat="identity") + scale_fill_manual(values = c("red","green","Blue"))

annotation<-annotation + ylim(range(c(0,100)))

annotation


annotation<-annotation +labs(y= "Sequenced CpG Methylation Sites (%)", x = "Genomic Context")
annotation
annotation<-annotation + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))
annotation


annotation<-annotation + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
annotation


annotation<-annotation + scale_y_continuous(limits = c(0,100), expand = c(0, 0))
annotation
annotation<-annotation + theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=14))


annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
annotation<-annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
annotation
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))

#tiff(file="CHH Genomic context.tiff", res = 800, width = 189, height = 100, units = "mm")
tiff(file="CpG Genomic context no sig.tiff", res = 800, width = 189, height = 100, units = "mm")
annotation
dev.off()



#Binomial tests for enrichment of age-associated sites in genomic ranges

#CpG 

##island (age associated vs background)
binom.test(125, 533, p = 0.1347,
           alternative = c("two.sided"),
           conf.level = 0.95) ##enriched in islands

#shore
binom.test(91, 533, p = 0.1386,
           alternative = c("two.sided"),
           conf.level = 0.95) ##enriched in shores

#shelf
binom.test(16, 533, p = 0.063,
           alternative = c("two.sided"),
           conf.level = 0.95) ##depleted in shelves


#open sea
binom.test(301, 533, p = 0.6634,
           alternative = c("two.sided"),
           conf.level = 0.95) ##depleted in open sea

#so, Age-associated CPG cytosines were enriched in islands and shores and depleted in shelves and open seas


#significant sites
##island (age associated vs background)
binom.test(2, 10, p = 0.1347,
           alternative = c("two.sided"),
           conf.level = 0.95) ##enriched in islands

#shore
binom.test(3, 10, p = 0.1386,
           alternative = c("two.sided"),
           conf.level = 0.95) ##enriched in shores

#shelf
binom.test(0, 10, p = 0.063,
           alternative = c("two.sided"),
           conf.level = 0.95) ##depleted in shelves


#open sea
binom.test(5, 10, p = 0.6634,
           alternative = c("two.sided"),
           conf.level = 0.95) ##depleted in open sea


#no sig diff in sig cytosines


#CHH

####CHG
#start with 5x, sites across 80 percent of samples covered
#setting working directory
setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/ALL/")
mymat1<-read.csv(file="matCHG5desfalse80.csv", header = T, row.names="ID")
head(mymat1)


#transforming data so you can call variables and it is in the right order
lob1<- as.data.frame(t(mymat1))
#head(lob1)
dim(lob1)
# 24 35079 CpG sites


library(psych)

corrtest1 <- corr.test(lob1[-1], lob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval) ##63 that are significantly correlated
lowestpval<-sort(pval, decreasing = F)
lowestpval
#

decreasingcor<-sort(corrtest1$r, decreasing = F)
View(decreasingcor) 
#

increasingcor<-sort(corrtest1$r, decreasing = T)
View(increasingcor)
#

#generating histogram of sites correlated with Age
tiff(file="CHG5x80 age corr wlines des false.tiff", res = 800, width = 120, height = 130, units = "mm")
ZZ<-hist(corrtest1$r, main = "CHG5x80percent",xaxt = "n", xlab = "Cytosine Age correlation", col = c("green"), ylim = range(0:15000))
abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()


#CHG
#start with the 5x 80 percent ALL file

dim(lob1)
# 24 35079 sites prior to filtering

#identify near zero-variance variables

#By default, nearZeroVar will return the positions of the variables that are flagged to be problematic
library(caret)
#removing near zero variant sites
nzv <- nearZeroVar(lob1)
filteredlob1 <- lob1[, -nzv]
dim(filteredlob1)
#now 25501 sites (so removed 9577 sites, or ~27% of the sites)

nzv

write.csv(nzv, "CHGinvariantsites.csv")
#running new correlation with age

corrtest1 <- corr.test(filteredlob1[-1], filteredlob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#70 sites significantly correlated w/ age  (after filtering near zero variance sites)
#previously, CHG des false 5x 80 perc had 63 significant sites, so filtering added 7 (so CHG was more affected by filtering)

decreasingcor<-sort(corrtest1$r, decreasing = F)
View(decreasingcor) 
#Comparing with excel file, no change in numbers before vs after filtering, just affected P values

increasingcor<-sort(corrtest1$r, decreasing = T)
View(increasingcor)
#

a<-subset(corrtest1$r,corrtest1$r>0) 
length(a)
#getting total +
b<-subset(corrtest1$r,corrtest1$r<0) 
length(b)


#generating histogram of sites correlated with Age
tiff(file="CHG5x80 age corr wlines des false filtered.tiff", res = 800, width = 100, height = 120, units = "mm")
ZZ<-hist(corrtest1$r, main = " ",xaxt = "n", xlab = "Cytosine Age correlation", col = c("green"), ylim = range(0:15000))
abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()

#getting all correlated sites
a<-subset(corrtest1$r, corrtest1$r > 0.5)
b<-subset(corrtest1$r, corrtest1$r < -0.5)
c<-subset(corrtest1$p.adj, corrtest1$r > 0.5)
d<-subset(corrtest1$p.adj, corrtest1$r < -0.5)

corrsites<-rbind(a,b)
corrsites

corrpval<-rbind(c,d)

correlated<-data.frame(corrsites, corrpval)

correlated

write.csv(correlated, "CHGCorrelatedsites.csv")

#manually made a new csv file of the top 35 sites in excel
CHGtop35<-read.csv(file = "CHGCorrelatedsitestop35.csv")
list<-as.matrix(CHGtop35[c(1)])
list
setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth")

CHGmeth<-read.csv(file = "methCHG5desfalse80.csv")

CHGmeth[c(list),][,c(1:3)]

dim(methfile[c(list),][,c(1:3)])

top35cHGtable<-CHGmeth[c(list),][,c(1:3)]
top35cHGtable
write.csv(top35cHGtable, file = "Top35CHG.csv")

dim(filteredlob1)
#now 25502 sites 

#assessing genomic ranges for CHG sites 


setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")


library(GenomicRanges)

#island bed file
cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)


#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)


#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks


setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "CHGinvariantsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

#dim(filepathpeaks[c(inv),])

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)


filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts1=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts1)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.

findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks
dim(filepathpeaks)

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CHGCorrelatedsites.csv")#need CHG correlated sites
dim(corrsites)

Corrsites<-corrsites[,1]
length(Corrsites)

#use Agesites to pull out sites from filepathpeaks 
#filepathpeaks[2360,]
newcorrsites<-filepathpeaks[c(Corrsites),]
#length(filepathpeaks)
#newcorrsites<-filepathpeaks[c(Corrsites),]
#length(newcorrsites)
dim(newcorrsites)

#newagesites[,2]

sitenames<-newcorrsites[,2]
sitenames

startsites<-as.numeric(newcorrsites[,3])
endsites<-as.numeric(newcorrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps[1:125]

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps
View(overlaps)

corrsites[c(overlapping),]

Islandsites<-corrsites[c(overlapping),]

length(Islandsites$corrsites[Islandsites$corrsites>0]) # 24
length(Islandsites$corrsites[Islandsites$corrsites<0]) # 105


write.csv(Islandsites, file = "CHG Age associated Island sites.csv")

#now get the sites significantly correlated with age: (10 of these), following inv. filtering:
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CHG5x80percinvfilt.csv")#need significantly correlated sites
dim(agesites)
agesites


Agesites<-agesites[,1]
Agesites

newagesites<-filepathpeaks[c(Agesites),]


dim(newagesites)


newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


counts3=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts3)


findOverlaps(cpgisites.gr,cpgi.gr)


overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps
View(overlaps)


overlapping<-overlaps@from
overlapping

corrsites[c(overlapping),]

Islandsites<-agesites[c(overlapping),]
dim(Islandsites)
Islandsites

length(Islandsites$r[Islandsites$r>0]) # 0
length(Islandsites$r[Islandsites$r<0]) # 21


write.csv(Islandsites, file = "CHG Age Sig Island sites.csv")

#Now moving to shores:

setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)


#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks
# read the peaks from meth file

dim(filepathpeaks)

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "CHGinvariantsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

#filepathpeaks<-(filepathpeaks[-c(inv),])
filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)


filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)

findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CHGCorrelatedsites.csv")#need chg correlated sites
corrsites

Corrsites<-corrsites[,1]
Corrsites

setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

#use Agesites to pull out sites from filepathpeaks 

Corrsites<-filepathpeaks[c(Corrsites),]

Corrsites

newagesites[,2]

sitenames<-Corrsites[,2]
sitenames

startsites<-as.numeric(Corrsites[,3])
endsites<-as.numeric(Corrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

findOverlaps(cpgisites.gr,cpgi.gr)

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps
View(overlaps)

overlapping<-overlaps@from
overlapping

corrsites[c(overlapping),]
dim(corrsites[c(overlapping),])
#125

Shoresites<-corrsites[c(overlapping),]

length(Shoresites$corrsites[Shoresites$corrsites>0]) # 24
length(Shoresites$corrsites[Shoresites$corrsites<0]) # 105


write.csv(Shoresites, file = "CHG Age associated Shore sites.csv")

#now trying to get the sites correlated with age: (70 of these), following inv. filtering:
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CHG5x80percinvfilt.csv")#need significant correlated sites
agesites

Agesites<-agesites[,1]
Agesites

#use Agesites to pull out sites from filepathpeaks 
dim(filepathpeaks)
newagesites<-filepathpeaks[c(Agesites),]

newagesites

newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#6 of our sig sites fall into shores, so 6 / 70 = 

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)


findOverlaps(cpgisites.gr,cpgi.gr)

cpgisites.gr

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps
View(overlaps)


overlaps[]
overlapping<-overlaps@from
overlapping

agesites[c(overlapping),]


Shoresites<-agesites[c(overlapping),]
Shoresites

length(Shoresites$r[Shoresites$r>0]) # 1
length(Shoresites$r[Shoresites$r<0]) # 5


write.csv(Shoresites, file = "CHG Age sig Shore sites.csv")

#now for shelves:
setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)


#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
dim(filepathpeaks)
# read the peaks from meth file


dim(filepathpeaks)

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "CHGinvariantsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)


filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)

#this is the one we are interested in
findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CHGCorrelatedsites.csv")#need chg correlated sites
corrsites

Corrsites<-corrsites[,1]
Corrsites

setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

#use Agesites to pull out sites from filepathpeaks 
Corrsites<-filepathpeaks[c(Corrsites),]

Corrsites

dim(Corrsites)
newagesites[,2]

sitenames<-Corrsites[,2]
sitenames

startsites<-as.numeric(Corrsites[,3])
endsites<-as.numeric(Corrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)


overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps


overlappingsites<-overlaps@from
dim(corrsites)

dim(corrsites[c(overlappingsites),])

Shelfsites<-corrsites[c(overlappingsites),]
Shelfsites

length(Shelfsites$corrsites[Shelfsites$corrsites>0]) # 1
length(Shelfsites$corrsites[Shelfsites$corrsites<0]) # 5

write.csv(Shelfsites, file = "CHG age associated shelf sites.csv")

#now trying to get the sites correlated with age: (70 of these), following inv. filtering:
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

agesites<-read.csv(file="CHG5x80percinvfilt.csv")#need significant correlated sites
agesites

Agesites<-agesites[,1]
Agesites

#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[2360,]
newagesites<-filepathpeaks[c(Agesites),]

newagesites

newagesites[,2]

sitenames<-newagesites[,2]
sitenames

startsites<-as.numeric(newagesites[,3])
endsites<-as.numeric(newagesites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)

overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps


overlappingsites<-overlaps@from
dim(corrsites)

dim(agesites[c(overlappingsites),])

Shelfsites<-agesites[c(overlappingsites),]
Shelfsites

write.csv(Shelfsites, file = "CHG age sig shelf sites.csv")

#getting pos and neg sites for open sea regions
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

ageislands<-read.csv(file="CHG Age associated Island sites.csv")
ageislands
ageshores<-read.csv(file = "CHG Age associated Shore sites.csv")
ageshelves<-read.csv(file="CHG age associated shelf sites.csv")
ageshores
ageshelves

Correlatedsites<-read.csv(file = "CHGCorrelatedsites.csv")
dim(Correlatedsites)
head(Correlatedsites)
head(ageislands)

ageislands<-ageislands[c(2,3,4)]
ageshores<-ageshores[c(2,3,4)]
ageshelves<-ageshelves[c(2,3,4)]
ageshores
ageislands
ageshelves
head(ageislands)

island<-ageislands[, 1]
island

head(Correlatedsites)
dim(Correlatedsites)

island
subset(Correlatedsites, (Site==246))
list<-setdiff(Correlatedsites$X, ageislands$X)   # elements in a$x NOT in b$y
list

length(list)

dim(Correlatedsites[Correlatedsites$X %in% list, ])

Correlatedminusisland<-Correlatedsites[Correlatedsites$X %in% list, ]
dim(Correlatedminusisland)

Correlatedminusisland
#removing shores
list<-setdiff(Correlatedminusisland$X, ageshores$X)   # elements in a$x NOT in b$y
list

list<-setdiff(Correlatedminusisland$X, ageshelves$X)   # elements in a$x NOT in b$y
list

length(list)# going from 740 to 680, means  60 sites coming out

dim(Correlatedminusisland[Correlatedminusisland$X %in% list, ])

Correlatedminusislandandshelf<-Correlatedminusisland[Correlatedminusisland$X %in% list, ]
dim(Correlatedminusislandandshelf)

list<-setdiff(Correlatedminusislandandshelf$X, ageshores$X)   # elements in a$x NOT in b$y
list

length(list) #went from 680 to 594, so 86 sites came out. 
View(ageshores)


write.csv(Correlatedminusislandandshelf, file = "CHG age associated minus islands and shelves.csv")

Correlatedminusislandandshore<-Correlatedminusisland[Correlatedminusisland$X %in% list, ]
dim(Correlatedminusislandandshore)

#removing shelves
list<-setdiff(Correlatedminusislandandshore$X, ageshelves$X)   # elements in a$x NOT in b$y
list

length(list)# going from 317 to 301, means  16 sites coming out

dim(Correlatedminusislandandshore[Correlatedminusislandandshore$X %in% list, ])

Correlatedminusislandandshoreandshelf<-Correlatedminusislandandshore[Correlatedminusislandandshore$X %in% list, ]
dim(Correlatedminusislandandshoreandshelf)

#now that all other sites removed, look at pos and neg values for open sea regions!

length(Correlatedminusislandandshoreandshelf$corrsites[Correlatedminusislandandshoreandshelf$corrsites>0]) # 149
#149 open sea sites that are positive with age 
length(Correlatedminusislandandshoreandshelf$corrsites[Correlatedminusislandandshoreandshelf$corrsites<0]) # 152
#152 CpG sites falling in open seas that are negative with age

write.csv(Correlatedminusislandandshoreandshelf, file = "CHG age associated open sea sites.csv")

#significant sites
ageislands<-read.csv(file="CHG Age Sig Island sites.csv")
ageislands
ageshores<-read.csv(file = "CHG Age sig Shore sites.csv")
ageshelves<-read.csv(file="CHG age sig shelf sites.csv")
ageshores
ageshelves

Correlatedsites<-read.csv(file = "CHG5x80percinvfilt.csv")
Correlatedsites
Correlatedsites
ageislands

ageislands<-ageislands[c(2,3,4)]
ageshores<-ageshores[c(2,3,4)]
ageshelves<-ageshelves[c(2,3,4)]
ageshelves
head(ageislands)

island<-ageislands[, 1]
island

head(Correlatedsites)
dim(Correlatedsites)

#Correlatedsites[Correlatedsites$Site]

#dim(Correlatedsites[ , -c(island)])
Correlatedsites$X
#Correlatedsites$Site!=246
island
subset(Correlatedsites, (Site==246))
list<-setdiff(Correlatedsites$X, ageislands$X)   # elements in a$x NOT in b$y
list
Correlatedsites$X
length(list)

dim(Correlatedsites[Correlatedsites$X %in% list, ])

Correlatedminusisland<-Correlatedsites[Correlatedsites$X %in% list, ]
dim(Correlatedminusisland)

#removing shores
list<-setdiff(Correlatedminusisland$X, ageshores$X)   # elements in a$x NOT in b$y
list

length(list)# going from 408 to 317, means 91 sites coming out

dim(Correlatedminusisland[Correlatedminusisland$X %in% list, ])

Correlatedminusislandandshore<-Correlatedminusisland[Correlatedminusisland$X %in% list, ]
dim(Correlatedminusislandandshore)
Correlatedminusislandandshore

#removinog shelves
list<-setdiff(Correlatedminusislandandshore$X, ageshelves$X)   # elements in a$x NOT in b$y
list
Correlatedsites$X
length(list)

dim(Correlatedminusislandandshore[Correlatedminusislandandshore$X %in% list, ])

Correlatedminusislandandshoreandshelf<-Correlatedminusislandandshore[Correlatedminusislandandshore$X %in% list, ]
dim(Correlatedminusislandandshoreandshelf)

write.csv(Correlatedminusislandandshoreandshelf, file = "CHG sig age open sea sites.csv")
#Correlatedminusislandandshoreandshelf$
length(Correlatedminusislandandshoreandshelf$r[Correlatedminusislandandshoreandshelf$r>0]) # 149
length(Correlatedminusislandandshoreandshelf$r[Correlatedminusislandandshoreandshelf$r<0]) # 149

#now that all other sites removed, look at pos and neg values for open sea regions!


###Graphing results

#Now that I have the results:

#CpG islands: 651,745 islands that exist across the loblolly genome
#background: 5027 / 30350 = 16.56% of the CpG sites we have are falling into islands 
#age associated: 3/10 = 30%




#CpG Shores: 1,144,719 shores that exist across the loblolly genome (in the scaffolds
#background: 4427 / 30350 = 14.58% of the sites we have fell into shores
#age-associated: 2/10 = 20%



#CPG shelves: 700,544 shelves that exist across the loblolly genome
#background: 1816 / 30350 =  0.0598 or 5.98% 
#age-associated: 1 / 10 = 10% 


#Open sea (by process of elimination)
#Genome:?
#Background: (Islands + shores = shelves = 11,270. 30350 - 11270 = 19080.) So, 19080 / 30350 = 0.6287 or 62.87%
#so 62.87% fall into open seas
#Age-associated: 4 / 10 = 40%. So 40% of the age-associated sites fall in open seas. 

#now need to make graph 
#vs ggplot
Region <- c(rep("Islands" , 2) , rep("Shores" , 2) , rep("Shelves" , 2) , rep("Open Seas" , 2) )
Sites <- rep(c("Background" , "Age-Associated") , 4)
value <- c(10.35, 14.84, 13.68, 10.01, 7.02, 6.9, 68.95, 68.24)
data <- data.frame(Region,Sites,value)
data


#setting the order for how I want it to show up in the plot
data$Region <- factor(data$Region, levels = c("Islands", "Shores", "Shelves", "Open Seas"))
data$Sites <- factor(data$Sites, levels = c("Background", "Age-Associated"))
data$Sites
data$Region

# Grouped
annotation<-ggplot(data, aes(fill=Sites, y=value, x=Region)) + 
  geom_bar(position="dodge", size = 0.5, color = "black", stat="identity")+ scale_fill_manual(values = c("black","green"))
dev.off()
annotation

#cytosine<-ggplot(data, aes(fill=context, y=value, x=status)) + 
#  geom_bar(position="stack", size = 0.5, color = "black", stat="identity") + scale_fill_manual(values = c("red","green","Blue"))

annotation<-annotation + ylim(range(c(0,100)))

annotation


annotation<-annotation +labs(y= "Sequenced CHG Methylation Sites (%)", x = "Genomic Context")
annotation

annotation<-annotation + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))
annotation


annotation<-annotation + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 10, colour = "black"))
annotation


annotation<-annotation + scale_y_continuous(limits = c(0,100), expand = c(0, 0))
annotation

annotation<-annotation + theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=14))

annotation

#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))


#tiff(file="CHH Genomic context.tiff", res = 800, width = 189, height = 100, units = "mm")


tiff(file="CHG Genomic context.tiff", res = 800, width = 189, height = 100, units = "mm")
annotation
dev.off()



#Binomial tests for enrichment of age-associated sites in genomic ranges

#CHG


##island (age associated vs background)
binom.test(129, 870, p = 0.1035,
           alternative = c("two.sided"),
           conf.level = 0.95) ##enriched in islands

#shore
binom.test(87, 870, p = 0.1368,
           alternative = c("two.sided"),
           conf.level = 0.95) ##depleted in shores

#shelf
binom.test(60, 870, p = 0.0702,
           alternative = c("two.sided"),
           conf.level = 0.95) ##similar in shelves


#open sea
binom.test(594, 870, p = 0.6895,
           alternative = c("two.sided"),
           conf.level = 0.95) ##sim in open sea


#so, Age-associated CHG cytosines were enriched in islands and depleted in shores


#significant sites
##island (age associated vs background)
binom.test(21, 70, p = 0.1035,
           alternative = c("two.sided"),
           conf.level = 0.95) ##enriched in islands

#shore
binom.test(6, 70, p = 0.1368,
           alternative = c("two.sided"),
           conf.level = 0.95) ##sim in shores

#shelf
binom.test(1, 70, p = 0.0702,
           alternative = c("two.sided"),
           conf.level = 0.95) ##sim in shelves


#open sea
binom.test(42, 70, p = 0.6895,
           alternative = c("two.sided"),
           conf.level = 0.95) ##sim in open sea


# sig cytosines enriched in islands


#####################################
#CHH
#start with 5x, sites across 80 percent of samples covered
#setting working directory
setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/ALL/")
mymat1<-read.csv(file="matCHH5desfalse80.csv", header = T, row.names="ID")
head(mymat1)


#transforming data so you can call variables and it is in the right order
lob1<- as.data.frame(t(mymat1))
#head(lob1)
dim(lob1)
# 24 81491


library(psych)

corrtest1 <- corr.test(lob1[-1], lob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval) ##5 that are significantly correlated
lowestpval<-sort(pval, decreasing = F)
lowestpval
#0 significant

decreasingcor<-sort(corrtest1$r, decreasing = F)
View(decreasingcor) 
#

increasingcor<-sort(corrtest1$r, decreasing = T)
View(increasingcor)
#

#generating histogram of sites correlated with Age
tiff(file="CHH5x80 age corr wlines des false.tiff", res = 800, width = 120, height = 130, units = "mm")
ZZ<-hist(corrtest1$r, main = "CHH5x80percent",xaxt = "n", xlab = "Cytosine Age correlation", col = c("red"), ylim = range(0:15000))
abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()

#start with the 5x 80 percent ALL file
dim(lob1)
# 24 81491 sites prior to filtering

#identify near zero-variance variables

#By default, nearZeroVar will return the positions of the variables that are flagged to be problematic

#removing near zero variant sites
nzv <- nearZeroVar(lob1)
filteredlob1 <- lob1[, -nzv]
dim(filteredlob1)
#now  33152 sites (so removed 48,339 sites, or ~59% of the sites) (far more affected compared to CpG or CHH)

write.csv(nzv, "CHHinvariantsites.csv")

#running new correlation with age

corrtest1 <- corr.test(filteredlob1[-1], filteredlob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#0 sites significantly correlated w/ age  (after filtering near zero variance sites)
#No change from before filtering 

decreasingcor<-sort(corrtest1$r, decreasing = F)
View(decreasingcor) 
# no change in numbers before vs after filtering, just affected P values

increasingcor<-sort(corrtest1$r, decreasing = T)
View(increasingcor)
#
a<-subset(corrtest1$r, corrtest1$r>0)
length(a)
b<-subset(corrtest1$r, corrtest1$r<0)
length(b)

#generating histogram of sites correlated with Age
tiff(file="CHH5x80 age corr wlines des false filtered.tiff", res = 800, width = 100, height = 120, units = "mm")
ZZ<-hist(corrtest1$r, main = " ",xaxt = "n", xlab = "Cytosine Age correlation", col = c("red"), ylim = range(0:15000))
abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()


#getting all correlated sites
a<-subset(corrtest1$r, corrtest1$r > 0.5)
b<-subset(corrtest1$r, corrtest1$r < -0.5)
c<-subset(corrtest1$p.adj, corrtest1$r > 0.5)
d<-subset(corrtest1$p.adj, corrtest1$r < -0.5)


corrsites<-rbind(a,b)
corrsites

corrpval<-rbind(c,d)


correlated<-data.frame(corrsites, corrpval)

correlated

write.csv(correlated, "CHHCorrelatedsites.csv")

#Assessing genomic ranges for CHH sites


setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

library(GenomicRanges)

cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CHH/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHH5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks
# read the peaks from meth file


setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "CHHinvariantsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)

filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts1=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts1)


findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/CHH/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHH5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)


#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks


setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CHHCorrelatedsites.csv")#need CHG correlated sites
dim(corrsites)
corrsites

Corrsites<-corrsites[,1]
length(Corrsites)

#use Agesites to pull out sites from filepathpeaks 
#filepathpeaks[2360,]
newcorrsites<-filepathpeaks[c(Corrsites),]

dim(newcorrsites)


sitenames<-newcorrsites[,2]
sitenames

startsites<-as.numeric(newcorrsites[,3])
endsites<-as.numeric(newcorrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)


overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlaps


overlapping<-overlaps@from
overlapping

dim(corrsites[c(overlapping),])

Islandsites<-corrsites[c(overlapping),]

length(Islandsites$corrsites[Islandsites$corrsites>0])
length(Islandsites$corrsites[Islandsites$corrsites<0])

write.csv(Islandsites, file = "CHH age associated island sites.csv")


#Now moving to shores:

setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")


cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)


#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CHH/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHH5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks
# read the peaks from meth file

dim(filepathpeaks)

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "CHHinvariantsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)

filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)

#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)

#this is the one we are interested in
findOverlaps(cpgisites.gr,cpgi.gr)

#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CHHCorrelatedsites.csv")#need chg correlated sites
corrsites

Corrsites<-corrsites[,1]
Corrsites


setwd("~/Parrott lab/Research/loblolly pine experiment/CHH/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHH5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)


#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

#use Agesites to pull out sites from filepathpeaks 
#filepathpeaks[2360,]
Corrsites<-filepathpeaks[c(Corrsites),]

Corrsites

newagesites[,2]

sitenames<-Corrsites[,2]
sitenames

startsites<-as.numeric(Corrsites[,3])
endsites<-as.numeric(Corrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr


# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)


overlapping<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping

overlapping<-overlapping@from

Shoresites<-corrsites[c(overlapping),]
Shoresites

length(Shoresites$corrsites[Shoresites$corrsites>0])
length(Shoresites$corrsites[Shoresites$corrsites<0])


write.csv(Shoresites, file = "CHH Shore age associated sites.csv")

#now for shelves:
setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

#storing scaffold names
scaffoldnames<-cpgi.df[,1]
scaffoldnames

#getting start and end sections
cpgi.df[,2]#start
cpgi.df[,3]#end

#storing as genomic range
cpgi.gr=GRanges(seqnames = scaffoldnames,
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

cpgi.gr

setwd("~/Parrott lab/Research/loblolly pine experiment/CHH/meth")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHH5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
dim(filepathpeaks)
# read the peaks from meth file

dim(filepathpeaks)

setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

inv<-read.csv(file = "CHHinvariantsites.csv")

inv<-inv[,2]
head(inv)

head(filepathpeaks)

filepathpeaks<-filepathpeaks[-c(inv),]
filepathpeaks  
dim(filepathpeaks)


filepathpeaks[,2]

sitenames<-filepathpeaks[,2]
sitenames

startsites<-as.numeric(filepathpeaks[,3])
endsites<-as.numeric(filepathpeaks[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts)

#this is the one we are interested in
findOverlaps(cpgisites.gr,cpgi.gr)


#getting sites correlated (absolute values of correlation coefficients)
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

corrsites<-read.csv(file="CHHCorrelatedsites.csv")#need chg correlated sites
corrsites
dim(corrsites)
Corrsites<-corrsites[,1]
Corrsites


setwd("~/Parrott lab/Research/loblolly pine experiment/CHH/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHH5desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

#use Agesites to pull out sites from filepathpeaks 
#filepathpeaks[2360,]
Corrsites<-filepathpeaks[c(Corrsites),]
Corrsites

dim(Corrsites)
newagesites[,2]

sitenames<-Corrsites[,2]
sitenames

startsites<-as.numeric(Corrsites[,3])
endsites<-as.numeric(Corrsites[,4])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)


overlapping<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping


overlappingsites<-overlapping@from

Shelfsites<-corrsites[c(overlappingsites),]
Shelfsites

length(Shelfsites$corrsites[Shelfsites$corrsites>0])
length(Shelfsites$corrsites[Shelfsites$corrsites<0])

write.csv(Shelfsites, file = "CHH shelf age associated sites.csv")

#Open sea regions
#getting age associated open sea pos and neg:

#getting pos and neg sites for open sea regions
setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")

ageislands<-read.csv(file="CHH age associated island sites.csv")
ageislands
ageshores<-read.csv(file = "CHH Shore age associated sites.csv")
ageshelves<-read.csv(file="CHH shelf age associated sites.csv")
ageshores
ageshelves

Correlatedsites<-read.csv(file = "CHHCorrelatedsites.csv")
dim(Correlatedsites)
head(Correlatedsites)
head(ageislands)

ageislands<-ageislands[c(2,3,4)]
ageshores<-ageshores[c(2,3,4)]
ageshelves<-ageshelves[c(2,3,4)]
ageshores
ageislands
ageshelves
head(ageislands)

island<-ageislands[, 1]
island

head(Correlatedsites)
dim(Correlatedsites)


island
subset(Correlatedsites, (Site==246))
list<-setdiff(Correlatedsites$X, ageislands$X)   # elements in a$x NOT in b$y
list

length(list)


dim(Correlatedsites[Correlatedsites$X %in% list, ])

Correlatedminusisland<-Correlatedsites[Correlatedsites$X %in% list, ]
dim(Correlatedminusisland)

Correlatedminusisland

#removing shelves

list<-setdiff(Correlatedminusisland$X, ageshelves$X)   # elements in a$x NOT in b$y
list

length(list) #went from 287 to 264, so lost 23 sites, whih is correct

Correlatedminusislandandshelf<-Correlatedminusisland[Correlatedminusisland$X %in% list, ]
dim(Correlatedminusislandandshelf)

write.csv(Correlatedminusislandandshelf, file = "CHH age associated minus island and shelf sites.csv")

#removing shores
list<-setdiff(Correlatedminusisland$X, ageshores$X)   # elements in a$x NOT in b$y
list

length(list)# going from 408 to 317, means 91 sites coming out

dim(Correlatedminusisland[Correlatedminusisland$X %in% list, ])

Correlatedminusislandandshore<-Correlatedminusisland[Correlatedminusisland$X %in% list, ]
dim(Correlatedminusislandandshore)

#removing shelves
list<-setdiff(Correlatedminusislandandshore$X, ageshelves$X)   # elements in a$x NOT in b$y
list

length(list)# going from 317 to 301, means  16 sites coming out

dim(Correlatedminusislandandshore[Correlatedminusislandandshore$X %in% list, ])

Correlatedminusislandandshoreandshelf<-Correlatedminusislandandshore[Correlatedminusislandandshore$X %in% list, ]
dim(Correlatedminusislandandshoreandshelf)

#now that all other sites removed, look at pos and neg values for open sea regions!

length(Correlatedminusislandandshoreandshelf$corrsites[Correlatedminusislandandshoreandshelf$corrsites>0]) # 149
#149 open sea sites that are positive with age 
length(Correlatedminusislandandshoreandshelf$corrsites[Correlatedminusislandandshoreandshelf$corrsites<0]) # 152
#152 CpG sites falling in open seas that are negative with age

write.csv(Correlatedminusislandandshoreandshelf, file = "CHH age associated open sea sites.csv")

###Graphing results

#Now that I have the results:

#CpG islands: 651,745 islands that exist across the loblolly genome
#background: 5027 / 30350 = 16.56% of the CpG sites we have are falling into islands 
#age associated: 3/10 = 30%


#CpG Shores: 1,144,719 shores that exist across the loblolly genome (in the scaffolds
#background: 4427 / 30350 = 14.58% of the sites we have fell into shores
#age-associated: 2/10 = 20%

#CPG shelves: 700,544 shelves that exist across the loblolly genome
#background: 1816 / 30350 =  0.0598 or 5.98% 
#age-associated: 1 / 10 = 10% 


#Open sea (by process of elimination)
#Genome:?
#Background: (Islands + shores = shelves = 11,270. 30350 - 11270 = 19080.) So, 19080 / 30350 = 0.6287 or 62.87%
#so 62.87% fall into open seas
#Age-associated: 4 / 10 = 40%. So 40% of the age-associated sites fall in open seas. 

#graphing

#vs ggplot
Region <- c(rep("Islands" , 2) , rep("Shores" , 2) , rep("Shelves" , 2) , rep("Open Seas" , 2) )
Sites <- rep(c("Background" , "Age-Associated") , 4)
value <- c(8.08, 6.82, 13.43, 11.04, 6.65, 7.47, 71.85, 74.68)
data <- data.frame(Region,Sites,value)
data


#setting the order for how I want it to show up in the plot
data$Region <- factor(data$Region, levels = c("Islands", "Shores", "Shelves", "Open Seas"))
data$Sites <- factor(data$Sites, levels = c("Background", "Age-Associated"))
data$Sites
data$Region

# Grouped
annotation<-ggplot(data, aes(fill=Sites, y=value, x=Region)) + 
  geom_bar(position="dodge", size = 0.5, color = "black", stat="identity")+ scale_fill_manual(values = c("black","red"))

annotation

#cytosine<-ggplot(data, aes(fill=context, y=value, x=status)) + 
#  geom_bar(position="stack", size = 0.5, color = "black", stat="identity") + scale_fill_manual(values = c("red","green","Blue"))

annotation<-annotation + ylim(range(c(0,100)))

annotation


annotation<-annotation +labs(y= "Sequenced CHH Methylation Sites (%)", x = "Genomic Context")
annotation

annotation<-annotation + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))
annotation


annotation<-annotation + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
annotation


annotation<-annotation + scale_y_continuous(limits = c(0,100), expand = c(0, 0))
annotation
annotation<-annotation + theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=14))


annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
annotation<-annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
annotation
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))

tiff(file="CHH Genomic context.tiff", res = 800, width = 189, height = 100, units = "mm")
annotation
dev.off()


#Binomial tests for enrichment of age-associated sites in genomic ranges

#CHH meth

##island (age associated vs background)
binom.test(21, 309, p = 0.0808,
           alternative = c("two.sided"),
           conf.level = 0.95) 

#shore
binom.test(34, 309, p = 0.1343,
           alternative = c("two.sided"),
           conf.level = 0.95) 

#shelf
binom.test(23, 309, p = 0.0665,
           alternative = c("two.sided"),
           conf.level = 0.95) 


#open sea
binom.test(231, 309, p = 0.7185,
           alternative = c("two.sided"),
           conf.level = 0.95) 


#so, no enrichment or depletions in genomic regions vs background

