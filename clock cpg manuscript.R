setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/CPG/DestrandF/mat files/")

library(psych)
library(caret)
#check to see how to load in libraries on cluster
library(methylKit)
library(impute)

#CpG
#start with 10x, sites across 80 percent of samples covered
mymat1<-read.csv(file="matCpG10desfalse80.csv", header = T, row.names="ID")
head(mymat1)


#transforming data so you can call variables and it is in the right order
lob1<- as.data.frame(t(mymat1))
#head(lob1)
dim(lob1)
# 24 26521


library(psych)

#removing near zero variant sites
nzv <- nearZeroVar(lob1)
filteredlob1 <- lob1[, -nzv]
dim(filteredlob1)
#now 18845 sites (so removed 8784 sites, or ~29% of the sites)
nzv

write.csv(nzv, "InvariantCpGsites10x80.csv")

#running new correlation with age

corrtest1 <- corr.test(filteredlob1[-1], filteredlob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p
#View(pval)
lowestpval<-sort(pval, decreasing = F)
View(lowestpval)
#6 sites significantly correlated compared to 10 with 5x 80

decreasingcor<-sort(corrtest1$r, decreasing = F)
View(decreasingcor) 
# strongly neg correlated (between -0.7 and -0.8) (same)
# between -0.6 and -0.7 (same)
# between -0.5 and -0.6  So no change in numbers before vs after filtering, just affected P values

increasingcor<-sort(corrtest1$r, decreasing = T)
View(increasingcor)
# between 0.8 and 0.9 (same)
# between 0.7 and 0.8 (same)
# between 0.6 and 0.7 (same)
# between 0.5 and 0.6 (same)




#generating histogram of sites correlated with Age
tiff(file="CpG10x80 age corr wlines des false filtered.tiff", res = 800, width = 100, height = 120, units = "mm")
#ZZ<-hist(corrtest1$r, main = "CpG5x80percent filtered",xaxt = "n", xlab = "Cytosine Age correlation", col = c("blue"), ylim = range(0:15000))
ZZ<-hist(corrtest1$r, main = " ",xaxt = "n", xlab = "Cytosine Age correlation", col = c("blue"), ylim = range(0:15000))

abline(h=0)
#if you dont want lines, hashtag the next line
abline(v=c(-0.5, 0.5), lty = 2)
ticks<-c(-0.5, 0.0, 0.5)
axis(side = 1, at = ticks, labels = c("-0.5", "0", "0.5"), line = -0.6)
dev.off()


#now imputing


library(impute)

bclockt <- as.matrix(t(filteredlob1))

#View(bclockt)

dim(bclockt)

#rowmax will be the MAX amount of missing data allowed in any row (the default is 50%).
#Note that any rows that have more than the value I set as missing will be imputed using the overall mean per sample

#col max is the max percent of missing data allowed in the columns

#rows are sites, columns are individuals (1-24)

#maxp is the largest block of genes that are imputed using the knn algorithm
imputelob1 <- impute.knn(bclockt,k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069) ##choosing  k=2 because smallest age group has 3 samples in it

imputedlob1 <- as.matrix(imputelob1$data)

#View(imputedlob1)

##Works!!
is.na(imputedlob1)
is.matrix(imputedlob1)

Lob1<-as.data.frame(t(imputedlob1))
dim(Lob1)


corrtest1 <- corr.test(Lob1[-1], Lob1$AGE, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj 
#View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#9 sites with pearsons correlations that are significant

#elastic net 



#running with 5 individuals in test set, 19 in train


dim(Lob1)

lob_age <- Lob1[c(1)]
#View(lob_age) #pulled age for each sample
dim(lob_age)
#lob_age[c(4, 5, 6), c(1)]


dim(Lob1)


lob_CpG <- Lob1[c(2:18845)]#just pulled all the CpG sites (now columns), for every individual (1-28), and it gives the 
#amount of methylation present at each site (need to see what values mean) (go back to what generated "meth")

#View(lob_CpG)
dim(lob_CpG)
#so, it looks like its going to run model on how values at all 146918 sites vary with age prediction

##subset x to make training and test matrices (randomly chose first 3 for test)
#(need to check whats going on here, why)
#so what this did was pull out 3 of the 24 individual trees, (c(4,5,6) which was individuals 4,5 , and 6, 
# each with chronological ages of 37 years (lob_age[c(4, 5, 6), c(1)]), and kept all of their methylations at all 1991 sites)
test.CpG <- lob_CpG[c(14, 16, 20, 22, 24), c(1:18844)]#2615 if 5x file
#View(test.CpG)
dim(test.CpG)
#so the test matrix was 3 trees, all of the same age? May play with this to see what happens when changing

#will come back and redo with trees of differing ages


#this is taking the remaining individuals of varying ages 

train.CpG <- lob_CpG[c(1:13, 15, 17:19, 21, 23), c(1:18844)]#2615 if 5 x
#View(train.CpG)
dim(train.CpG)
#

##subset ages
##subset x to make training and test matrices (randomly chose first 3 for test)
l_age <- as.matrix(lob_age) ##need to make it a matrix first so that the ID names come along with the ages
#head(l_age)
test.age <- l_age[c(14,16,20,22,24), c(1)]
#View(test.age)
dim(test.age)
#so, test age and train age are being put in as matrices
train.age <- l_age[c(1:13, 15, 17:19, 21, 23), c(1)]
#View(train.age)
length(train.age)


##making data readable by glmnet (it doesn't like the character variables in a matrix)
CpG.train <- as.matrix(as.data.frame(lapply(train.CpG, as.numeric)))
CpG.test <- as.matrix(as.data.frame(lapply(test.CpG, as.numeric)))


#install.packages("glmnet")
library(glmnet)
set.seed(123) #setting seed to get similar results

#using 10 fold cross validation to estimate the lambda parameter 
#note that the nfolds argument needs to be set to the highest number of samples I am running this on
#"leave one out approach"
lob.cv.day <- cv.glmnet(as.matrix(CpG.train), 
                        as.matrix(train.age),
                        alpha=0.5, nfolds=19, family="gaussian")
#plot(lob.cv.day)


#define the lambda parameter
best_lambda <- lob.cv.day$lambda.min
best_lambda #best lambda of CpG imputed data set is 0.8965374 when using n folds = 21 but using the differing ages


#I need to see what this means and where and how it comes into play

#in downstream analysis, this shows up at the lambda value for the last coefficient that is still explaining any of
#the variation. 
#fit the elastic net model to the data

lob.model <- glmnet(as.matrix(CpG.train), 
                    as.matrix(train.age),
                    alpha=0.5, family="gaussian", nlambda=100)
#plot the coordinate descent path
#plot(lob.model, label = T)

#print glmnet path
#print(lob.model)
#so in my model, df is the number of non zero coefficients, so does that mean that if I have 88 coefficients,
#then this model used 88 of the 1991 CpG sites found in the Bam files to produce this age model?

#If so, then I need a way to go back and then find which sites these corresponded to, then annotate them

#plot(lob.model, xvar = "lambda", label = TRUE)
lob.model$dev.ratio



#coefficients for fitted model
age_coeff.lob <- coef(lob.model, s=best_lambda)
dim(age_coeff.lob)
#dimensions are 21567

#tells me which sites have non zero coefficient  values, so I can find out which sites were used
which(age_coeff.lob!=0)   
age_coeff.lob
#36 sites

list<-which(age_coeff.lob!=0) 
list

age_coeff.lob[c(list),]



#fitted values
fitted.lob <- predict(lob.model, CpG.train, s=best_lambda)

#plot actual by predicted values for training data
actual.fitted <- data.frame(fitted.lob, train.age)
actual.fitted
plot(actual.fitted$s1~actual.fitted$train.age, xlab="chronological age", ylab="predicted age")
abline(lm(actual.fitted$s1~actual.fitted$train.age), col="red") # regression line (y~x)

results<-lm(actual.fitted$s1~actual.fitted$train.age)
summary(results)


top5trainage<-ggplot(actual.fitted, aes(x = train.age, y = s1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

top5trainage



top5trainage<-top5trainage +labs(y= "Predicted Age", x = "Chronological Age")
top5trainage

top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
top5trainage


top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage


top5trainage<-top5trainage + scale_y_continuous(limits = c(0,120), expand = c(0, 0))
top5trainage
top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

top5trainage<-top5trainage + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=14))
top5trainage<-top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage+ theme(plot.margin = unit(c(0.25,1,0.25,0.25), "cm"))
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="testmodelcpgtop5.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()







#calculate R2 for fitted data
R2.train <- cor(fitted.lob, train.age)^2
R2.train
#mine is 0.9998667 when nfolds is 21
#so better when using these samples 

actual.fitted$s1


##retrieve predicted ages
actual.fitted$X1

#residuals(actual_predictedage)

##get the error associated with each predicted ages
train_predicted_error <- (actual.fitted$s1 - actual.fitted$train.age)
View(train_predicted_error)
train_predicted_error

trainerror<-c(0.4846985, 0.9038177, 0.3103076, 0.3099947, 0.1504449, 0.4655714,  0.1925104,  0.9315863,  0.1819918,
              0.7685617,  0.7091500,  0.7951598, 0.8707412, 0.7466938,  0.8873453,  0.9614224,  0.4390744, 0.1686975,
              1.4558348)

df<-as.data.frame(trainerror)
View(df)

View(trainerror)
mean(trainerror)
#mean is 0.6175581 years
Se<-sd(trainerror)/sqrt(19)
Se#se is 0.08, so error is 0.61755(+/- 0.08)



#predicted values for test set
predicted <- predict(lob.model, CpG.test, s=best_lambda)
predicted
test.age

#plot actual by predicted values for training data
actual.predicted <- data.frame(predicted, test.age)
actual.predicted
plot(actual.predicted$s1~actual.predicted$test.age, xlab="chronological age", ylab="predicted age")
##fit line
abline(lm(actual.predicted$s1~actual.predicted$test.age), col="red") # regression line (y~x)
#note that this model is predicting that the 3 test trees that were all chronologically 37 years old
#however this is predicting that they are 72, 48, 38 years old

results<-lm(actual.predicted$s1~actual.predicted$test.age)
summary(results)

top5trainage<-ggplot(actual.predicted, aes(x = test.age, y = s1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

top5trainage


top5trainage<-top5trainage +labs(y= "Predicted Age", x = "Chronological Age")
top5trainage

top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
top5trainage


top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage

top5trainage<-top5trainage + scale_y_continuous(expand = c(0, 0))
top5trainage

#top5trainage<-top5trainage + scale_y_continuous(limits = c(0,120), expand = c(0, 0))
top5trainage
#top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

top5trainage<-top5trainage + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=14))
#top5trainage<-top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage+ theme(plot.margin = unit(c(0.25,1,0.25,0.25), "cm"))
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="testmodelcpgtop5.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()







#calculate R2 for fitted data
R2.test <- cor(predicted, test.age)^2
R2.test
#R2 was 0.7905 (so not great) 



##get the error associated with each predicted ages
test_predicted_error <- (actual.predicted$s1 - actual.predicted$test.age)
View(test_predicted_error)

test_predicted_error

testerror<-c(12.782652,  1.551413,  3.346544, 29.131230, 66.587677)
df<-as.data.frame(testerror)
View(df)
mean(testerror)
#mean is 22.679
Se<-sd(testerror)/sqrt(5)
Se
#Se is 12.016, so error is 22.679 (+/- 12.016)



#Doing pearson correlations and linear models

#need to make a train set to obtain the correlated sites, then run on the test individuals
train<-c(1:13,15,17:19,21,23)
#test<-c(14,16,20,22,24)
dim(Lob1[c(train),])

Lobtrain<-Lob1[c(train),]

corrtest1 <- corr.test(Lobtrain[-1], Lobtrain$AGE, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#3 sites significant
r<-subset(corrtest1$r, corrtest1$p<0.05)
p<-subset(pval,pval<0.05)
sigsites<-data.frame(r,p)
sigsites

a<-subset(corrtest1$r, corrtest1$r > 0.7)
b<-subset(corrtest1$r, corrtest1$r < -0.7)
c<-subset(corrtest1$p, corrtest1$r > 0.7)
d<-subset(corrtest1$p, corrtest1$r < -0.7)


a
b
d

corrsites<-rbind(a,b)
corrsites

corrpval<-rbind(c,d)


correlated<-data.frame(corrsites, corrpval)
correlated
dim(correlated)
correlated

highestcorr<-sort(correlated$corrsites, decreasing = T)
View(highestcorr)

write.csv(correlated, file = "traincorrelatedcpg.csv")

#using 5 sites to estimate age:

#testing if doing pearson on only train model is same as test model:



Lob2<-Lob1[-1]
#Lob2[c(list)]
a<-Lobtrain$'21133'
b<-Lobtrain$'2070'
c<-Lobtrain$'15746'
d<-Lobtrain$'12447'
#e<-Lob2$'15746'
#f<-Lob2$'15748'
g<-Lobtrain$'13865'
#h<-Lob2$'15751'
#i<-Lob2$'15753'
dim(Lob2)
pearsonmodel<-data.frame(a,b,c,d,g)
pearsonmodel


age<-Lobtrain[c(1)]
age


pearsonmodel<-data.frame(pearsonmodel, age)
pearsonmodel#need to divide this into train and test values and use model from train to predict test
#divide the pearsonmodel into train (1-13, 15, 17-19, 21,23)
#divide into test (14,16,20,22,24)
#train<-c(1:13,15,17:19,21,23)
#test<-c(14,16,20,22,24)
#train
#test
#trainmodel<-pearsonmodel[c(train),]
#dim(trainmodel)
#testmodel<-pearsonmodel[c(test),]
#dim(testmodel)
#testmodel<-testmodel[-c(10)]
#trainmodel
#testmodel
results<-lm(AGE ~ a+b+c+d+g, data = pearsonmodel)
summary(results)

predict(results, newdata = pearsonmodel)
pearsonmodel$AGE

chrono<-c(46, 46, 46, 37, 37, 37, 19, 19, 19, 10, 10, 10, 55, 55,  1,  1, 28, 28, 97)
pred<-c(45.79964, 31.74183, 50.00420, 39.31835, 36.31914, 22.83154, 13.17457, 29.04555, 19.79061, 16.70574, 12.00511, 11.69973, 57.06210, 55.98719,  9.49734,  8.67715, 33.23731, 12.86361, 95.23929)

results<-lm(pred~chrono)
summary(results)





plot(AGE ~ a+b+c+d+g, data = trainmodel)
#a,b,c,d,g are sites 2070, 12447, 13240, 15143, 15749

trainmodel$a
is.data.frame(trainmodel)
#vs ggplot
#Error <- c(trainmodel$a, trainmodel$b, trainmodel$c, trainmodel$d, trainmodel$g)
#Error
#Context <- c(rep("2070", 19), rep("12447", 19), rep("13240", 19), rep("15143", 19), rep("15749", 19))
#Age <- c()
#data <- data.frame(Error, Context, value)
#data
trainmodel

#setting the order for how I want it to show up in the plot
data$Context <- factor(data$Context, levels = c("CpG", "CHG", "CpG and CHG"))
#data$Sites <- factor(data$Sites, levels = c("CpG", "CHG", "CpG and CHG"))
#data$Sites

pearsonmodel$AGE
trainmodel$AGE
Methylation <- c(pearsonmodel$a, pearsonmodel$b, pearsonmodel$c, pearsonmodel$d, pearsonmodel$g)
Error
Age <- c(pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE)
Site <- c(rep("21133", 19), rep("2070", 19), rep("15746", 19), rep("12447", 19), rep("13865", 19))
data<-data.frame(Age, Methylation, Site)
#Age <- c()
#data <- data.frame(Error, Context, value)
data
#length(c(rep(trainmodel$AGE), 5))
#length(c(trainmodel$AGE, trainmodel$AGE))

trainplot<-ggplot(data, aes(x = Methylation, y = Age)) 
trainplot<-trainplot+ geom_point(aes(colour = Site), size = 1.5) 
trainplot

trainplot<-trainplot +labs(y= "Age", x = "Methylation %")
trainplot

trainplot<-trainplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))
trainplot


trainplot<-trainplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
trainplot

trainmodel$AGE
trainplot<-trainplot + scale_y_continuous(limits = c(0,100), expand = c(0, 0))
trainplot
trainplot<-trainplot + theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=14))


#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
trainplot<-trainplot + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
trainplot
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


trainplot + geom_abline(intercept = 3.8663)

coef(lm(AGE ~ a + b + c + d + g, data = trainmodel))
trainplot<-trainplot + geom_smooth(method = "lm", color = "red")
trainplot

#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="ActualTrainmodel 5 corr sites.tiff", res = 800, width = 189, height = 100, units = "mm")
trainplot
dev.off()



#predicting test model

test<-c(14,16,20,22,24)
dim(Lob1[c(test),])

Lobtest<-Lob1[c(test),]


a<-Lobtest$'21133'
b<-Lobtest$'2070'
c<-Lobtest$'15746'
d<-Lobtest$'12447'
#e<-Lob2$'15746'
#f<-Lob2$'15748'
g<-Lobtest$'13865'
#h<-Lob2$'15751'
#i<-Lob2$'15753'
dim(Lob2)
pearsontestmodel<-data.frame(a,b,c,d,g)
pearsontestmodel


age<-Lobtest[c(1)]
age


pearsontestmodel<-data.frame(pearsontestmodel, age)
pearsontestmodel#need to divide this into train and test values and use model from train to predict test

results<-lm(AGE ~ a+b+c+d+g, data = pearsonmodel)
summary(results)


#model R2 was 0.8831 (multiple), adjusted was 0.8381


#getting r2 for train set 
predict(results, newdata = pearsonmodel)
pearsonmodel$AGE

chrono<-c(46, 46, 46, 37, 37, 37, 19, 19, 19, 10, 10, 10, 55, 55,  1,  1, 28, 28, 97)
pred<-c(45.79964, 31.74183, 50.00420, 39.31835, 36.31914, 22.83154, 13.17457, 29.04555, 19.79061, 16.70574, 12.00511, 11.69973, 57.06210, 55.98719,  9.49734,  8.67715, 33.23731, 12.86361, 95.23929)

results<-lm(pred~chrono)
summary(results)


#MAe for train set 
Error<-chrono - pred
Error
Error<-c(0.20036,  14.25817,  4.00420,  2.31835,   0.68086,  14.16846,   5.82543, 10.04555,  0.79061,  6.70574,
         2.00511,  1.69973,  2.06210,  0.98719,  8.49734,  7.67715,  5.23731,  15.13639,   1.76071)
df<-as.data.frame(Error)
View(df)
mean(Error)#5.48 years
SE<-sd(Error)/sqrt(19)
SE#1.13, so MAE for train set is 5.47 (+/- 1.13) years
#so R2 for prediction of age vs chrono age in train model was 0.8976

pearsonmodel<-data.frame(pearsonmodel,pred)
pearsonmodel<-data.frame(pearsonmodel,pred,Error)
pearsonmodel
#results<-lm(Error~AGE, data = Pearsonmodel)
#summary(results)

top5trainage<-ggplot(pearsonmodel, aes(x = AGE, y = pred)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

top5trainage



top5trainage<-top5trainage +labs(y= "Predicted Age", x = "Chronological Age")
top5trainage

top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
top5trainage


top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage


top5trainage<-top5trainage + scale_y_continuous(limits = c(0,120), expand = c(0, 0))
top5trainage
top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

top5trainage<-top5trainage + theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14))
top5trainage<-top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage+ theme(plot.margin = unit(c(0.25,1,0.25,0.25), "cm"))
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="trainmodelcpgtop5.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()









#getting r2 for test set 

predict(results, newdata = pearsontestmodel)

chrono<-c(55, 1, 28, 82, 119)
pred<-c(29.343,   8.827,  19.488,  65.344, 57.534)

error<-chrono - pred
error

error<-c(25.657, 7.827,  8.512, 16.656, 61.466)
View(as.data.frame(error))

mean(error)
#so mean error was 24.0236
SE<-sd(error)/sqrt(5)
SE#st err was 9.904. So MAE is 24.02 (+/- 9.904) years. elastic net was 22.679 (+/- 12.016)


pearsontestmodel<-data.frame(pearsontestmodel,pred, error)
pearsontestmodel
results<-lm(pred~AGE, data = pearsontestmodel)
summary(results)

top5trainage<-ggplot(pearsontestmodel, aes(x = AGE, y = pred)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

top5trainage



top5trainage<-top5trainage +labs(y= "Predicted Age", x = "Chronological Age")
top5trainage

top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
top5trainage


top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage


top5trainage<-top5trainage + scale_y_continuous(limits = c(0,120), expand = c(0, 0))
top5trainage
top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

top5trainage<-top5trainage + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=14))
top5trainage<-top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage+ theme(plot.margin = unit(c(0.25,1,0.25,0.25), "cm"))
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="testmodelcpgtop5.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



#doing for top 10

#Doing pearson correlations and linear models

#need to make a train set to obtain the correlated sites, then run on the test individuals
train<-c(1:13,15,17:19,21,23)
#test<-c(14,16,20,22,24)
dim(Lob1[c(train),])

Lobtrain<-Lob1[c(train),]

corrtest1 <- corr.test(Lobtrain[-1], Lobtrain$AGE, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p
View(pval)
lowestpval<-sort(pval, decreasing = F)
View(lowestpval)
#3 sites significant
r<-subset(corrtest1$r, corrtest1$p<0.05)
p<-subset(pval,pval<0.05)
sigsites<-data.frame(r,p)
sigsites

a<-subset(corrtest1$r, corrtest1$r > 0.7)
b<-subset(corrtest1$r, corrtest1$r < -0.7)
c<-subset(corrtest1$p, corrtest1$r > 0.7)
d<-subset(corrtest1$p, corrtest1$r < -0.7)

corrsites<-rbind(a,b)
corrsites

corrpval<-rbind(c,d)


correlated<-data.frame(corrsites, corrpval)
correlated
dim(correlated)
correlated

highestcorr<-sort(correlated$corrsites, decreasing = T)
View(highestcorr)

write.csv(correlated, file = "traincorrelatedcpg.csv")

#using 5 sites to estimate age:

#testing if doing pearson on only train model is same as test model:



Lob2<-Lob1[-1]
#Lob2[c(list)]
a<-Lobtrain$'21133'
b<-Lobtrain$'2070'
c<-Lobtrain$'15746'
d<-Lobtrain$'12447'
e<-Lobtrain$'13865'
f<-Lobtrain$'4695'
g<-Lobtrain$'953'
h<-Lobtrain$'9795'
i<-Lobtrain$'13140'
j<-Lobtrain$'8564'
dim(Lob2)
pearsonmodel<-data.frame(a,b,c,d,e,f,g,h,i,j)
pearsonmodel


age<-Lobtrain[c(1)]
age


pearsonmodel<-data.frame(pearsonmodel, age)
pearsonmodel#need to divide this into train and test values and use model from train to predict test

results<-lm(AGE ~ a+b+c+d+e+f+g+h+i+j, data = pearsonmodel)
summary(results)


trainmodel$a
is.data.frame(trainmodel)

trainmodel

#setting the order for how I want it to show up in the plot
data$Context <- factor(data$Context, levels = c("CpG", "CHG", "CpG and CHG"))
#data$Sites <- factor(data$Sites, levels = c("CpG", "CHG", "CpG and CHG"))
#data$Sites

pearsonmodel$AGE
trainmodel$AGE
Methylation <- c(pearsonmodel$a, pearsonmodel$b, pearsonmodel$c, pearsonmodel$d, pearsonmodel$e, pearsonmodel$f, pearsonmodel$g, pearsonmodel$h, pearsonmodel$i, pearsonmodel$j)
#Error
Age <- c(pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE,pearsonmodel$AGE)
Site <- c(rep("21133", 19), rep("2070", 19), rep("15746", 19), rep("12447", 19), rep("13865", 19), rep("4695", 19), rep("953", 19), rep("9795", 19), rep("13140", 19), rep("8564", 19))
data<-data.frame(Age, Methylation, Site)

data


trainplot<-ggplot(data, aes(x = Methylation, y = Age)) 
trainplot<-trainplot+ geom_point(aes(colour = Site), size = 1.5) 
trainplot

trainplot<-trainplot +labs(y= "Age", x = "Methylation %")
trainplot

trainplot<-trainplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
trainplot


trainplot<-trainplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
trainplot

trainmodel$AGE
trainplot<-trainplot + scale_y_continuous(limits = c(0,100), expand = c(0, 0))
trainplot
trainplot<-trainplot + theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=14))


#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
trainplot<-trainplot + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
trainplot
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


trainplot + geom_abline(intercept = 3.8663)

coef(lm(AGE ~ a + b + c + d + g, data = trainmodel))
trainplot<-trainplot + geom_smooth(method = "lm", color = "red")
trainplot

#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="ActualTrainmodel 10 corr sites.tiff", res = 800, width = 189, height = 100, units = "mm")
trainplot
dev.off()



#predicting test model

test<-c(14,16,20,22,24)
dim(Lob1[c(test),])

Lobtest<-Lob1[c(test),]


a<-Lobtest$'21133'
b<-Lobtest$'2070'
c<-Lobtest$'15746'
d<-Lobtest$'12447'
e<-Lobtest$'13865'
f<-Lobtest$'4695'
g<-Lobtest$'953'
h<-Lobtest$'9795'
i<-Lobtest$'13140'
j<-Lobtest$'8564'
dim(Lob2)
pearsontestmodel<-data.frame(a,b,c,d,e,f,g,h,i,j)
pearsontestmodel
dim(Lob2)



age<-Lobtest[c(1)]
age


pearsontestmodel<-data.frame(pearsontestmodel, age)
pearsontestmodel#need to divide this into train and test values and use model from train to predict test

results<-lm(AGE ~ a+b+c+d+e+f+g+h+i+j, data = pearsonmodel)
summary(results)


#model R2 was 0.8831 (multiple), adjusted was 0.8381


#getting r2 for train set 
predict(results, newdata = pearsonmodel)
pearsonmodel$AGE

chrono<-c(46, 46, 46, 37, 37, 37, 19, 19, 19, 10, 10, 10, 55, 55,  1,  1, 28, 28, 97)
pred<-c(43.715412, 47.448405, 43.435286, 34.713916, 34.650175, 23.470360, 14.923763, 29.730759, 17.587069, 16.555494, 10.708692, 14.492738, 53.965193, 56.342523,  3.601402,  3.754067, 35.130151, 18.869763, 97.904830)

results<-lm(pred~chrono)
summary(results)


#MAe for train set 
Error<-chrono - pred
Error
Error<-c(2.284588,  1.448405,   2.564714,   2.286084,   2.349825,  13.529640,   4.076237, 10.730759,   1.412931,
         6.555494,  0.708692,  4.492738,   1.034807,  1.342523,  2.601402,  2.754067,  7.130151,   9.130237,
         0.904830)
View(as.data.frame(Error))
mean(Error)#4.07 myears
SE<-sd(Error)/sqrt(19)
SE#.842, so MAE for train set is 4.14 (+/- 1.02) years
#so R2 for prediction of age vs chrono age in train model was 0.9308


pearsonmodel<-data.frame(pearsonmodel,pred,Error)
pearsonmodel
#results<-lm(Error~AGE, data = Pearsonmodel)
#summary(results)

top5trainage<-ggplot(pearsonmodel, aes(x = AGE, y = pred)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

top5trainage



top5trainage<-top5trainage +labs(y= "Predicted Age", x = "Chronological Age")
top5trainage

top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
top5trainage


top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage


top5trainage<-top5trainage + scale_y_continuous(limits = c(0,120), expand = c(0, 0))
top5trainage
top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

top5trainage<-top5trainage + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=14))
top5trainage<-top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage+ theme(plot.margin = unit(c(0.25,1,0.25,0.25), "cm"))
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="trainmodelcpgtop10.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



#getting r2 for test set 

predict(results, newdata = pearsontestmodel)

chrono<-c(55, 1, 28, 82, 119)
pred<-c(31.026358, 28.321702, 40.326131, 68.311643,  5.617678 )

error<-chrono - pred
error

error<-c(23.97364, 27.32170, 12.32613,  13.68836, 113.38232)
View(as.data.frame(error))
mean(error)
#so mean error was 38.13843
SE<-sd(error)/sqrt(5)
SE#st err was 19.03. So MAE is 38.14 (+/- 19.03) years. elastic net was 22.679 (+/- 12.016)


pearsontestmodel<-data.frame(pearsontestmodel,pred, error)
pearsontestmodel
results<-lm(error~AGE, data = pearsontestmodel)
summary(results)

top5trainage<-ggplot(pearsontestmodel, aes(x = AGE, y = pred)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

top5trainage



top5trainage<-top5trainage +labs(y= "Predicted Age", x = "Chronological Age")
top5trainage

top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
top5trainage


top5trainage<-top5trainage + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.title.x = element_text(vjust=-1), axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage


top5trainage<-top5trainage + scale_y_continuous(limits = c(0,120), expand = c(0, 0))
top5trainage
top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

top5trainage<-top5trainage + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=14))
top5trainage<-top5trainage + scale_x_continuous(limits = c(0,120), expand = c(0, 0))

#annotation + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage + theme(axis.text.x = element_text(size = 12, colour = "black"), axis.text.y = element_text(size = 12, colour = "black"))
#annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
#annotation<-annotation + theme(axis.title.y = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 12, colour = "black"))
top5trainage<-top5trainage+ theme(plot.margin = unit(c(0.25,1,0.25,0.25), "cm"))
#annotation<-annotation + theme(axis.text=element_text(size=12),
#                              axis.title=element_text(size=14))


#abline(lm(actual.fitted$X1~actual.fitted$train.age), col="red")
tiff(file="testmodelcpgtop10.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()





#genomic ranges 



#CpG elastic net sites
library(GenomicRanges)

#Getting CpG ELastic net sites

# read CpGi data set
#filePath=system.file("Loblolly_cpg_islands.bed")
#filePath
setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
colnames(Pearsonmodel)

#putting in sites for CpG elastic net 
Agesites<-c(953,
            1452,
            1467,
            1834,
            2070,
            3831,
            4315,
            4695,
            5028,
            5733,
            6692,
            7141,
            7690,
            8564,
            10119,
            10155,
            11268,
            12446,
            12447,
            13865,
            13976,
            14545,
            15273,
            15746,
            15749,
            15751,
            15753,
            19760,
            20247,
            20699,
            21133,
            22152,
            23443,
            24505,
            24970)
Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#2 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Islandsites<-newagesites[c(overlapping),]
Islandsites
#site 12447 falls into a shore!!
is.data.frame(Shelfsites)

write.csv(Islandsites, file = "islands.csv")




setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
#Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
#colnames(Pearsonmodel)

#putting in sites for CpG elastic net 

Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#6 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)
newagesites
Shoresites<-newagesites[c(overlapping),]
Shoresites
#site 12447 falls into a shore!!
#is.data.frame(Shelfsites)

write.csv(Shoresites, file = "CpGelasticnetshores.csv")





setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
#Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
#colnames(Pearsonmodel)

#putting in sites for CpG elastic net 


Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#0 elastic net sites fell into shelves!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Shelfsites<-newagesites[c(overlapping),]
Shelfsites
#site 12447 falls into a shore!!
#is.data.frame(Shelfsites)

write.csv(Shelfsites, file = "CpGelasticnetshores.csv")

Islandsites
Shoresites
Shelfsites



#Additional sites for combined CpG and CHG e net model

setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)



#putting in sites for CpG elastic net 
Agesites<-c(9487, 10374, 14243, 23787)
Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#2 elastic net both sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Islandsites<-newagesites[c(overlapping),]
Islandsites
#site 12447 falls into a shore!!
is.data.frame(Shelfsites)

write.csv(Islandsites, file = "islands.csv")




setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
#Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
#colnames(Pearsonmodel)

#putting in sites for CpG elastic net 
#Agesites<-c(2070,
#            12447,
#            13865,
#            15746,
#            21133)
Agesites


#filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
#filepathpeaks[c(Agesites),]
#newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#0 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Shoresites<-newagesites[c(overlapping),]
Shoresites
#site 12447 falls into a shore!!
#is.data.frame(Shelfsites)

write.csv(Shoresites, file = "CpGelasticnetshores.csv")





setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
#Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
#colnames(Pearsonmodel)

#putting in sites for CpG elastic net 
#Agesites<-c(2070,
#            12447,
#            13865,
#            15746,
#            21133)
Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#0 elastic net sites fell into shelves!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Shelfsites<-newagesites[c(overlapping),]
Shelfsites
#site 12447 falls into a shore!!
#is.data.frame(Shelfsites)

write.csv(Shelfsites, file = "CpGelasticnetshores.csv")

Islandsites
Shoresites
Shelfsites



#CpG top 5 model sites
#Getting CpG Pearson sites

# read CpGi data set
#filePath=system.file("Loblolly_cpg_islands.bed")
#filePath
setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)



#putting in sites for CpG elastic net 
Agesites<-c(2070,
            12447,
            13865,
            15746,
            21133)
Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#0 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Islandsites<-newagesites[c(overlapping),]
Islandsites
#site 12447 falls into a shore!!
is.data.frame(Shelfsites)

write.csv(Islandsites, file = "islands.csv")




setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
#Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
#colnames(Pearsonmodel)

#putting in sites for CpG elastic net 
Agesites<-c(2070,
            12447,
            13865,
            15746,
            21133)
Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#1 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Shoresites<-newagesites[c(overlapping),]
Shoresites
#site 12447 falls into a shore!!
#is.data.frame(Shelfsites)

write.csv(Shoresites, file = "CpGelasticnetshores.csv")





setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks

dim(filepathpeaks)


#setwd("~/Parrott lab/Research/loblolly pine experiment/ALL")
#Pearsonmodel
#agesites<-read.csv(file="CpGCorrelatedsites.csv")
#agesites
#colnames(Pearsonmodel)

#putting in sites for CpG elastic net 
Agesites<-c(2070,
            12447,
            13865,
            15746,
            21133)
Agesites


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#0 elastic net sites fell into shelves!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Shelfsites<-newagesites[c(overlapping),]
Shelfsites
#site 12447 falls into a shore!!
#is.data.frame(Shelfsites)

write.csv(Shelfsites, file = "CpGelasticnetshores.csv")

Islandsites
Shoresites
Shelfsites


#Additional 2 sites needed for CpG top 10 model 

#final 2 sites for c[g top 10 model 
setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_islands.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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
CpGmeth<-read.csv(file = "methCpG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks


Agesites<-c(9795, 13140)


filepathpeaks<-filepathpeaks[c(-1)]
#use Agesites to pull out sites from filepathpeaks 
filepathpeaks[c(Agesites),]
newagesites<-filepathpeaks[c(Agesites),]


newagesites


newagesites[,1]

sitenames<-newagesites[,1]
sitenames

startsites<-as.numeric(newagesites[,2])
endsites<-as.numeric(newagesites[,3])
endsites
cpgi.df

#making sites to check overaps
cpgisites.gr=GRanges(seqnames = sitenames,
                     ranges=IRanges(start=startsites,
                                    end=endsites))
cpgi.gr



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#2 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Islandsites<-newagesites[c(overlapping),]
Islandsites

setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Loblolly_cpg_shores.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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



# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#2 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Islandsites<-newagesites[c(overlapping),]
Islandsites




setwd("~/Parrott lab/Research/loblolly pine experiment/Annotation file")

cpgi.df = read.table("Lobolly_cpg_shelves.bed", header = FALSE,stringsAsFactors=FALSE)
cpgi.df

head(filePath)

# remove chr names:
#cpgi.df[grep("_",cpgi.df[,1],invert=TRUE),]


#cpgi.df =filePath [grep("_",cpgi.df[,1],invert=TRUE),]

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

# get the peaks that overlap with CpG islands
subsetByOverlaps(cpgisites.gr,cpgi.gr)


#2 elastic net sites fell into islands!!

#the below will tell you if sites are falling in ranges. Many will have zeros
counts2=countOverlaps(cpgisites.gr,cpgi.gr)
head(counts2)

#the find)verlaps function tells me which sites (the queryHits integers) that are falling in the islands,
# and the subjectHits are telling me which island they are falling in.


#this is the one we are interested in
overlaps<-findOverlaps(cpgisites.gr,cpgi.gr)
overlapping<-overlaps@from

dim(newagesites)

Islandsites<-newagesites[c(overlapping),]
Islandsites


