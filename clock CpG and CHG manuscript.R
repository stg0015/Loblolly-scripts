setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/CpG and CHG clock/")

library(psych)
library(caret)
#check to see how to load in libraries on cluster
library(methylKit)
library(impute)

#CpG
#start with 10x, sites across 80 percent of samples covered
mymat1<-read.csv(file="matCpGandCHG10desfalse80.csv", header = T, row.names="ID")
head(mymat1)


#transforming data so you can call variables and it is in the right order
lob1<- as.data.frame(t(mymat1))
#head(lob1)
dim(lob1)
# 24 57249


library(psych)

#removing near zero variant sites
nzv <- nearZeroVar(lob1)
filteredlob1 <- lob1[, -nzv]
dim(filteredlob1)
#now 41108 sites (so removed 8784 sites, or ~29% of the sites)
nzv

write.csv(nzv, "InvariantCpGsites10x80.csv")

#running new correlation with age

corrtest1 <- corr.test(filteredlob1[-1], filteredlob1$AGE, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
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


imputelob1 <- impute.knn(bclockt,k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069) ##choosing  k=2 because smallest age group has 3 samples in it

imputedlob1 <- as.matrix(imputelob1$data)


is.na(imputedlob1)
is.matrix(imputedlob1)

Lob1<-as.data.frame(t(imputedlob1))
dim(Lob1)


corrtest1 <- corr.test(Lob1[-1], Lob1$AGE, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
#View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#25 sites with pearsons correlations that are significant

#elastic net 



#running with 5 individuals in test set, 19 in train


dim(Lob1)

lob_age <- Lob1[c(1)]
#View(lob_age) #pulled age for each sample
dim(lob_age)
#lob_age[c(4, 5, 6), c(1)]


dim(Lob1)


lob_CpG <- Lob1[c(2:41108)]#just pulled all the CpG sites (now columns), for every individual (1-28), and it gives the 
#amount of methylation present at each site (need to see what values mean) (go back to what generated "meth")

#View(lob_CpG)
dim(lob_CpG)
#so, it looks like its going to run model on how values at all 146918 sites vary with age prediction


test.CpG <- lob_CpG[c(14, 16, 20, 22, 24), c(1:41107)]#2615 if 5x file

dim(test.CpG)


train.CpG <- lob_CpG[c(1:13, 15, 17:19, 21, 23), c(1:41107)]#2615 if 5 x
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

#View(CpG.train)
#View(CpG.test)

#install.packages("glmnet")
library(glmnet)
set.seed(123) #setting seed to get similar results

#using 10 fold cross validation to estimate the lambda parameter 

lob.cv.day <- cv.glmnet(as.matrix(CpG.train), 
                        as.matrix(train.age),
                        alpha=0.5, nfolds=19, family="gaussian")



#define the lambda parameter
best_lambda <- lob.cv.day$lambda.min
best_lambda #best lambda of CpG imputed data set is 0.9839 when using n folds = 21 but using the differing ages


lob.model <- glmnet(as.matrix(CpG.train), 
                    as.matrix(train.age),
                    alpha=0.5, family="gaussian", nlambda=100)
#plot the coordinate descent path
#plot(lob.model, label = T)

#print glmnet path
#print(lob.model)



#coefficients for fitted model
age_coeff.lob <- coef(lob.model, s=best_lambda)
dim(age_coeff.lob)
#dimensions are 21567

#tells me which sites have non zero coefficient  values, so I can find out which sites were used
which(age_coeff.lob!=0)   
age_coeff.lob
#48 sites

list<-which(age_coeff.lob!=0) 
list

age_coeff.lob[c(list),][c(1)]



#fitted values
fitted.lob <- predict(lob.model, CpG.train, s=best_lambda)

#plot actual by predicted values for training data
actual.fitted <- data.frame(fitted.lob, train.age)
View(actual.fitted)

tiff(file = "CPGandCHGelasticnettrain.tiff", res = 800, width = 100, height = 120, units = "mm")
plot(actual.fitted$s1~actual.fitted$train.age, xlab="chronological age", ylab="predicted age")
abline(lm(actual.fitted$s1~actual.fitted$train.age), col="red") # regression line (y~x)
dev.off()



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


##retrieve predicted ages
actual.fitted$s1


##get the error associated with each predicted ages
train_predicted_error <- (actual.fitted$s1 - actual.fitted$train.age)
#View(train_predicted_error)
train_predicted_error

trainerror<-c(0.42657050, 0.71836324, 0.28493131, 0.28159030, 0.08358966, 0.04722401,  0.05615927,  0.84306411,
              0.28941865,  0.70492330,  0.71739142,  0.16925885, 0.63424705, 0.66600783,  0.85870162,  1.24317272,
              0.05311685,  0.13732383, 1.82377302)
View(as.data.frame(trainerror))
mean(trainerror)
#mean is 0.6175581 years
Se<-sd(trainerror)/sqrt(21)
Se#se is 0.077, so error is 0.61755(+/- 0.077)



#predicted values for test set
predicted <- predict(lob.model, CpG.test, s=best_lambda)
predicted
test.age

#plot actual by predicted values for training data
actual.predicted <- data.frame(predicted, test.age)
actual.predicted

tiff(file = "CPGandCHGelasticnettest.tiff", res = 800, width = 100, height = 120, units = "mm")
plot(actual.predicted$s1~actual.predicted$test.age, xlab="chronological age", ylab="predicted age")
##fit line
abline(lm(actual.predicted$s1~actual.predicted$test.age), col="red") # regression line (y~x)
#note that this model is predicting that the 3 test trees that were all chronologically 37 years old
#however this is predicting that they are 72, 48, 38 years old
dev.off()




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


top5trainage<-top5trainage + scale_y_continuous(limits = c(0,90))
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
tiff(file="trainmodelcpgtop10.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



#calculate R2 for fitted data
R2.test <- cor(predicted, test.age)^2
R2.test




##get the error associated with each predicted ages
test_predicted_error <- (actual.predicted$s1 - actual.predicted$test.age)

test_predicted_error

testerror<-c(14.062292,   9.285718,   4.461703, 24.651991, 60.180660)
View(as.data.frame(testerror))
mean(testerror)
#mean is 22.53
Se<-sd(testerror)/sqrt(5)
Se
#Se is 9.99



#Doing pearson correlations and linear models

#need to make a train set to obtain the correlated sites, then run on the test individuals
train<-c(1:13,15,17:19,21,23)
#test<-c(14,16,20,22,24)
dim(Lob1[c(train),])

Lobtrain<-Lob1[c(train),]
dim(Lobtrain)

corrtest1 <- corr.test(Lobtrain[-1], Lobtrain$AGE, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p.adj
View(pval)
lowestpval<-sort(pval, decreasing = F)
lowestpval
#0 sites significant


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
write.csv(correlated, file = "traincorrelatedCpg and CHG.csv")

#using 5 sites to estimate age:

#testing if doing pearson on only train model is same as test model:

dim(Lob1)
a<-Lob1$'21133'
b<-Lob1$'44243'
c<-Lob1$'2070'
d<-Lob1$'15746'
e<-Lob1$'32913'

#dim(Lob2)
pearsonmodel<-data.frame(a,b,c,d,e)
pearsonmodel


age<-Lob1[c(1)]
age


pearsonmodel<-data.frame(pearsonmodel, age)
pearsonmodel#need to divide this into train and test values and use model from train to predict test
#divide the pearsonmodel into train (1-13, 15, 17-19, 21,23)
#divide into test (14,16,20,22,24)
train<-c(1:13,15,17:19,21,23)
test<-c(14,16,20,22,24)
#train
#test
trainmodel<-pearsonmodel[c(train),]
dim(trainmodel)
testmodel<-pearsonmodel[c(test),]
dim(testmodel)
#testmodel<-testmodel[-c(10)]
trainmodel
testmodel
results<-lm(AGE ~ a+b+c+d+e, data = trainmodel)
summary(results)


predict(results, newdata = trainmodel)
trainmodel$AGE



chrono<-c(46, 46, 46, 37, 37, 37, 19, 19, 19, 10, 10, 10, 55, 55,  1,  1, 28, 28, 97)
pred<-c(39.775642, 33.175094, 49.607343, 41.083218, 39.571088, 35.162351, 15.687842, 22.070613, 19.547322, 10.755343,  7.929871, 7.743911, 55.977577, 53.753970,  7.918599,  7.985841, 20.598158, 34.416548, 98.239670)
error<-chrono - pred
error

error<-c(6.224358, 12.824906, 3.607343, 4.083218, 2.571088,  1.837649,  3.312158, 3.070613, 0.547322, 0.755343,  2.070129,
         2.256089, 0.977577,  1.246030, 6.918599, 6.985841,  7.401842, 6.416548, 1.239670)
View(as.data.frame(error))
mean(error)
#so mean error was 3.912
SE<-sd(error)/sqrt(19)
SE#st err was 1.70. 

error
length(error)
dim(trainmodel)
trainmodel<-data.frame(trainmodel, pred, error)
trainmodel


top5trainage<-ggplot(trainmodel, aes(x = AGE, y = pred)) + 
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
tiff(file="trainmodelcpGandCHGtop5r2.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



predict(results, newdata = testmodel)
testmodel$AGE


test.age
chrono<-c(55,   1,  28,  82, 119)
pred<-c(40.16550, 13.32552, 32.53266, 61.15256, 60.65425)
error<-chrono - pred
error

error<-c(14.83450, 12.32552,  4.53266,  20.84744,  58.34575)
View(as.data.frame(error))
mean(error)
#so mean error was 22.18
SE<-sd(error)/sqrt(5)
SE#st err was 9.41 


error
length(error)
dim(trainmodel)
testmodel<-data.frame(testmodel, pred, error)
testmodel

results<-lm(pred~AGE, data = testmodel)
summary(results)



top5trainage<-ggplot(testmodel, aes(x = AGE, y = pred)) + 
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
tiff(file="testmodelcpgandCHGtop5r2.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()


#top 10
a<-Lob1$'21133'
b<-Lob1$'44243'
c<-Lob1$'2070'
d<-Lob1$'15746'
e<-Lob1$'32913'
f<-Lob1$'12447'
g<-Lob1$'33652'
h<-Lob1$'42781'
i<-Lob1$'47299'
j<-Lob1$'31522'
#dim(Lob2)
pearsonmodel<-data.frame(a,b,c,d,e,f,g,h,i,j)
pearsonmodel


age<-Lob1[c(1)]
age


pearsonmodel<-data.frame(pearsonmodel, age)
pearsonmodel#need to divide this into train and test values and use model from train to predict test
#divide the pearsonmodel into train (1-13, 15, 17-19, 21,23)
#divide into test (14,16,20,22,24)
train<-c(1:13,15,17:19,21,23)
test<-c(14,16,20,22,24)
#train
#test
trainmodel<-pearsonmodel[c(train),]
dim(trainmodel)
testmodel<-pearsonmodel[c(test),]
dim(testmodel)
#testmodel<-testmodel[-c(10)]
trainmodel
testmodel
results<-lm(AGE ~ a+b+c+d+e+f+g+h+i+j, data = trainmodel)
summary(results)


predict(results, newdata = trainmodel)
trainmodel$AGE

#so this model is predicting that
#  55, 1, 28, 82, 119 yr old trees are
#35.23912, 4.560616, 21.063413, 86.93362, and 101.83 years old

chrono<-c(46, 46, 46, 37, 37, 37, 19, 19, 19, 10, 10, 10, 55, 55,  1,  1, 28, 28, 97)
pred<-c(44.146354, 41.449049, 48.526427, 38.320032, 40.076494, 37.344722, 16.240878, 22.717274, 17.611861, 10.741743,  6.553634, 2.799928, 53.254050, 55.557032,  7.688415,  6.635172, 27.764015, 26.595263, 96.977657)
error<-chrono - pred
error

error<-c(1.853646,  4.550951, 2.526427, 1.320032, 3.076494, 0.344722,  2.759122, 3.717274,  1.388139, 0.741743,  3.446366,
         7.200072,  1.745950, 0.557032, 6.688415, 5.635172,  0.235985,  1.404737,  0.022343)
View(as.data.frame(error))
mean(error)
#so mean error was 2.59
SE<-sd(error)/sqrt(19)
SE#st err was 0.5


predict(results, newdata = testmodel)



test.age
chrono<-c(55,   1,  28,  82, 119)
pred<-c(32.67019, 12.09077, 25.57573, 67.04944, 64.18188)
error<-chrono - pred
error
View(as.data.frame(error))
error<-c(22.32981, 11.09077, 2.42427, 14.95056, 54.81812)
mean(error)
#so mean error was 21.12
SE<-sd(error)/sqrt(5)
SE#st err was 9.01 

error
length(error)
dim(trainmodel)
testmodel<-data.frame(testmodel, pred, error)
testmodel

results<-lm(pred~AGE, data = testmodel)
summary(results)








