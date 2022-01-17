setwd("~/Documents/Parrott Lab/Research/Loblolly pine experiment/CHG/mat/")

library(psych)
library(caret)
#check to see how to load in libraries on cluster
library(methylKit)
library(impute)

#CpG
#start with 10x, sites across 80 percent of samples covered
mymat1<-read.csv(file="matCHG10desfalse80.csv", header = T, row.names="ID")
head(mymat1)


#transforming data so you can call variables and it is in the right order
lob1<- as.data.frame(t(mymat1))
#head(lob1)
dim(lob1)
# 24 30729


library(psych)

#removing near zero variant sites
nzv <- nearZeroVar(lob1)
filteredlob1 <- lob1[, -nzv]
dim(filteredlob1)
#now 22264 sites (so removed 8465 sites, or ~27.55% of the sites)
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
#17 sites with pearsons correlations that are significant


#running with 5 individuals in test set, 19 in train


dim(Lob1)

lob_age <- Lob1[c(1)]
#View(lob_age) #pulled age for each sample
dim(lob_age)
#lob_age[c(4, 5, 6), c(1)]


dim(Lob1)


lob_CpG <- Lob1[c(2:22264)]#just pulled all the CpG sites (now columns), for every individual (1-28), and it gives the 
#amount of methylation present at each site (need to see what values mean) (go back to what generated "meth")

#View(lob_CpG)
dim(lob_CpG)
#so, it looks like its going to run model on how values at all 146918 sites vary with age prediction


test.CpG <- lob_CpG[c(14, 16, 20, 22, 24), c(1:22263)]#2615 if 5x file
#View(test.CpG)
dim(test.CpG)



train.CpG <- lob_CpG[c(1:13, 15, 17:19, 21, 23), c(1:22263)]#2615 if 5 x
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



#define the lambda parameter
best_lambda <- lob.cv.day$lambda.min
best_lambda #best lambda of CpG imputed data set is 0.9627225 when using n folds = 21 but using the differing ages


lob.model <- glmnet(as.matrix(CpG.train), 
                    as.matrix(train.age),
                    alpha=0.5, family="gaussian", nlambda=100)


#coefficients for fitted model
age_coeff.lob <- coef(lob.model, s=best_lambda)
dim(age_coeff.lob)
#dimensions are 21567

#tells me which sites have non zero coefficient  values, so I can find out which sites were used
which(age_coeff.lob!=0)   
age_coeff.lob
#51 sites


list<-which(age_coeff.lob!=0)  
length(list)


#model:
age_coeff.lob[c(list),]

testtable<-age_coeff.lob[c(list),]
testtable
write.csv(testtable)
#so sites are: 433, 2252, 2255, 3237, 3791, 3895, 6243, 6244, 6343, 6393, 6547, 6572, 6756, 8659, 9791,
#11641, 13718, 14677, 15196, 16261, 16264, 16265, 16797, 17018, 17019, 17722, 17723, 17724, 18131, 18194, 
#18601, 18690, 20018, 20779, 21620, 22061, 23180, 24421, 24834, 25168, 25598, 25943, 25971, 26204, 28205, 
#28213, 28362, 28363, 28468, 30700

list<-c(433, 2252, 2255, 3237, 3791, 3895, 6243, 6244, 6343, 6393, 6547, 6572, 6756, 8659, 9791,
        11641, 13718, 14677, 15196, 16261, 16264, 16265, 16797, 17018, 17019, 17722, 17723, 17724, 18131, 18194, 
        18601, 18690, 20018, 20779, 21620, 22061, 23180, 24421, 24834, 25168, 25598, 25943, 25971, 26204, 28205, 
        28213, 28362, 28363, 28468, 30700)

#generating info for this list
setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")

methfile<-read.csv(file = "methCHG10desfalse80.csv")

dim(methfile)

methfile[c(list),][,c(1:3)]

dim(methfile[c(list),][,c(1:3)])
library(stringr) 
library(utils)
elasticnetcHGtable<-methfile[c(list),][,c(1:3)]
elasticnetcHGtable$chr
elasticnetcHGtable$start
write.csv(elasticnetcHGtable, file = "ElasticnetCHGtable.csv")

file.exists("C:/Users/Gardn/Documents/Parrott lab/Research/loblolly pine experiment/CHG/mat/")
#write(, "C:/myfile.txt")
write.csv(elasticnetcHGtable, file = "ELasticCHG.csv", "C:/Users/Gardn/Documents/Parrott lab/Research/loblolly pine experiment/CHG/mat/")
elasticnetcHGtable
write.csv(elasticnetcHGtable, file = "C:/Users/Gardn/Documents/Parrott lab/Research/loblolly pine experiment/ElasticnetCHGtable.csv")



setwd("C:/Users/Gardn/Documents/Parrott lab/Research/loblolly pine experiment/CHG/mat/")



#fitted values
fitted.lob <- predict(lob.model, CpG.train, s=best_lambda)

#plot actual by predicted values for training data
actual.fitted <- data.frame(fitted.lob, train.age)
View(actual.fitted)


tiff(file = "CHGelasticnettrainmodel.tiff", res = 800, width = 100, height = 120, units = "mm")
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

##retrieve predicted ages
actual.fitted$s1

#residuals(actual_predictedage)

##get the error associated with each predicted ages
train_predicted_error <- (actual.fitted$s1 - actual.fitted$train.age)
#View(train_predicted_error)
train_predicted_error

trainerror<-c(0.55482032, 0.56971335, 0.40875187, 0.24963356,  0.06358454,  0.21398983,  0.06865517,  0.72826909,
              0.41173880,  0.52799134,  0.63030980,  0.35782423, 0.79550685, 0.69739194,  1.07253767,  1.10285037,
              0.17379827,  0.09851597, 1.82665064)
View(as.data.frame(trainerror))

mean(trainerror)
#mean is 0.6175581 years
Se<-sd(trainerror)/sqrt(19)
Se#se is 0.1, so error is 0.553(+/- 0.1)



#predicted values for test set
predicted <- predict(lob.model, CpG.test, s=best_lambda)
predicted
test.age
train.age
#plot actual by predicted values for training data
actual.predicted <- data.frame(predicted, test.age)
actual.predicted

tiff(file = "CHGelasticnettestmodel.tiff", res = 800, width = 100, height = 120, units = "mm")
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
tiff(file="testmodelcpgtop5.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



#calculate R2 for fitted data
R2.test <- cor(predicted, test.age)^2
R2.test
#R2 was 0.7905 (so not great) 

##get the error associated with each predicted ages
test_predicted_error <- (actual.predicted$s1 - actual.predicted$test.age)
#View(test_predicted_error)

test_predicted_error

testerror<-c(10.9239249,  18.8683301,   0.3270453, 27.3921540, 67.4970080)
View(as.data.frame(testerror))

mean(testerror)
#mean is 22.679
Se<-sd(testerror)/sqrt(5)
Se
#Se is 11.52, 


#Doing pearson correlations and linear models

train<-c(1:13,15,17:19,21,23)
test<-c(14,16,20,22,24)
dim(Lob1)
Lob3<-Lob1[c(train),]
dim(Lob3)
#head(Lob3[c(2),])
Lob3[c(1)]
Lob3$AGE

corrtest1 <- corr.test(Lob3[-1], Lob3$AGE, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=FALSE)
pval <- corrtest1$p
View(pval)
lowestpval<-sort(pval, decreasing = F)
View(lowestpval)
#0 sites significant
a<-corrtest1$r
a
head(a)
c<-corrtest1$p


correlated<-data.frame(a,c)
correlated
dim(correlated)
correlated

View(correlated)


#getting all correlated sites with correlation coefficients greater than 0.6
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

write.csv(correlated, file = "cHGcorrelatedpearsontrain.csv")


#top 5 
dim(Lob1)

a<-Lob1$'17723'
b<-Lob1$'6393'
c<-Lob1$'16261'
d<-Lob1$'9791'
e<-Lob1$'20779'

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
train
test
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

#so this model is predicting that
#  55, 1, 28, 82, 119 yr old trees are
#35.23912, 4.560616, 21.063413, 86.93362, and 101.83 years old

chrono<-c(46, 46, 46, 37, 37, 37, 19, 19, 19, 10, 10, 10, 55, 55,  1,  1, 28, 28, 97)
#pred<-c(41.380891, 45.678699, 47.560323, 34.283344, 52.399068, 40.350283, 19.508582, 31.314640, 16.968398,  7.012066, 5.565504,  5.204868, 55.849606, 47.987515,  5.830662,  5.854652, 20.614820, 27.335139, 90.300939)
pred<-c(39.623508, 44.799734, 53.061423, 31.310933, 49.265915, 47.312615, 19.729306, 29.972876, 20.572915,  6.906010,  5.208994,  4.868483, 57.347150, 51.337835,  6.280990,  5.980829, 19.784001, 22.498368, 85.138114)


error<-chrono - pred
error

#error<-c(4.619109,   0.321301,  1.560323,   2.716656, 15.399068,  3.350283,  0.508582, 12.314640,   2.031602,
#         2.987934,  4.434496,   4.795132,  0.849606,   7.012485,  4.830662,  4.854652,   7.385180,   0.664861,
#         6.699061)
error<-c(6.376492,   1.200266,  7.061423,   5.689067, 12.265915, 10.312615,  0.729306, 10.972876,  1.572915,   3.093990,
         4.791006,   5.131517,  2.347150,   3.662165,  5.280990,  4.980829,   8.215999,  5.501632,  11.861886)
View(as.data.frame(error))

mean(error)
#so mean error was 4.59 (5.88)
SE<-sd(error)/sqrt(19)
SE#st err was 0.91. So MAE is 5.48 (+/- 1.13) years. 
#5.84 (+/- 0.81) 

error
length(error)
dim(trainmodel)
trainmodel<-data.frame(trainmodel, pred, error)
trainmodel

results<-lm(pred~AGE, data = trainmodel)
summary(results)

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
tiff(file="trainmodelcHGtop5r2.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()


results<-lm(AGE ~ a+b+c+d+e, data = trainmodel)
summary(results)

testmodel
trainmodel
#model R2 was 0.8831 (multiple), adjusted was 0.8381
predict(results, newdata = testmodel)

#so this model is predicting that
#  55, 1, 28, 82, 119 yr old trees are
#35.23912, 4.560616, 21.063413, 86.93362, and 101.83 years old

chrono<-c(55, 1, 28, 82, 119)
pred<-c(56.50930, 16.33473, 33.23193, 50.64411, 68.14478)

error<-chrono - pred
error

error<-c(1.50930, 15.33473,  5.23193,  31.35589,  50.85522)
View(as.data.frame(error))
mean(error)
#so mean error was 20.86
SE<-sd(error)/sqrt(5)
SE#st err was 9.1  


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
tiff(file="testmodelcHGtop5r2.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



#top 10

dim(Lob1)

a<-Lob1$'17723'
b<-Lob1$'6393'
c<-Lob1$'16261'
d<-Lob1$'9791'
e<-Lob1$'20779'
f<-Lob1$'3237'
g<-Lob1$'28363'
h<-Lob1$'21620'
i<-Lob1$'28205'
j<-Lob1$'3895'
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
train
test
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
pred<-c(42.576628, 40.635457, 47.296023, 33.113843, 42.713377, 40.423446, 13.022667, 26.284420, 15.019367,  8.934020, 11.548554,  8.883518, 51.930746, 55.461908,  4.421438,  5.938553, 24.721242, 29.979333, 98.095459)

error<-chrono - pred
error

error<-c(3.423372,  5.364543, 1.296023,  3.886157, 5.713377, 3.423446,  5.977333, 7.284420,  3.980633,  1.065980,
         1.548554,  1.116482,  3.069254, 0.461908, 3.421438, 4.938553,  3.278758, 1.979333, 1.095459)
View(as.data.frame(error))
mean(error)
#so mean error was 3.28
SE<-sd(error)/sqrt(19)
SE#st err was 0.91. So MAE is 5.48 (+/- 1.13) years. 

error
length(error)
dim(trainmodel)
trainmodel<-data.frame(trainmodel, pred, error)
trainmodel

results<-lm(pred~AGE, data = trainmodel)
summary(results)


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
tiff(file="trainmodelcHGtop10r2.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()


results<-lm(AGE ~ a+b+c+d+e+f+g+h+i+j, data = trainmodel)
summary(results)

testmodel
trainmodel
#model R2 was 0.8831 (multiple), adjusted was 0.8381
predict(results, newdata = testmodel)

#so this model is predicting that
#  55, 1, 28, 82, 119 yr old trees are
#35.23912, 4.560616, 21.063413, 86.93362, and 101.83 years old

chrono<-c(55, 1, 28, 82, 119)
pred<-c(58.22663,  4.65836, 24.21102, 54.75558, 46.51309)

error<-chrono - pred
error

error<-c(3.22663, 3.65836,  3.78898, 27.24442, 72.48691)
View(as.data.frame(error))
mean(error)
#so mean error was 22.08
SE<-sd(error)/sqrt(5)
SE#st err was 5.38 So MAE is 15.46 (+/- 5.38) years. 

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
tiff(file="testmodelcHGtop10r2.tiff", res = 800, width = 189, height = 100, units = "mm")
top5trainage
dev.off()



#getting distributions of the top 5 sites with islands, shores, shelves or open seas

#elastic net clock sites

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


setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

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

#putting in sites for CHG elastic net 
Agesites<-c(433,
            2252,
            2255,
            3237,
            3791,
            3895,
            6243,
            6244,
            6343,
            6393,
            6547,
            6572,
            6756,
            8659,
            9791,
            11641,
            13718,
            14677,
            15196,
            16261,
            16264,
            16265,
            16797,
            17018,
            17019,
            17722,
            17723,
            17724,
            18131,
            18194,
            18601,
            18690,
            20018,
            20779,
            21620,
            22061,
            23180,
            24421,
            24834,
            25168,
            25598,
            25943,
            25971,
            26204,
            28205,
            28213,
            28362,
            28363,
            28468,
            30700)
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


#8 elastic net sites fell into islands!!

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


setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
CpGmeth<-read.csv(file = "methCHG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

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


#8 fall in shores

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

write.csv(Shoresites, file = "CHGelasticnetshores.csv")





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


#4 elastic net sites fell into shelves!!

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

write.csv(Shelfsites, file = "CHGelasticnetshores.csv")

Islandsites
Shoresites
Shelfsites


#Additional sites for combined model
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

setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
CpGmeth<-read.csv(file = "methCHG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

#may not need to use the first column (site ID), but Ill try it. If not, remove and use only 2 - 4
#note that 2 - 4 are the scaffold ID, the start and end
filepathpeaks<-CpGmeth[,c(1:4)]
filepathpeaks<-filepathpeaks[-1,]#removing top row
filepathpeaks




#putting in sites for CpG elastic net 
Agesites<-c(5002,5003,7132, 12714, 25057)
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

write.csv(Islandsites, file = "CHGislands.csv")




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

write.csv(Shoresites, file = "CHGelasticnetshores.csv")





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

write.csv(Shelfsites, file = "CHGelasticnetshores.csv")

Islandsites
Shoresites
Shelfsites



#Top 5 pearson sites

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


setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
#reading in the CPG sites from the meth files
CpGmeth<-read.csv(file = "methCHG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

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
Agesites<-c(17723,
            6393,
            16261,
            9791,
            20779)
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


#2 top 5 pearson  sites fell into islands!!

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


setwd("~/Parrott lab/Research/loblolly pine experiment/CHG/meth/")
CpGmeth<-read.csv(file = "methCHG10desfalse80.csv", header = FALSE,stringsAsFactors=FALSE)

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


#0 fall in shores

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

write.csv(Shoresites, file = "CHGelasticnetshores.csv")





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

write.csv(Shelfsites, file = "CHGelasticnetshores.csv")

Islandsites
Shoresites
Shelfsites

