source("http://bioconductor.org/biocLite.R")
biocLite("CNVtools")
library(CNVtools)
data(A112)
head(A112)

# read in data with pure signal info
raw.signal <- as.matrix(A112[, -c(1,2)])
dimnames(raw.signal)[[1]] <- A112$subject # add subject names as rownames

# get mean and pca scores (or x)
mean.signal <- apply(raw.signal, MAR = 1, FUN = mean)
pca.signal <- apply.pca(raw.signal)

pr <- prcomp(raw.signal)
pr$rotation[,1]
hist(pr$x[,1], breaks = 50)


# plot
quartz()
hist(mean.signal, breaks = 50, main = "Mean signal",
     cex.lab = 1.3)
quartz()
hist(pca.signal, breaks = 50, main = "First PCA signal",
     cex.lab = 1.3)

plot(mean.signal, type = "p")
barplot(mean.signal)


# 3. Model selection using the Bayesian Information Criterion (BIC)
batches <- factor(A112$cohort)
sample <- factor(A112$subject)
set.seed(0)
results <- CNVtest.select.model(signal=pca.signal, 
                                batch = batches, 
                                sample = sample, 
                                n.H0 = 3, 
                                method="BIC", 
                                v.ncomp = 1:5, 
                                v.model.component = rep('gaussian',5),
                                v.model.mean = rep("~ strata(cn)",5),
                                v.model.var = rep("~1", 5))
ncomp <- results$selected

plot(-results$BIC, xlab="n comp", 
     ylab="-BIC", type="b", 
     lty=2, col="red", pch = '+')

# 4 Clustering the PCA transformed data
ncomp <- 3
batches <- factor(A112$cohort)
sample <- factor(A112$subject)
fit.pca <- CNVtest.binary(signal = pca.signal,
                          sample = sample,
                          batch = batches,
                          ncomp = ncomp, 
                          n.H0=3, 
                          n.H1=0,
                          model.var= '~ strata(cn)')

quartz()
cnv.plot(fit.pca$posterior.H0, batch = '58C', main = 'Cohort 58C', breaks = 50, col = 'red') 
cnv.plot(fit.pca$posterior.H0, batch = 'NBS', main = 'Cohort NBS', breaks = 50, col = 'red')


# 5 Assigning individuals to a copy number genotype
head(fit.pca$posterior.H0)
hist(fit.pca$posterior.H0[,6])

# 6 Improving the fit using the LDF procedure
ncomp <- 3
pca.posterior <- as.matrix((fit.pca$posterior.H0)[, paste('P',seq(1:ncomp),sep='')]) 
dimnames(pca.posterior)[[1]] <- (fit.pca$posterior.H0)$subject
ldf.signal <- apply.ldf(raw.signal, pca.posterior)
hist(pca.signal, breaks=50, main='First PCA signal', cex.lab=1.3)
hist(ldf.signal, breaks=50, main='LDF signal', cex.lab=1.3)


# 7 Testing for genetic association with a dichotomous disease trait
# 7.1 Some mathematical details
# 7.2 Example
ncomp <- 3
trait <- ifelse( A112$cohort == '58C', 0, 1)
fit.ldf <- CNVtest.binary (signal = ldf.signal, 
                           sample = sample, 
                           batch = batches, 
                           disease.status = trait, 
                           ncomp = ncomp, 
                           n.H0=3, n.H1=1, 
                           model.var = "~cn")
print(fit.ldf$status.H0)
print(fit.ldf$status.H1)

cnv.plot(fit.ldf$posterior.H0, batch = '58C', main = 'Cohort 58C', breaks = 50, col = 'red')
cnv.plot(fit.ldf$posterior.H0, batch = 'NBS', main = 'Cohort NBS', breaks = 50, col = 'red')

LR.statistic <- -2*(fit.ldf$model.H0$lnL - fit.ldf$model.H1$lnL)
print(LR.statistic)

# factor(cn)
fit.ldf <- CNVtest.binary (signal = ldf.signal, 
                           sample = sample, 
                           batch = batches, 
                           disease.status = trait, 
                           ncomp = 3, 
                           n.H0=3, 
                           n.H1=1, 
                           model.disease = " ~ as.factor(cn)")
print(fit.ldf$status.H0)
print(fit.ldf$status.H1)
LR.statistic <- -2*(fit.ldf$model.H0$lnL - fit.ldf$model.H1$lnL)
print(LR.statistic)


# 8 Testing for genetic association with a quantitative trait
batches <- rep("ALL", length(sample))
qt <- rnorm(length(sample), mean = 9.0 ,sd = 1.0)
fit.ldf <- CNVtest.qt(signal = ldf.signal,
                      sample = sample,
                      batch = batches,
                      qt = qt,
                      ncomp = ncomp,
                      n.H0 = 3,
                      n.H1 = 3,
                      model.var = "~strata(cn)")
print(fit.ldf$status.H0)
print(fit.ldf$status.H1)
LR.statistic <- -2*(fit.ldf$model.H0$lnL - fit.ldf$model.H1$lnL)
print(LR.statistic)
qt.plot(fit.ldf)




# - -=-------
library(dplyr)
library(ggplot2)
cont.5row <- read.table("~/Documents/R project/Memphis - Copy Number Variation/RIVUR/RIVUR aCGH/Cont_Case DNA/5000A_528524(5)_segMNT.txt",
                        header = T,
                        nrow = 5)
classes <- sapply(cont.5row, class)
cont <- read.table("~/Documents/R project/Memphis - Copy Number Variation/RIVUR/RIVUR aCGH/Cont_Case DNA/5000A_528524(5)_segMNT.txt",
                   header = T,
                   colClasses = classes)
cont <- tbl_df(cont)

hist(cont$GC)
plot(cont$REF_5000A_528524.5._635, cont$EXP_5000A_528524.5._532)
 
quartz()
ggplot(data = cont, aes(x = CHROMOSOME, y = RATIO_CORRECTED, col = CHROMOSOME))+
        geom_jitter()
hist(cont$GC)

library(reshape2)
melt(filter(cont, RATIO_CORRECTED == min(RATIO_CORRECTED)))
