library(dplyr)
library(reshape2)
library(ggplot2)
library(knitr)




# ~ 060614-Dong-Standard(Igll5-Znf80)+sample  ------------------------------------------------------------
# Read in data
data <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/Data/Raw data/060614-Dong-Standard(Igll5-Znf80)+samples.csv", skip = 1, sep = ",", header = T)
# name <- readLines("/Volumes/Untitled/cnv/arq2.txt", n=2)[2]
# names <- strsplit(name, split ="\t")[[1]]
# colnames(arq) <- make.names(colnames(arq))

## @knitr Data processing
s <- tbl_df(data)

s.sample <- filter(s, !grepl("[Ne][Ee][Gg]*", s$Sample.Name))
s.neg <- filter(s, grepl("[Ne][Ee][Gg]*", s$Sample.Name))
s <- s.sample

std <- filter(s, Sample.Name %in% c("1-1", "1-2", "1-3", "2-1", "2-2", "2-3"))

if (dim(std)[1] == 0){
        sample <- filter(s, !Sample.Name %in% c("1-1", "1-2", "1-3", "2-1", "2-2", "2-3") & !Sample.Name == c("ST"))
        
        sample.su <- group_by(sample, Pos, Sample.Name) %.% 
                summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %.%
                group_by(Sample.Name) %.%
                summarize(dCt.mean = mean(dCt), dCt.sem = sd(dCt)/sqrt(length(dCt)),
                          RelQ.mean = mean(RelQ), RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))
        
} else{
        std.su <- group_by(std, Pos, Sample.Name) %.% 
                summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %.%
                group_by(Sample.Name) %.%
                summarize(dCt.mean = mean(dCt), dCt.sem = sd(dCt)/sqrt(length(dCt)),
                          RelQ.mean = mean(RelQ), RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))
        std0 <- std.su[1:3,]
        # std1 <- std.su[4:6,]
        
        sample <- filter(s, !Sample.Name %in% c("1-1", "1-2", "1-3", "2-1", "2-2", "2-3") & !Sample.Name == c("ST"))
        
        sample.su <- group_by(sample, Pos, Sample.Name) %.% 
                summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %.%
                group_by(Sample.Name) %.%
                summarize(dCt.mean = mean(dCt), dCt.sem = sd(dCt)/sqrt(length(dCt)),
                          RelQ.mean = mean(RelQ), RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))
}

kable(sample.su, align = "c", digits = 2)


## @knitr Standard curve - table
if (!exists("std0")){
        std.curve0 <- data.frame(ratio = c("1-1","1-2","1-3"),
                         dCt = c(1.2000464, 0.9563856, 0.6797154),
                         RelQ = c(0.4354689, 0.5153470, 0.6243136),
                         copy.number = c(2,4,6))
        std.lm0 <- lm(copy.number ~ RelQ , data = std.curve0)
        } else {
                std.curve0 <- data.frame(ratio = std0$Sample.Name, dCt = std0$dCt.mean, RelQ = std0$RelQ.mean, copy.number = c(2,4,6))
                std.lm0 <- lm(copy.number ~ RelQ , data = std.curve0)
        }

kable(std.curve0, align = "c", digits = 2)


# std.curve1 <- data.frame(ratio = std1$Sample.Name, dCt = std1$dCt.mean, RelQ = std1$RelQ.mean, copy.number = c(1,2,3))
# std.lm1 <- lm(copy.number ~ RelQ , data = std.curve0)  


# plot(copy.number ~ RelQ, data = std.curve0)
# abline(std.lm0)

## @knitr Standard.curve.plot
ggplot(data = std.curve0, aes(y = RelQ, x = as.numeric(copy.number)))+
        geom_point(size = 3) +
        labs(x = "Copy Number", y = "Relative Amount", title = "Standard curve of IGLL5")+
        geom_smooth(method = 'lm') +
        ylim(0, 0.8)
#         geom_point(aes(x = round(predicted.CN, 2), y = 0))
#         xlim(0,0.7)+
#         ylim(0,7)

predicted.CN <- as.numeric(predict(std.lm0, data.frame(RelQ = sample.su$RelQ.mean)))


## @knitr Copy.number.prediction.table
sample.su$Copy.Nbs = round(predict(std.lm0, 
                                    data.frame(RelQ = sample.su$RelQ.mean)))

sample.su$Copy.Nbs.pred = predict(std.lm0, 
                                     data.frame(RelQ = sample.su$RelQ.mean))

kable(sample.su, align = "c", digits = 2)

diff <- abs(sample.su$Copy.Nbs.pred-sample.su$Copy.Nbs)

## @knitr Copy.number.prediction.plot
std.lm0.t <- lm(RelQ ~ copy.number, data = std.curve0)
ggplot(data = sample.su, aes(x = RelQ.mean, fill = Sample.Name))+
        geom_histogram()+
        geom_vline(xintercept = predict(std.lm0.t, 
                                        data.frame(copy.number = c(1,2,3,4,5,6))))+
        theme(legend.position = "bottom")























## @knitr
s060514.su <- group_by(s060514, Sample.Name, Type, Target.Name) %.% 
        summarize(Cp.mean = mean(Cp), RQ = 2^-(Cp.mean))

s060614.su <- filter(tbl_df(s060614), Sample.Name %in% c("1-1", "1-2", "1-3", "2-1", "2-2", "2-3")) %.%
        group_by(Sample.Name, Type, Target.Name) %.%
        summarize(Cp.mean = mean(Cp), RQ = 2^-Cp.mean) 


ggplot(data = s060514.su, aes(x = Sample.Name, y = RQ, fill = Type)) +
        geom_bar(stat = "identity", position = "dodge")

plot(s060514.std$Cp, s060614.std$Cp)

ggplot(data = s060614.su, aes(x = Sample.Name, y = RQ, fill = Type)) +
        geom_bar(stat = "identity", position = "dodge")


Tg.s5 <- filter(s060514.su, Type == "Reference Unknown")
Tg.s6 <- filter(s060614.su, Type == "Reference Unknown")

Tg.s5 <- filter(s060514.su, Type == "Target Unknown")
Tg.s6 <- filter(s060614.su, Type == "Target Unknown")

plot(Tg.s5$RQ, Tg.s6$RQ)
text(Tg.s5$RQ, Tg.s6$RQ, labels = Tg.s5$Sample.Name)

# compare btw different ratios
ggplot(data = s060514.su,
       aes(x = Sample.Name, y = RQ, fill = Target.Name)) +
        geom_bar(position = "fill",
                 stat = "identity") +
        geom_text(label = c("605", "605",
                            "590", "590",
                            "187", "187",
                            "180", "180", 
                            "471", "471",
                            "487", "487"))


ggplot() +
        geom_bar(data = Tg.s5, position = "stack",
                 aes(x = Sample.Name, y = RQ),
                 stat = "identity")


+
        geom_bar(stat = "identity", 
                 position = "dodge",
                 data = Tg.s6, 
                 aes(x = Sample.Name, y = log10(RQ),
                     fill = Sample.Name,
                     alpha = 0.3)) +
        
        geom_text(label = round(Tg.s5$RQ,3), vjust = 5)
        
        
        


# ~ test000 ~  ===========
# Read in data ------------------------------------------------------------
arq <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/Data/Raw data/arq2.csv", skip = 1, sep = ",", header = T)
# name <- readLines("/Volumes/Untitled/cnv/arq2.txt", n=2)[2]
# names <- strsplit(name, split ="\t")[[1]]
# colnames(arq) <- make.names(colnames(arq))

# Data preprocessing ------------------------------------------------------

arq <- tbl_df(arq)


arq.su1 <- group_by(arq, Pos, Sample.Name) %.% 
        summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %.%
        group_by(Sample.Name) %.%
        summarize(dCt.mean = mean(dCt), dCt.sem = sd(dCt)/sqrt(length(dCt)),
                  RelQ.mean = mean(RelQ), RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))
arq.su1 <- arq.su1[1:3,]

# Basic plot --------------------------------------------------------------
quartz()
par(mfrow = c(1,2))
par(mar = c(10,4,4,4))
barplot(arq.su2[-4,]$dCt, ylim = c(0, 2),
        col = unclass(arq.su2$Sample.Name) + 1,
        ylab = "dCt",
        names.arg = c("Amplified DNA (50ng)", "Original DNA (50ng)", "Amplified DNA (150ng)"),
        las = 3)

barplot(arq.su2[-4,]$Relative.Amount, ylim = c(0, 0.5),
        col = unclass(arq.su2$Sample.Name) + 1,
        ylab = "2^(-dCt)",
        names.arg = c("Amplified DNA (50ng)", "Original DNA (50ng)", "Amplified DNA (150ng)"),
        las = 3)

# ggplot ------------------------------------------------------------------
# Ct
data <- arq.su1[-4,]
data1 <- melt(data) 
ggplot(data = data,
       aes(x = Sample.Name,
           y = Ct.mean, fill = Sample.Name)) +
        geom_bar(stat = "identity")+
        geom_errorbar(aes(ymin = Ct.mean - Ct.sem,
                          ymax = Ct.mean + Ct.sem), 
                      width = 0.4) +
        theme_bw()+
        labs(x = "" , y = "dCt value", title = "Copy number of IGLL5") +
        theme(axis.text.x = element_text(angle = 45, size = 10, 
                                         hjust = 1.1, vjust = 1.1))+
        theme(legend.text = element_text(size = 10))
        
# RelQ
ggplot(data = data,
       aes(x = Sample.Name,
           y = RelQ.mean, fill = Sample.Name)) +
        geom_bar(stat = "identity")+
        geom_errorbar(aes(ymin = RelQ.mean - RelQ.sem,
                          ymax = RelQ.mean + RelQ.sem), 
                      width = 0.4)+
        labs(x = "" , y = "Relative Quantification", title = "Copy number of IGLL5") +
        ylim()
        theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1.1, vjust = 1.1))+
        theme(legend.text = element_text(size = 10))

# Ct and RelQ
 

# Copy number prediction --------------------------------------------------
arq.su1$Copy.Number.pred = predict(std.lm0, 
                                   data.frame(RelQ = arq.su1$RelQ.mean))

arq.su1$Copy.Number = round(predict(std.lm0, 
                                    data.frame(RelQ = arq.su1$RelQ.mean)))







# ~ 060514-Dong-Standard(Igll5-Znf80) 
# Read in data ------------------------------------------------------------
standard <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/Data/Raw data/060514-Dong-Standard(Igll5-Znf80).csv", skip = 1, sep = ",", header = T)
s060514 <- standard
# name <- readLines("/Volumes/Untitled/cnv/arq2.txt", n=2)[2]
# names <- strsplit(name, split ="\t")[[1]]
# colnames(arq) <- make.names(colnames(arq))

# cps <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/Data/Raw data/060514-Dong-Standard(Igll5-Znf80).csv", skip = 1, sep = ",", header = T)
# cp <- filter(cps, Target.Name == "ZNF80") %.%
#         group_by(Sample.Name) %.%
#         summarize(cp.mean = mean(Cp))
# cp$conc = c(605, 590, 187, 180, 471, 293)
# plot(x = cp[1:3,]$cp.mean, y = cp[1:3,]$conc)
# plot(x = cp[4:6,]$cp.mean, y = cp[4:6,]$conc, 
#      xlim = c(5,10),
#      col = unclass(cp[4:6,]$Sample.Name)+2)
# text(labels=cp[4:6,]$Sample.Name, pos=4,
#      x = cp[4:6,]$cp.mean, 
#      y = cp[4:6,]$conc)
# legend(legend = as.character(cp[4:6,]$Sample.Name),
#        "topright",
#        pch = c(1,3,2),
#        col = unclass(cp[4:6,]$Sample.Name)+2)

CairoSVG("a.svg")
ggplot(data = cp[4:6,],
       aes(x = cp.mean, y = conc)) +
        geom_point(size = 3, aes(col = Sample.Name,
                   shape = Sample.Name))+
        theme_bw()+
        geom_line(linetype = 1, col = "blue")
dev.off()
library(Cairo)


# Data preprocessing ------------------------------------------------------
std <- tbl_df(standard)

std.su <- group_by(std, Pos, Sample.Name) %.% 
        summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %.%
        group_by(Sample.Name) %.%
        summarize(dCt.mean = mean(dCt), dCt.sem = sd(dCt)/sqrt(length(dCt)),
                  RelQ.mean = mean(RelQ), RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))
std0 <- std.su
# std1 <- std.su[1:3,]
# std2 <- std.su[4:6,]

# standard curve and copy number prediction
# std.curve1 <- data.frame(dCt = std1$dCt.mean, RelQ = std1$RelQ.mean, copy.number = c(2,4,6))
# std.lm1 <- lm(copy.number ~ RelQ , data = std.curve1)  
# plot(copy.number ~ RelQ, data = std.curve1)
# abline(std.lm1)
# 
# 
# std.curve2 <- data.frame(dCt = std2$dCt.mean, RelQ = std2$RelQ.mean, copy.number = c(1,2,3))
# std.lm2 <- lm(copy.number ~ dCt , data = std.curve2)  
# plot(copy.number ~ dCt, data = std.curve2)
# abline(std.lm2)

# Standard curve ----------------------------------------------------------
std.curve0 <- data.frame(ratio = std0$Sample.Name, dCt = std0$dCt.mean, RelQ = std0$RelQ.mean, copy.number = c(2,4,6))
std.lm0 <- lm(copy.number ~ RelQ , data = std.curve0)  
plot(copy.number ~ RelQ, data = std.curve0)
abline(std.lm0)

ggplot(data = std.curve0, aes(x = RelQ, y = copy.number))+
        geom_point()+
        labs(y = "Copy Number", title = "Copy typing of IGLL5")+
        geom_smooth(method = 'lm')+
        geom_point(data = s, aes(x = RelQ, y = copy.number, col = 'red'))+
        geom_smooth(data = s, method = "lm", col = "red")
#         xlim(0,0.7)+
#         ylim(0,7)
ggsave("standard-Igll5.svg")