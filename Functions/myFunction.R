
f <- function(data = data, keyMap_std = "default"){
        # library
        library(knitr)
        library(xtable)
        library(dplyr)
        library(ggplot2)
        library(reshape2)
        
        # key-map
        map <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/IGLL5/Map/RIVUR_DILUTION_PLATES.csv",sep = ",", header = T)
        p384 <- sapply(LETTERS[1:16], function(x) paste0(x, 1:24))%>%t
        p96 <- sapply(LETTERS[1:8], function(x) paste0(x, 1:12)) %>%t
        
        v96 = vector("character")
        v96_l <- lapply(seq_along(p96[,1]), function(i) {
                a <- paste(paste(p96[i,],p96[i,]), p96[i,]) 
                b <- strsplit(a, split = " ") %>% unlist
                b[37:48] = NA
                v96 = c(v96, b)
        }) 
        ## make a 384-96 map 
        keyMap <- data.frame(p96 = unlist(v96_l), p384 = as.vector(t(p384)))
        if (keyMap_std == "default"){
                keyMap_std <- data.frame(pos = c("B13", "B14", "B15",
                                                 "D13", "D14", "D15",
#                                                  "I17", "I18", "I16",                        
                                                 "F13", "F14", "F15",
                                                 "H13", "H14", "H15",
                                                 "J13", "J14", "J15",
                                                 "L13", "L14", "L15",
                                                 "N13", "N14", "N15"),
                                         sample = c("K106957-", "K106957-", "K106957-",
                                                    "K44781-", "K44781-", "K44781-",
                                                    "K96097-", "K96097-", "K96097-",
                                                    "K94705-", "K94705-", "K94705-",
                                                    "K84083-", "K84083-", "K84083-",
                                                    "K75881-", "K75881-", "K75881-",
                                                    "K94705B-", "K94705B-", "K94705B-"))
        }
       
        ## retrieve info form key map
        map = tbl_df(map)
        map$POSITION = strsplit(as.character(map$POSITION), split = " ")%>%unlist # clean whitespace
        plate_code = unique(data$plate_ID)
        map_key = filter(map, Plate_Box_Inv_Code %in% plate_code)

        # test 384-well plate data info
        data = tbl_df(data)
        data_su = group_by(data, Pos, Sample.Name, plate_ID) %>%
                summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)
        
        # merge key map and plate data
        ## data info
        data_384 <- merge(keyMap, data_su, by.x = "p384", by.y = "Pos", all = T)%>% tbl_df
        
        data_merge <- merge(data_384, map_key, by.x = "p96", by.y = "POSITION", all.x = T, sort = T ) %>% tbl_df
        
        data_m_su <- group_by(data_merge, RUID, plate_ID) %>% 
                summarize(ddCt.mean = mean(ddCt), 
                          ddCt.sem = sd(ddCt)/sqrt(length(ddCt)), 
                          RelQ.mean = mean(RelQ), 
                          RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
                          conc = mean(Concentration_ng._L))  %>%
                filter(RUID != "NA") %>%
                select(-conc)
        
        ## standard info
        std_merge <- merge(keyMap_std, data_su, by.x = "pos", by.y = "Pos", all.x = T)
        
        std_m_su <- group_by(std_merge, sample, plate_ID) %>% 
                summarize(ddCt.mean = mean(ddCt), 
                          ddCt.sem = sd(ddCt)/sqrt(length(ddCt)), 
                          RelQ.mean = mean(RelQ), 
                          RelQ.sem = sd(RelQ)/sqrt(length(RelQ))) %>%
                arrange(plate_ID)
        
        # linear regression and copy number estimation
        std_m_su$copy.Nbs = c(NA, 1, NA, 2, NA, NA, NA)
        std.fit <- lm(copy.Nbs ~ RelQ.mean, data = std_m_su)
        
        data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
        data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
        diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
        data_m_su$diff <- diff
#         file <- paste(c("IGLL5-", plate_code, ".csv"), collapse = "")
#         write.csv(file = file, data_m_su)
        
        result = list()
        result$data_su = tbl_df(data_m_su)
        result$std_su = tbl_df(std_m_su)
        result$diff = diff
        result$sampleLoc = select(data_merge, c(1,2,4,8))
        re = result
}

# f(data)$sampleLoc


checkLoc <- function(RUID){
        map <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/IGLL5/Map/RIVUR_DILUTION_PLATES.csv",sep = ",", header = T)
        p384 <- sapply(LETTERS[1:16], function(x) paste0(x, 1:24))%>%t
        p96 <- sapply(LETTERS[1:8], function(x) paste0(x, 1:12)) %>%t
        
        v96 = vector("character")
        v96_l <- lapply(seq_along(p96[,1]), function(i) {
                a <- paste(paste(p96[i,],p96[i,]), p96[i,]) 
                b <- strsplit(a, split = " ") %>% unlist
                b[37:48] = NA
                v96 = c(v96, b)
        }) 
        ## make a 384-96 map 
        keyMap <- data.frame(p96 = unlist(v96_l), p384 = as.vector(t(p384)) )
        km = keyMap
}

### example #######################################
# file_list <- list()
# fileloc = "/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/IGLL5/Raw data/"
# fi <- list.files(fileloc)
# for (file in fi){
#         file_list[[length(file_list)+1]] <- read.table(paste0(fileloc, file),skip = 1, sep = ",", header = T)
# } 
# 
# for(i in 1:length(fi)){
#         file_list[[i]]$plate_ID = strsplit(fi[i], split = "-")[[1]][3]
#         attr(file_list[[i]], "plate_ID") = strsplit(fi[i], split = "-")[[1]][3]
# }
# 
# 
# data_batch <- lapply(file_list, function(x) f(x)$data_su)
# std_batch <- lapply(file_list, function(x) f(x)$std_su)
# 
# data_all <- do.call(rbind, data_batch)
# std_all <- do.call(rbind, std_batch)
# table(data_batch$Copy.Nbs.round)
# plot(table(data_batch$Copy.Nbs.round))
### example #######################################

makeForestPlotForRCTs <- function(mylist, referencerow=2){
        require("rmeta")
        numstrata <- length(mylist)
        # make an array "ntrt" of the number of people in the exposed group, in each stratum
        # make an array "nctrl" of the number of people in the unexposed group, in each stratum
        # make an array "ptrt" of the number of people in the exposed group that have the disease,
        # in each stratum
        # make an array "pctrl" of the number of people in the unexposed group that have the disease,
        # in each stratum
        ntrt <- vector()
        nctrl <- vector()
        ptrt <- vector()
        pctrl <- vector()
        if (referencerow == 1) { nonreferencerow <- 2 }
        else                   { nonreferencerow <- 1 }
        for (i in 1:numstrata)
        {
                mymatrix <- mylist[[i]]
                DiseaseUnexposed <- mymatrix[referencerow,1]
                ControlUnexposed <- mymatrix[referencerow,2]
                totUnexposed <- DiseaseUnexposed + ControlUnexposed
                nctrl[i] <- totUnexposed
                pctrl[i] <- DiseaseUnexposed
                DiseaseExposed <- mymatrix[nonreferencerow,1]
                ControlExposed <- mymatrix[nonreferencerow,2]
                totExposed <- DiseaseExposed + ControlExposed
                ntrt[i] <- totExposed
                ptrt[i] <- DiseaseExposed
        }
        names <- as.character(seq(1,numstrata))
        myMH <- meta.MH(ntrt, nctrl, ptrt, pctrl, conf.level=0.95, names=names)
        print(myMH)
        tabletext<-cbind(c("","Study",myMH$names,NA,"Summary"),
                         c("Disease","(exposed)",ptrt,NA,NA),
                         c("Disease","(unexposed)",pctrl, NA,NA),
                         c("","OR",format(exp(myMH$logOR),digits=2),NA,format(exp(myMH$logMH),digits=2)))
        print(tabletext)
        m<- c(NA,NA,myMH$logOR,NA,myMH$logMH)
        l<- m-c(NA,NA,myMH$selogOR,NA,myMH$selogMH)*2
        u<- m+c(NA,NA,myMH$selogOR,NA,myMH$selogMH)*2
        forestplot(tabletext,m,l,u,zero=0,is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
                   clip=c(log(0.1),log(2.5)), xlog=TRUE,
                   col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))

#  example
# mymatrix1 <- matrix(c(198,728,128,576),nrow=2,byrow=TRUE)
# mymatrix2 <- matrix(c(96,437,101,342),nrow=2,byrow=TRUE)
# mymatrix3 <- matrix(c(1105,4243,1645,6703),nrow=2,byrow=TRUE)
# mymatrix4 <- matrix(c(741,2905,594,2418),nrow=2,byrow=TRUE)
# mymatrix5 <- matrix(c(264,1091,907,3671),nrow=2,byrow=TRUE)
# mymatrix6 <- matrix(c(105,408,348,1248),nrow=2,byrow=TRUE)
# mymatrix7 <- matrix(c(138,431,436,1576),nrow=2,byrow=TRUE)
# mylist <- list(mymatrix1,mymatrix2,mymatrix3,mymatrix4,mymatrix5,mymatrix6,mymatrix7)
# 
# makeForestPlotForRCTs(mylist)
}




#  forestplot ####################
forestplot <- function(d, xlab="Odds Ratio", ylab="Study"){
        require(ggplot2)
        p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi)) + 
                geom_pointrange() + 
                coord_flip() +
                geom_hline(aes(x=0), lty=2) +
                ylab(xlab) +
                xlab(ylab) +#switch because of the coord_flip() above
                theme_science
        return(p)
}

# d is a data frame with 4 columns
# d$x gives variable names
# d$y gives center point
# d$ylo gives lower limits
# d$yhi gives upper limits


# Create some dummy data.
d <- data.frame(x = toupper(letters[1:10]),
                y = rnorm(10, .05, 0.1))
d <- transform(d, ylo = y-1/10, yhi=y+1/10)

forestplot(d)





