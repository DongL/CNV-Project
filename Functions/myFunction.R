
f <- function(data = data, keyMap_std = "default", calib ){
#         print(keyMap_std)
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
                                                 "C13", "C14", "C15",
                                                 "F13", "F14", "F15",
                                                 "H13", "H14", "H15",
                                                 "J13", "J14", "J15",
                                                 "L13", "L14", "L15",
                                                 "N13", "N14", "N15"),
                                         sample = c(rep("K106957-", 3), 
                                                    rep("K44781-", 3),
                                                    rep("K96097-", 3),
                                                    rep("K94705-", 3),
                                                    rep("K84083-", 3),
                                                    rep("K75881-", 3),
                                                    rep("K94705B-", 3)),
                                         copy_Nbs = rep(c(NA, 1, NA, 2, NA, NA, NA), times = 3))
        }
       
        ## retrieve info form key map
        map = tbl_df(map)
        map$POSITION = strsplit(as.character(map$POSITION), split = " ")%>%unlist # clean whitespace
        plate_code = unique(data$plate_ID); print(plate_code)
        map_key = filter(map, Plate_Box_Inv_Code %in% plate_code)

        # test 384-well plate data info
        data = tbl_df(data)
        
        calib_df = filter(data, Pos %in% calib) %>% group_by(Pos, Sample.Name, plate_ID) %>% summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)
        calib_dCt = mean(calib_df$dCt)
        
        data_su = group_by(data, Pos, Sample.Name, plate_ID) %>%
                summarize(ddCt = Cp[1]-Cp[2]-calib_dCt, RelQ = 2^-ddCt)
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
#         print(std_m_su)
        # linear regression and copy number estimation
        t = length(keyMap_std$copy_Nbs)/3
        std_m_su$copy.Nbs = keyMap_std$copy_Nbs[1:t]
        std.fit <- lm(copy.Nbs ~ RelQ.mean, data = std_m_su)
        
        data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
        data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
        diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
        data_m_su$diff <- diff
        file <- paste(c("IGLL5-", plate_code, ".csv"), collapse = "")
        write.csv(file = file, data_m_su)
        
        result = list()
        result$data_su = tbl_df(data_m_su)
        result$std_su = tbl_df(std_m_su)
        result$diff = diff
        result$sampleLoc = select(data_merge, c(1,2,4,8))
        re = result
}



plt_triple <- function (va) {
        data = rva1a3_array[!is.na(rva1a3_array[va]), ]
        p = list()
        lapply(seq_along(va), function (i) plt_box(data, "DEFA1.A3.CN", va)) -> p[1]
        lapply(seq_along(va), function (i) plt_prop(data, "DEFA1.A3.CN", va)) -> p[2]
        lapply(seq_along(va), function (i) plt_ecdf(data, "DEFA1.A3.CN", va)) -> p[3]
        do.call(grid.arrange, c(p, list(nrow = 1, ncol =3)))
}

