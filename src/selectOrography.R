#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameters --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wkpath <- "C:/Users/dival/Documents/Tesis/Structural change/"
setwd(wkpath)
outdir <- "dmel_5kb/"
species <- "D. melanogaster"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Libraries ---------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(dplyr)
library(tidyverse)

# for vis:
library(plotly)
library(ggpubr)
library(scico)
library(moonBook)
library(webr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Auxiliary functions------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

look_left <- function(vect, m){

  if(m == 1){return(0)}
          
  i <- m -1
  flag <- FALSE
  
  while(vect[i] == "Majority" & i >= 1){
    i <- i - 1
  }
  
  if(vect[i] == "Valley"){
    flag <- TRUE
  }
  
  ifelse(flag,
         return(i),
         return(0))
}

look_right <- function(vect, m, n){
  i <- 1
  flag <- FALSE
  nReal <- n - m
  
  while(vect[i] == "Majority" & i <= nReal){
    i <- i + 1
  }
  
  if(vect[i] == "Valley"){
    flag <- TRUE  
  }
  
  ifelse(flag,
         return(m + i - 1),
         return(0))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading and wrangling ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orography_ <- read.csv(paste0(outdir, "orography_out.bed"), 
                       header = FALSE, sep = "\t") 

#head(orography_)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Select mountains ---------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chrs <- unique(orography_$V1)
newOrography <- data.frame()
interestRegions <- data.frame()


for (chr in chrs){
    
    orography <- orography_ %>% filter(V1 == chr) 
    orography
    
    mPos <- which(orography$V4 == "Majority")
    mPosTypes <- c()
    n <- dim(orography)[1]
    
    for(m in mPos){
        type <- 0
        L <- look_left(orography$V4[1:m], m)
        ifelse(L > 0, 
               type <- type + 1,
               L <- m)
        
        R <- look_right(orography$V4[m:n], m, n)
        ifelse(R > 0,
               type <- type + 2, 
               R <- m)
        
        if (type > 0){
            newRegion <- orography %>% slice(L:R)
            interestRegions <- bind_rows(interestRegions, newRegion) %>%
                                   unique()
        }
        
        # if (type == 3){
        #   newRegion <- orography %>% slice(L:R)
        #   interestRegions <- bind_rows(interestRegions, newRegion) %>% 
        #     unique()
        # }
        
        
        mPosTypes <- c(mPosTypes, type)
    }
    
    orography$mType <- -1
    orography$mType[mPos] <- mPosTypes
    orography$mType <- factor(orography$mType, 
                              levels = c("-1","0", "1","2","3"),
                              labels = c("Valley|Hill","Between hills", 
                                         "Left valley", "Right valley",
                                         "Between valleys"))
    
    newOrography <- bind_rows(newOrography, orography)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Split interest regions --------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


write_tsv(interestRegions, 
          paste0(outdir,"wd_interestRegions.bed"),
          col_names = FALSE, quote = "none")

wellDefinedValleys <- interestRegions %>% filter(V4 == "Valley")
wellDefinedMountains <- interestRegions %>% filter(V4 == "Majority")

write_tsv(wellDefinedValleys,
          paste0(outdir,"wd_valleys.bed"),
          col_names = FALSE, quote = "none")

write_tsv(wellDefinedMountains,
          paste0(outdir,"wd_mountains.bed"),
          col_names = FALSE, quote = "none")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualize statistics ----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newOrography$V4 <- factor(newOrography$V4, 
                       levels = c("Valley", "Hills", "Majority"),
                       labels = c("Valleys", "Hills", "MV-Mountains"))


ggplot(newOrography, aes(x = V4 , fill = V4)) +
  geom_bar(position = "dodge") +
  labs(x = "Final orography classification",
       y = "Number of regions") +
  scale_fill_manual(values = c("#97b792",
                               "#875f79",
                               "#190d33")) +
  theme_pubclean() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8))


df <- newOrography %>% 
                mutate(size = V3 - V2) %>% 
                group_by(V4, mType) %>% 
                summarise(s = sum(size))


PieDonut(df, aes(pies = V4, donuts = mType, count = s),
         labelposition = 1,
         title = "% of genome in each class",
         r0 = 0.4,
         star = 3*pi/6
)


