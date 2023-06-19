#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameters --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wkpath <- "/homes/bierdepot/dulce/Documents/MasterDul/master-evo-in-3D/TRIFORCE/TAD-triforce/test"
setwd(wkpath)
outdir <- "/homes/bierdepot/dulce/Documents/MasterDul//master-evo-in-3D/TRIFORCE/TAD-triforce/test/out"
species <- "D. melanogaster"
args <- commandArgs(trailingOnly = TRUE)

# Working directories
wkpath <- args[1]
outdir <- args[2]
setwd(wkpath)
species <- args[3]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Libraries ---------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(dplyr)
library(tidyverse)

# for vis:
library(plotly)
library(ggpubr)
library(scico)
#library(moonBook)
#library(webr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Auxiliary functions------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

look_left <- function(vect, m){

  if(m == 1){return(0)}
          
  i <- m -1
  flag <- FALSE
  
  while(vect[i] == "Mountain_2" & i >= 1){
    i <- i - 1
  }
  
  if(vect[i] == "Valley_1"){
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
  
  while(vect[i] == "Mountain_2" & i <= nReal){
    i <- i + 1
  }
  
  if(vect[i] == "Valley_1"){
    flag <- TRUE  
  }
  
  ifelse(flag,
         return(m + i - 1),
         return(0))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading and wrangling ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orography_ <- read.csv(paste0(outdir, "/orography_out.bed"), 
                       header = FALSE, sep = "\t")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Select mountains ---------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chrs <- unique(orography_$V1)
newOrography <- data.frame()
interestRegions <- data.frame()

for (chr in chrs){
    orography <- orography_ %>% filter(V1 == chr) 
    orography
     
    #mPos <- which(orography$V4 == "Majority")
    mPos <- which(orography$V4 == "Mountain_2")
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

# All regions
write_tsv(interestRegions,
          paste0(outdir,"/wd_interestRegions.bed"),
          col_names = FALSE, quote = "none")

# Split into out-of-tads and high-confidence tads
wellDefinedValleys <- interestRegions %>% filter(V4 == "Valley_1")
wellDefinedMountains <- interestRegions %>% filter(V4 == "Majority_3")

write_tsv(wellDefinedValleys,
         paste0(outdir,"/wd_valleys.bed"),
         col_names = FALSE, quote = "none")

write_tsv(wellDefinedMountains,
         paste0(outdir,"/wd_mountains.bed"),
         col_names = FALSE, quote = "none")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualize statistics ----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newOrography$V4 <- factor(newOrography$V4, 
                       levels = c("Valley_1", "Hill_3", "Mountain_2"),
                       labels = c("Out of TAD", "Fuzzy", "Majority vote TAD"))


countPlot<- ggplot(newOrography, aes(x = V4 , fill = V4)) +
                geom_bar(position = "dodge", color = "black") +
                labs(x = "Region class",
                     y = "Number of regions") +
                scale_fill_manual(values = c("#039399ff","#e25baeff","#ffbd59ff")) +
                coord_flip() + 
                theme_pubclean() +
                theme(legend.position = "none",
                      axis.title.x = element_text(size = 8))

ggsave(countPlot, 
       filename = paste0(outdir, "/region_type_count.pdf"),
       device = "pdf",
       width = 6.93,
       height = 2.07)

df <- newOrography %>%
                mutate(size = V3 - V2,
                       mTypeLabel = ifelse(mType == "Between hills",
                                           "B",
                                           ifelse(mType == "Valley|Hill",
                                                  V4, "H"))
                       ) %>% 
        group_by(V4, mTypeLabel) %>% 
        summarise(s = sum(size))

#df
total <- sum(df$s)

sizePlot <- ggplot(df, aes(x = V4, y = s*100/total, fill = mTypeLabel)) +
                geom_bar(stat = "identity", position = "stack", color = "black") +
                scale_fill_manual(values = c("#039399ff","#e25baeff",
                                             "#fff1ddff","#ffddabff"),
                                  labels = c("Out of TAD", "Fuzzy",
                                             "Between fuzzy regions", "High confidence TADs")) +
                labs(x = "Region class", 
                     y = paste0("% genome in class \n (Total length ", as.character(total/1e6)," Gb)"),
                     fill = "Class refinement") + 
                coord_flip() +
                theme_pubclean() + 
                theme(legend.position = "right")
#sizePlot
ggsave(sizePlot, 
       filename = paste0(outdir, "/size_class.pdf"),
       device = "pdf",
       width = 6.93,
       height = 2.07)
