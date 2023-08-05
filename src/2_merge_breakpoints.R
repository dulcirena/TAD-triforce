#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Arguments  --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
args <- commandArgs(trailingOnly = TRUE)
#print("hello wordl!")

# Working directories
wkpath <- args[1]
outdir <- args[2]
setwd(wkpath)

# Path to the TAD separation score (tad_score.bm) file
species <- args[3]
discardChr <- args[4]

print("Running step 2/5: MERGE BREAKPOINTS")
print("Parameters:")
cat(args, sep = "\n")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameters --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# wkpath <- "/homes/bierdepot/dulce/Documents/MasterDul/master-evo-in-3D/TRIFORCE/TAD-triforce/Replicas/Dmel/5kb"
# setwd(wkpath)
# outdir <- "/homes/bierdepot/dulce/Documents/MasterDul/master-evo-in-3D/TRIFORCE/TAD-triforce/Replicas/Dmel/5kb/out"
# species <- "Dmel"
# discardChr <- "NC_004353.4"
            # Dsim: "NC_029796.1"

# # for Dsim:
# # chrRefSeq <- c("NT_479533.1","NT_479534.1","NT_479535.1",
# #                "NT_479536.1","NC_029796.1","NC_029795.1") 
# # for Dmel:
# chrRefSeq <- c("NT_033779.5", "NT_033778.4"
#              ,"NT_037436.4" ,"NT_033777.3"
#              ,"NC_004353.4" ,"NC_004354.4")
# #chrMuller <- c("2L", "2R", "3L", "3R","4","X")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Libraries ---------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(plotly)
library(ggpubr)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Auxiliary functions------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
get_chr_name_refseq <- function(x){
  xSplit <- unlist(str_split(unlist(str_split(x,"myData_"))[2],"[.]"))
  chrName <- paste0(xSplit[1], ".", xSplit[2])
  
  return(chrName)
}

get_chr_name <- function(x){
  x <- outDFiles[j]
  xSplit <- tail(unlist(strsplit(x, "/")), n=1) %>%
                str_replace("myData_", "") %>%
                str_replace(".tsv", "")
  
  return(xSplit[1])
  # xSplit <- unlist(str_split(unlist(str_split(x,"myData_"))[2],"[.]"))
  # #chrName <- paste0(xSplit[1], ".", xSplit[2])
  # 
}

make_range <- function(df, chr){
  
  start <- df$position[1]
  n <- dim(df)[1]
  end <- df$positionEnd[n]
  avg <- mean(df$metric)
  
  return(cbind(chr, start, end, avg))
}

as_bed <- function(df, type){
  
  if(type == "Valley"){
    color <- "151,183,146"
  }else{
    ifelse(type == "Mountain",
           color <- "25,13,51",
           color <- "135,95,121")
  }
  
  df <- df %>% 
    select(!avg) %>%
    mutate(orography = type,
           x = "0",
           strand = ".",
           s2 = df$start,
           e2 = df$end,
           rgb = color)
  
  return(df)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading  -----------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outBFiles <- list.files(outdir, pattern = "mybreakpoints_", full.names = TRUE)
outDFiles <- list.files(outdir, pattern = "myData_", full.names = TRUE)
f <- length(outBFiles)

# Obtain the names of the chromosomes: 
# chrRefSeq <- unlist(strsplit(outBFiles, "/"))
# n <- length(chrRefSeq)/f
# chrRefSeq <- str_replace(chrRefSeq[seq(1:(f*n)) %% n == 0], "mybreakpoints_", "")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create genomic ranges from breakpoints ----------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

structureRanges <- data.frame()
chrRefSeq <- c()

for (j in 1:f){
  breakpoints <- read.csv(outBFiles[j])
  data <- read_tsv(outDFiles[j], show_col_types = FALSE)
  chr <- get_chr_name(outDFiles[j])
  if(!(chr %in% chrRefSeq))
    chrRefSeq <- c(chrRefSeq, chr)
  
  n <- dim(breakpoints)[1]
  
  for (i in 1:(n-1)){
    start <- breakpoints$x[i] + 1
    end <- breakpoints$x[i+1]
    
    newdf <- data %>% slice(start:end)
    structureRanges <- rbind(structureRanges, 
                             make_range(newdf, chr))
  }
}

colnames(structureRanges) <- c("chr", "start", "end", "avg")

structureRanges$start <- as.numeric(structureRanges$start)
structureRanges$end <- as.numeric(structureRanges$end)
structureRanges$avg <- as.numeric(structureRanges$avg)

SRlength <- sum(structureRanges$end - structureRanges$start)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get thresholds for "orography" classification ---------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thresholds <- structureRanges %>% 
                  filter(chr != discardChr) %>% 
                  summarize(dep = quantile(avg, 0.25),
                            mount = quantile(avg, 0.50))

tDepresions <- thresholds$dep
tMountain <- thresholds$mount

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Classify regions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mountains_raw <- structureRanges %>% 
                    filter(avg >= tMountain) 
depresions_raw <- structureRanges %>% 
                    filter(avg <= tDepresions)
hills_raw <- structureRanges %>%
                    filter(avg > tDepresions, 
                           avg < tMountain)
# Write output
write_tsv(as_bed(mountains_raw, "Mountain"), 
          paste0(outdir,"/mountains_raw.bed"),
          col_names = FALSE, quote = "none")

write_tsv(as_bed(depresions_raw, "Valley"),
          paste0(outdir, "/valleys_raw.bed"),
          col_names = FALSE, quote = "none")

write_tsv(as_bed(hills_raw, "Hills"),
          paste0(outdir,"/hills_raw.bed"),
          col_names = FALSE, quote = "none")

# Check size
# bind_rows(mountains_raw, hills_raw, depresions_raw) %>% 
#       mutate(size = end - start) %>% 
#       summarise(s = sum(size))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization of avg. contacts per region/chr   -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
structureRanges$chr <- as.factor(structureRanges$chr)
structureRanges$chr <- factor(structureRanges$chr, 
                              levels = chrRefSeq,
                              labels = chrRefSeq)

thresholdsBox <- ggplotly(
                    
  ggplot(structureRanges, aes(x = "All", y = avg)) +
      geom_boxplot(fill = "#A5A1EB", color = "black") +
      labs(x = "", y = "") + 
      coord_flip() +
      theme_pubclean() +
      labs(title = paste0("Average TAD sep. score per SC-region \nin ",
                          species)) +    
      geom_boxplot(data = structureRanges %>% filter(chr !=  "4"),
                   aes(x = "All except chr4", y = avg), 
                   fill = "#D2CEF6", color = "black")  +
      geom_boxplot(aes(x = chr, y = avg), 
                  fill = "#184E60", color = "black") 
                        
)

htmlwidgets::saveWidget(thresholdsBox
                        , file = paste0(outdir,
                                        "/avgSCregion_boxplot.html"))


#NEXT:---------------------------
#       majority_vote.R         %
#       finalOrography.sh       %
#                               %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
