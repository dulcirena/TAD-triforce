#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Arguments  --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  args <- commandArgs(trailingOnly = TRUE)
  #print("hello wordl!")

  # Working directories
  wkpath <- args[1]
  outdir <- args[2]
  setwd(wkpath)

  pathA <- args[3]
  pathH <- args[4]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameters --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #  
 # wkpath <- "/homes/bierdepot/dulce/Documents/MasterDul/master-evo-in-3D/TRIFORCE/TAD-triforce/Replicas/Dmel/5kb"
 # setwd(wkpath)
 # outdir <- "/homes/bierdepot/dulce/Documents/MasterDul/master-evo-in-3D/TRIFORCE/TAD-triforce/Replicas/Dmel/5kb/out"
 # pathA <- "arrowhead_tads"
 # pathH <- "dmel_Merge_TAD_5kb_domains.bed" # *domains.bed

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Libraries ---------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(dplyr)
library(tidyverse)
library(plotly)
library(ggpubr)
library(scico)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Auxiliary functions------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modify_voters <- function(value){
  return(1*value[1] + 2*value[2] + 4*value[3])
} 

as_bed_MV <- function(df){
  
  color <- "212,147,71"
  df <- df %>% 
    select(!maxVotes) %>%
    mutate(orography = "Majority",
           x = "0",
           strand = ".",
           s2 = df$start,
           e2 = df$end,
           rgb = color)
  
  return(df)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading and wrangling ----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Upload data
A <- read_tsv(pathA) %>%
          select(chr1, x1, x2) %>% 
          rename("chr" = chr1,"start" = x1, "end" = x2) 

# H: HiCExplorer
H <- read.csv(pathH,
              header = FALSE, sep = "\t") %>%
              select(V1, V2, V3) %>%
              rename("chr" = V1,"start" = V2, "end" = V3)            

pathS <- paste0(outdir,"/mountains_raw.bed")
# S: Structure change 
S <-  read.csv(pathS,
               header = FALSE, sep = "\t") %>%
               select(V1, V2, V3) %>%
               rename("chr" = V1,"start" = V2, "end" = V3)            

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Start elections --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create table:
# It should have this format: 
# <Voter> <Side> <Position> <Chr?> 

electionsOK <- bind_rows(A %>% mutate(voter = "A"),
                       H %>% mutate(voter = "H"),
                       S %>% mutate(voter = "S")) 

chrs <- unique(electionsOK$chr)

allMajority <- data.frame()

for (k in chrs){
    curChr <- k 
    print(curChr)
    
    elections <- electionsOK %>% 
                        filter(chr == curChr) %>%
                        pivot_longer(!c(chr, voter), values_to = "position") %>%
                        arrange(position,name) #%>%   # ends will appear first in case of
                                                     # duplicated positions
                        #slice(1:15)
    
    # length:
    #sizeTest <- elections$position[15] - elections$position[1]
    print("Data wranglin: Done!")
    
    # Let's vote!
    n <- dim(elections)[1]
    
    votesCount <- 0
    votes <- c()
    
    whoVotedVect <- c()
    voterValues <- list("A" = 1, "H" = 2, "S" = 3)
    votersState <- c(FALSE, FALSE, FALSE)
    print("Votation: starting...")
    
    for (i in 1:n){
      voter <- elections$voter[i]
      
      if(elections$name[i] == "start"){
        votesCount <- votesCount + 1    
        votersState[voterValues[[voter]]] <- TRUE
        
      }else{
        votesCount <- votesCount - 1
        votersState[voterValues[[voter]]] <- FALSE
      }
      
      votes <- c(votes, votesCount)
      whoVotedVect <- c(whoVotedVect, modify_voters(votersState))
    }
    
    elections$votes <- votes
    elections$who <- whoVotedVect
    
    print("Votation: done!")
    # For each repeated position, get max(votes)
    print("Counting votes")
    finalCount <- elections %>% 
      group_by(position) %>%
      summarise(max(votes), max(who))
    
    # Create ranges
    print("Creating ranges..")
    notSingles <- c(1,2,4)
    
    finalCount$end <- c(finalCount$position[-1], 0) 
    
    finalCount <- finalCount %>%
      mutate(start = position,
             votes = `max(votes)`,
             who = `max(who)`) %>%
      select(start, end, votes, who) %>% 
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      mutate(won = (votes>=2 & !(who %in% notSingles))) 
#      filter(votes >= 2, !(who %in% notSingles))
      # do not merge ranges!
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # 
    
    # Colapse ranges with >=2 votes
    # Each range should be the first start with won==TRUE and the
    # last end with won==TRUE, for each set of consecutive won==TRUE
    
    print("Colapsing ranges...")
    n <- dim(finalCount)[1]
    inRange <- FALSE
    majorityRanges <- data.frame()

    for (i in 1:n){
      if(finalCount$won[i] == TRUE){          # Starting a new range to save
        if(!inRange){                         # Is it a new one (!inRange)
          # or are we inside a range (inRange)?
          inRange <- TRUE                     # Change flag status
          start <- finalCount$start[i]        # Initialize new range
          end <- finalCount$end[i]
          maxVotes <- finalCount$votes[i]
        }else{
          end <- finalCount$end[i]                # Update end and maxVotes for the
          if(finalCount$votes[i] > maxVotes){     # current range
            maxVotes <- finalCount$votes[i]
          }
        }
      }else{
        if(inRange){    # The new range has ended so we save it and change
          # flag value
          majorityRanges <- rbind(majorityRanges, c(start, end, maxVotes))
          inRange <- FALSE
        }
      }
    }
    
    colnames(majorityRanges) <- c("start", "end", "maxVotes")
    
    majorityRanges <- majorityRanges %>%
      mutate(chr = curChr) %>%
      select(chr, start, end, maxVotes)
    print("Writing file")
    
    write_tsv(majorityRanges,
              paste0(outdir, 
                     "/majorityRanges_",
                     curChr,".bed"))
    allMajority <- rbind(allMajority, majorityRanges)
}

# Check size:
 # allMajority %>% 
 #       mutate(size = end - start) %>% 
 #       summarise(s = sum(size))

allMajority <- allMajority %>% 
                  arrange(chr,start) %>% 
                  unique()
write_tsv(as_bed_MV(allMajority),
          paste0(outdir, 
                 "/majorityRanges.bed"),
          quote = "none", col_names = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualize election statistics -------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all <- bind_rows(A %>% mutate(voter = "A"),
                 H %>% mutate(voter = "H"),
                 S %>% mutate(voter = "S"),
                 allMajority %>%
                   select(chr, start, end) %>% 
                   mutate(voter = "M"))

all$sizes <- all$end - all$start
all$voter <- factor(all$voter, 
                    levels = c("M", "S", "H", "A"),
                    labels =  c("Majority vote",
                                "Structural change",
                                "HiCExplorer",
                                "Arrowhead"))

thresholdBox <-
  ggplot(all, aes(sizes/1e3, x = voter, fill = voter)) +
    geom_boxplot(color = "black", outlier.shape = NA) +
    labs(y = "Length (kb)",
         x = paste0("MV: n = ", dim(allMajority)[1], "; ",
                    "SC: n = ", dim(S)[1], "; ",
                    "HiCExp: n = ", dim(H)[1], "; ",
                    "Arrowh: n = ", dim(A)[1], ""),
         title = "Domain sizes per tool")+
    scale_fill_scico_d(palette = "batlow", alpha = 0.5, 
                       begin = 1, end = 0) +
    theme_minimal() + 
    theme(legend.position = "none",
          axis.title.x = element_text(size = 7))

htmlwidgets::saveWidget(ggplotly(thresholdBox)
                        , file = paste0(outdir,
                                        "/domainSizes.html"))

#all %>% filter(sizes >= 1500*1e3)

# save svg for docs !

# NEXT:----------------------------
#       finalOrography.sh         |
#       selectOrography.sh        |
# **********************************
