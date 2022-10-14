#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                           (THE) TRIFORCE                                #
# autor: Dulce I. Valdivia                                                #
# contact: dulce.valdivia@cinvestav.mx                                    #
# github: dulcirena                                                       #
#                                                                         #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Usage from command line -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Arguments  --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# wkpath <- "C:/Users/dival/Documents/Triforce/test"
# setwd(wkpath)
# outdir <- "C:/Users/dival/Documents/Triforce/test/out/"
# pathTSS <- "C:/Users/dival/Documents/Triforce/test/data/eye/Dmel_eye_5kb_thres0.01_delta0.01_fdr_tad_score.bm"
# hvalue <- 20 # set to 10 for 10kb resolution or to 20 for 5kb res.

args <- commandArgs(trailingOnly = TRUE)
#print("hello wordl!")


# Working directories
wkpath <- args[1]
outdir <- args[2]
setwd(wkpath)

# Path to the TAD separation score (tad_score.bm) file
pathTSS <- args[3]

# The parameter h is the minimum of bins that should
# represent one partition. Since TADs are approx. ~100 kb,
# We adjust h = 100 / r, where r is the resolution of the
# matrix. E.g. If the TADs where computed in a matrix with
# resolution r = 5 kb, then h = 20; or if they were called
# in a r = 10 kb resolution, then h = 10.

hvalue <- floor(100/as.numeric(args[4]))

print("Running step 1/5: STRUCTURAL CHANGE")
print("Parameters:")
cat(args, sep = "\n")
cat(paste("h value:", hvalue))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Libraries ---------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(dplyr)
library(strucchangeRcpp)
library(plotly)
library(lubridate)
library(tictoc)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data loading  -----------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input: *_tad_score.bm
dfTADSep <- read.csv(pathTSS,
                     sep = "\t", header = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Formatting --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dfTADSep <- dfTADSep[-1,]
dfTADSep$metric <- apply(dfTADSep[,c(4:10)],1, sum)/7

tic()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Change point dectection -------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# chrList <- c("NT_033779.5", "NT_033778.4"
#              ,"NT_037436.4" ,"NT_033777.3"
#              ,"NC_004353.4" ,"NC_004354.4")

chrList <- unique(dfTADSep$V1)

for (chr in chrList){
  #chr <- chrList[1]
  print(paste(chr, " ..."))

    myData <- dfTADSep %>%
                  filter(V1 == chr) %>%
                  mutate(chr = V1, position = V2, positionEnd = V3) %>%
                  select(chr, position, positionEnd, metric)


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Change point dectection -------------------------------------------------
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    print(">>> Computing breaks:")
    # If we want average values over periods
    myData_breakpoints <- breakpoints(metric ~ 1
                                      , h = hvalue
                                      , data = myData
    )
    print(">>> DONE!")

    myBreakpoins <- c(0, myData_breakpoints$breakpoints, dim(myData)[1])

    # Save info:
    #models[[iplot]] <- myData_breakpoints

    write.csv(myBreakpoins
              ,paste0(outdir
                      ,"/mybreakpoints_"
                      ,chr)
              , quote = FALSE, row.names = FALSE)
    write_tsv(myData
              , paste0(outdir
                       ,"/myData_"
                       , chr
                       , ".tsv"),
              quote = "none")

    print(">>> Main output saved!")

    # Identify the number of break points
    sss <- summary(myData_breakpoints)$RSS[2,]

    optimum_num_breaks <- which.min(sss) - 1

    myData$avg <- as.numeric(fitted(myData_breakpoints
                                    , breaks = optimum_num_breaks)
    )

    xxx <- confint(myData_breakpoints
                   , breaks = optimum_num_breaks
    )
    xxx$confint
    
    # Inpute C.I. NA values to the breakpoint (CI 50%)
    idx_na <- which(is.na(xxx$confint[,1]) == TRUE)
    # Print idx_na
    xxx$confint[idx_na,1] <- xxx$confint[idx_na,2]
    xxx$confint[idx_na,3] <- xxx$confint[idx_na,2]
    
    write_tsv(as.data.frame(xxx$confint)
              , paste0(outdir
                       ,"/confidenceInterval"
                       , chr
                       , ".tsv"),
              quote = "none")

    print(">>> CIs saved!")

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Visualizations ----------------------------------------------------------
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #xxx$confint2 <- as.tibble(xxx$confint) %>% drop_na()

    YAXIS.TITLE_ <- "TAD sep. score"

    # Prepare the output plot
    # Decoration for the time series
    # scale_factor <- 0.50
    # Y0 <- mean(myData$metric, na.rm = TRUE)*(1-scale_factor)
    # Y1 <- mean(myData$metric, na.rm = TRUE)*(1+scale_factor)

    SHAPES_ <- list()

    for(j in 1:optimum_num_breaks){
        SHAPES_[[j]] <-  list(type = "rect"
                              , fillcolor = "green"
                              , line = list(color = "green")
                              , opacity = 0.2
                              , x0 = myData$position[xxx$confint[j,2]]
                              , x1 = myData$position[xxx$confint[j,3]]
                              , xref = "x"
                              , y0 = -1
                              , y1 = 2
                              , yref = "y"
        )
    }

    # Main plot
    output <- plot_ly(data = myData
                      , x = ~position
                      , y = ~metric
                      , type = "scatter"
                      , mode = "lines"
                      , name = "Observed"
                      , line = list(width = 4)
    ) %>%
        add_trace(y = ~avg
                  , line = list(dash = "dash", color = "black", width = 2)
                  , name = "Trend"
        )  %>%
        layout(shapes = SHAPES_
               , yaxis = list(title = YAXIS.TITLE_
                              , range = c(-1, 2)
               )
               , xaxis = list(title = "")
        )


    htmlwidgets::saveWidget(output
                            , file = paste0(outdir,
                                            "/structure_",chr,".html"))
    print(">>> Visualization saved!")
}

print("STRUCTURAL CHANGE ANALYSIS DONE --------------")
toc()


#NEXT:---------------------------
#       merge_breakpoints.R       |
#       majority_vote.R           |
#       finalOrography.sh         |
#       makeTracks.sh             |
#---------------------------------