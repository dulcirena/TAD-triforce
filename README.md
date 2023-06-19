# TRIFORCE
The Triforce consists of an automatized pipeline for processing Hi-C data with different methods. 

![image](https://user-images.githubusercontent.com/112836459/188355356-43d1fcc4-7199-4245-91e3-20517e4a4985.png)

Developed by: Dulce I. Valdivia + Luis Delaye + Kasia Oktaba

## Installation

Clone this repository in your working environment:

```
git clone https://github.com/dulcirena/TAD-triforce.git
```

## Software Requirements
1. [Juicer](https://github.com/aidenlab/juicer)
2. [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/)
3. [R](https://cran.r-project.org/)
4. R packages:
	- tidyverse
	- dplyr
	- strucchangeRcpp
	- plotly		
	- ggpubr
	- scico
5. [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
		
## Usage:
### From fastq raw data to topologically associated domains using Arrowhead and HiCExplorer
*Work in progress* 

### Compute consensus TADs using _the_ TRIFORCE
Once you have ran Arrowhead and HiCExplorer to obtain the their corresponding TAD annotation files, it is time to use TRIFORCE to obtain a consensus set of high-confidence TADs. 

To run the complete workflow for TAD annotation:

```
./triforce.sh  <WORK_DIR> \
		<FILES_DIR> <TAD_SEP_SCORE> \
		<RESOLUTION_KB> <PROJECT_NAME> <DISSMISS_CHR> \
		<FILE_ARROWHEAD_TADS> <HICEXPLORER_TADS>    
```

*Description:*

- **WORK_DIR**: Directory where the out/ directory will be created.
- **FILED_DIR**: Directory where the input files are.
- **TAD_SEP_SCORE**: TAD separation score file computed with HiCExplorer. Should end with tad_score.bm
- **RESOLUTION_KB**: Resolution of the matrix in kb.
- **PROJECT_NAME**: ID for the project (do not use spaces)
- **DISSMISS_CHR**: The name of the chromosome you want to dissmiss during the analysis.
- **FILE_ARROWHEAD**: Arrowhead's TAD calling file (e.g. 10000_blocks)
- **HICEXPLORER_ARROWHEAD**: HiCExplorer's TAD calling file. Should end with domains.bed

See the file src/run_test.sh for a working example.

*Output:*

All the outputs are stored in the directory WORK_DIR/out/

- **structure_CHR.html**: An interactive file per chromosome (CHR) showing the breakpoints of the TAD separation score (TAD-SS) according to the structural change analysis (SCA). SCA breaks the TAD-SS in genomic regions that exhibit similar contact trends. The width of each breakpoint (lightgreen) represent its confidence interval. 
- **confidenceIntervalCHR.tsv**: The coordinates of the 5% and 95% CI for each breakpoint for each CHR. The coordinate of the 50% CI is used in the downstream analysis.
- **avgSCregion_boxplot.html**: Boxplot showing the distribution of the average TAD separation score in each SCA-region. Only the regions above the overall median are kept for downstream analysis.
- **domainSizes_files**: Distribution of TAD legth in the different steps of the analysis. Usually the TADs computed by the majority vote script produces longer TADs because it merges the regions of consecutive TADs.
- **orography_out.bed**: All classified regions. Includes a color column for visualization in HiCExplorer.
- **mountains_out.bed**: High-condifence TADs
- **hills_out.bed**: Fuzzy regions
- **valleys_out.bed**: Out of TAD regions
- **region_type_count.pdf** : A plot showing the number of regions in each class (Majority Vote, Fuzzy or Out of TAD). This is made before the refinement of majority vote areas as high-confidence TADs. 
- **size_class.pdf**: A plot showing the total length of the genome classified as High-confidence TADs, TADs between fuzzy regions, Fuzzy regions and Out of TAD regions. 



