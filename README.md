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
*From fastq raw data to topologically associated domains annotation*. 

To run the complete workflow for TAD annotation:

```
./triforce.sh  <WORK_DIR> \
		<FILES_DIR> <TAD_SEP_SCORE> \
		<RESOLUTION_KB> <PROJECT_NAME> <DISSMISS_CHR> \
		<FILE_ARROWHEAD_TADS> <HICEXPLORER_TADS>    
```

Description:

- **WORK_DIR**: Directory where the out/ directory will be created.
- **FILED_DIR**: Directory where the input files are.
- **TAD_SEP_SCORE**: TAD separation score file computed with HiCExplorer. Should end with tad_score.bm
- **RESOLUTION_KB**: Resolution of the matrix in kb.
- **PROJECT_NAME**: ID for the project (do not use spaces)
- **DISSMISS_CHR**: The name of the chromosome you want to dissmiss during the analysis.
- **FILE_ARROWHEAD**: Arrowhead's TAD calling file (e.g. 10000_blocks)
- **HICEXPLORER_ARROWHEAD**: HiCExplorer's TAD calling file. Should end with domains.bed
