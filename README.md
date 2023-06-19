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

Run the following scripts within the repository folder:
```
./setupTriforce.sh
./runJuicer.sh <pathToYourFastaFiles>
./runHiCExplorer.sh <pathToYourFastaFiles>
./triforce.sh
```
