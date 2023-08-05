#/usr/bin/bash


# ARGUMENTS -----------------------------------------------------
# Principal (root)  working directory
wkdir=$1
# Folder with raw sequencing files
datadir=$2

# Arguments for the structural analysis
tssFile=$3
res=$4
projectName=$5
discardChr=$6
arrowFile=$7
hicfFile=$8

# Arguments for majority vote:

# Check ---------------------------------------------------------
echo "Your parameters:"
echo "wkdir:" $wkdir
echo "datadir:" $datadir
echo "resolution:" $res
echo "h value": $hvalue

# Avoid relative/absolut paths issues
curdir=$PWD
echo "Running at:" $curdir
cd $wkdir
wkdir=$PWD
cd $curdir
cd $datadir
datadir=$PWD
cd $curdir

echo "wkdir:" $wkdir
echo "datadir:" $datadir
echo

#----------------------------------------------------------------
echo " ▲   TRIFORCE -------------------------------------- "
echo "▲ ▲  Step 0: Dir setup"                            |
echo "----------------------------------------------------"

#./0_set_dir_tree.sh $wkdir $datadir

echo "Outputs can be found at:"  ${wkdir}/out
echo "Input data taken from:" ${wkdir}/data

echo 
echo " ▲   TRIFORCE -------------------------------------- "
echo "▲ ▲  Step 1: Structural change analysis"
echo "----------------------------------------------------"

#Rscript 1_structure_analysis.R $wkdir \
#				${wkdir}/out \
#				$tssFile \
#				$res

echo
echo " ▲   TRIFORCE -------------------------------------- "
echo "▲ ▲  Step 2: Initial categorization of regions"
echo "----------------------------------------------------"

#Rscript 2_merge_breakpoints.R $wkdir ${wkdir}/out \
#				$projectName $discardChr

echo " ▲   TRIFORCE -------------------------------------- "
echo "▲ ▲  Step 3: Majority vote"
echo "----------------------------------------------------"

#Rscript 3_majority_vote.R $wkdir ${wkdir}/out \
#				$arrowFile $hicfFile

echo " ▲   TRIFORCE -------------------------------------- "
echo "▲ ▲  Step 4: High confidence TADs"
echo "----------------------------------------------------"

echo "▲---- Re-classify acccording to majority voting"

#./4_finalOrography.sh ${wkdir}

# Select high confidence TADs, i.e., those bounded by at least
# one out of TAD region:

echo "▲---- Select high-confidence TADs"

Rscript 5_selectOrography.R $wkdir \
			    ${wkdir}/out \
                            $projectName

echo "🪰 Hey! Listen!"
echo " ... If TRIFORCE was useful for your research please cite:"
