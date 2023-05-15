curdir=$PWD # Save working dir
outdir=$1/out/ # Set output dir
cd $outdir

VALLEY=valleys_raw.bed
HILL=hills_raw.bed
MV=majorityRanges.bed

# REFINE VALLEYS AND HILLS -----------------------------------------------------
	bedtools subtract -a $VALLEY -b $MV > final_valleys.bed
        bedtools subtract -a $HILL -b $MV > final_hills.bed

#-------------------------------------------------------------------------------

# MAKE FINAL OROGRAPHY.BED -----------------------------------------------------------
	awk '{$4="Valley_1"; $5="0"; $6="."; $7=$2; $8=$3; $9="151,183,146"}1' \
		final_valleys.bed | sed 's/ /\t/g' > valleys_out.bed

	awk '{$4="Hill_3"; $5="0"; $6="."; $7=$2; $8=$3; $9="135,95,121"}1' \
		final_hills.bed | sed 's/ /\t/g' > hills_out.bed

	awk '{$4="Mountain_2"; $5="0"; $6="."; $7=$2; $8=$3; $9="25,13,51"}1' \
		$MV | sed 's/ /\t/g' > mountains_out.bed

	cat valleys_out.bed hills_out.bed mountains_out.bed \
		 | sort -k 1,1 -k2,2 > orography_out.bed

#       cat valleys_out.bed hills_out.bed mountains_out.bed \
#                 | grep -v "chr" | sort -k 1,1 -k2,2 > orography_out.bed

#---- Sizes:
statsFile=classification_stats.txt

echo "FINAL SIZES" > $statsFile
echo "Out of TAD regions (bp):" >> $statsFile
awk  '{$4=$3-$2}1' valleys_out.bed | awk '{sum+=$4} END {print sum}' \
	>> $statsFile
#42244000

echo "Majority vote regions (bp):" >> $statsFile
awk  '{$4=$3-$2}1' mountains_out.bed | awk '{sum+=$4} END {print sum}' \
	>> $statsFile
#80088000

echo "Fuzzy regions (bp):" >> $statsFile
awk  '{$4=$3-$2}1' hills_out.bed | awk '{sum+=$4} END {print sum}' \
	>> $statsFile

cd $curdir
