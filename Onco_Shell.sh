#!/usr/bin/bash

# NEEDED FILES: 
# 1. POINT MUTATIONS
# 2. COPY NUMBER DATA
#	(EITHER THE INTERMED LIST OR THE RAW TABLE AND PROCESSING LIST)
# 3. INDEL DATA
# 4. TARGET GENES
# 5. ANNOTATIONS

# IF COPY NUMBER DATA HAS NOT BEEN PROCESSED, PROCESS IT
if [ -e "CNV_Input.txt" ]
then 
	echo "Using Copy-Number data from CNV_Input.txt"
else 
# CAREFUL!!! THIS IS HARD-CODED FOR THE NEJM DATA
	echo "Compiling Copy-Number data into file, CNV_Input.txt"
	perl Process_CNV.pl NEJM.cnv.events.tsv.txt cnv_genes.lst > CNV_Input.txt
fi
exit

# ASSUME COPY NUMBER DATA EXISTS
# CREATE R INPUT (also creates file Mut_Count.anno)
perl Make_R_Input.pl $1 CNV_Input.txt $3 $4 > Onco_Input.tsv


#################### MELANOMA SPECIFIC ###########################
# DEFINE HOTSPOT INPUT AND SUB-TYPES LIST
# CREATES Melanoma_Modfifications.tsv AND Melanoma_Subtypes.anno
perl HotSpotter.pl hotspots.lst $1 $3

# RUN SCRIPT TO MODIFY R INPUT
perl Modify_R_Input.pl Onco_Input.tsv Melanoma_Modifications.mods > temp
mv temp Onco_Input.tsv
#################### MELANOMA SPECIFIC ##############################


wait;
# CREATE ANNOTATIONS FILE
ANNOS=($(ls *.anno))
ANNOLIST=${ANNOS[*]}

perl Make_Annotations.pl $5 $ANNOLIST > testAnnotations.txt


# RUN ONCOSCRIPT TO PRODUCE ONCOPRINT AS WELL AS TABULATED DATA
Rscript OncoScript.R Onco_Input.tsv testAnnotations.txt testoncoprint $4


# RUN FISHER SCRIPTS AND PRODUCE BOXPLOTS
Rscript BoxPlots.R








