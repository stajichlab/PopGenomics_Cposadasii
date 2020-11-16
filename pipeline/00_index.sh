#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 --mem 3gb --out logs/index.log
module load samtools/1.9
module load bwa/0.7.17
if [ -f config.txt ]; then
	source config.txt
fi

# Genome was manually installed here
# From Local project work see 
# /bigdata/stajichlab/stajichcollab/Tgen/Silveira_pac/annotation/Pilon/annotate/Coccidioides_posadasii_Silveira/annotate_results/
# THIS IS EXAMPLE CODE FOR HOW TO DOWNLOAD DIRECT FROM FUNGIDB
FASTAFILE=$REFGENOME
if [[ ! -f $FASTAFILE.fai || $FASTAFILE -nt $FASTAFILE.fai ]]; then
	samtools faidx $FASTAFILE
fi
if [[ ! -f $FASTAFILE.bwt || $FASTAFILE -nt $FASTAFILE.bwt ]]; then
	bwa index $FASTAFILE
fi

DICT=genome/$(basename $FASTAFILE .fa)".dict"
echo "$DICT"
if [[ ! -f $DICT || $FASTAFILE -nt $DICT ]]; then
	samtools dict $FASTAFILE > $DICT
	ln -s $(basename $DICT) $FASTAFILE.dict 
fi
