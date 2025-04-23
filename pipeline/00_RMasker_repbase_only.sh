#!/usr/bin/bash
#SBATCH -p short -N 1 -n 32 --mem 24gb --out logs/rmasker_repbase.log


module load RepeatMasker/4-0-7
module unload ncbi-rmblast
module load ncbi-rmblast/2.9.0-p2
module unload miniconda2
module load miniconda3
module load mcclintock
source activate mcclintock
GENOME=Coccidioides_posadasii_Silveira.scaffolds.fa
BASE=$(basename $GENOME .fa)
RepeatMasker -pa 32 -s -e ncbi -lib lib/Cpos_C735_refTEs.fa.classified genome/$GENOME >& logs/$BASE.RM.out

perl scripts/repeatmasker_to_gff3.pl genome/$GENOME.out > genome/$BASE.repeat_masker.gff3
