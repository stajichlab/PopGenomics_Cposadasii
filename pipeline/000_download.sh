#!/usr/bin/bash
#SBATCH -p batch --out  logs/sra_fetch.%A.log -n 4 --mem 2gb

module load aspera
module load sratoolkit
mkdir -p input
pushd input
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

tail -n +2 ../samples.csv | cut -d, -f2 | sed -n ${N}p | while read NAME
do
	if [ ! -f ${NAME}_1.fastq.gz ]; then
		fastq-dump  --gzip --split-e $NAME
	fi
done
