#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=16G
#$ -l rmem=16G
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 6
#$ -v OMP_NUM_THREADS=6
#$ -j y

#####################################################################################
#	Script Name: 	EDTA.sh
#	Description: 	Annotate TEs with EDTA in a reference genome
#	Author:		LPereiraG
#	Last updated:	28/06/2022, LPG
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and files
wd=/mnt/fastdata/bo1lpg/Aristida-v3
input=${wd}/A03_results/Aristida_nuclear.asm.bp.p_ctg.fa

#### Parameters

#### Step 1: run EDTA
cd ${wd}
mkdir A06_results
cd A06_results
source activate EDTA
EDTA.pl --genome ${input} -t 6 --anno 1
