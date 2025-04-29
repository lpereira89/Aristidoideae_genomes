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

#### Step 1: run RepeatMasker
cd ${wd}
module load apps/python/conda
source activate /shared/dunning_lab/Shared/conda_env/repeatmasker
RepeatMasker -pa 6 -s -lib A06_results/Aristida_nuclear.asm.bp.p_ctg.fa.mod.EDTA.TElib.fa -dir A07_results -e ncbi ${input}
