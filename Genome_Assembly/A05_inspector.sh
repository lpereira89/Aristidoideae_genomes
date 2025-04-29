#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=8G
#$ -l rmem=8G
#$ -pe openmp 8
#$ -v OMP_NUM_THREADS=8
#$ -j y

#####################################################################################
#	Script Name: 	inspector.sh
#	Description: 	Evaluate assembly made from pacbio reads
#	Author:		LPereiraG
#	Last updated:	21/07/2022
#####################################################################################

module load apps/python/anaconda2-4.2.0
source activate /shared/dunning_lab/Shared/conda_env/inspector/

#### Parameters
sample='Aristida_nuclear.asm.bp.p_ctg.fa'

#### Directories and files
wd=/mnt/fastdata/bo1lpg/Aristida-v3
assembly=/mnt/fastdata/bo1lpg/Aristida-v3/A03_results/${sample}
reads=/mnt/fastdata/bo1lpg/Aristida-v3/A03_results/Aristida_nuclear_reads.fa

#### Step 1: evaluate misassemblies
cd ${wd}
mkdir A05_results
cd A05_results
inspector.py -c ${assembly} -r ${reads} -o ${sample} --datatype hifi -t 8
