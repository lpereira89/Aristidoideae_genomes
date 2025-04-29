#!/bin/bash
#$ -l h_rt=96:00:00
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 2
#$ -v OMP_NUM_THREADS=2
#$ -l mem=12G
#$ -l rmem=12G
#$ -j y

#####################################################################################
#       Script Name:    busco.sh
#       Description:    Check completeness of nuclear assembly
#       Author:         LukeTDunning
#       Last updated:   01/11/2022
#####################################################################################

source /usr/local/extras/Genomics/.bashrc
alias python=python3

i=$(expr $SGE_TASK_ID)

#### Directories and files
assemblies=/mnt/fastdata/bo1lpg/Aristida-v3/A03_results/Aristida_nuclear.asm.bp.p_ctg.fa
lineages=/shared/dunning_lab/Shared/programs/BUSCO-db/poales_odb10
tempfiles=/mnt/fastdata/bo1lpg/Aristida-v3/A04_results/temp
wd=/mnt/fastdata/bo1lpg/Aristida-v3/A04_results

#### Parameters
sample="Aristida_nuclear"

cd ${wd}
# cp -r /usr/local/extras/Genomics/apps/augustus/current/config .
# export AUGUSTUS_CONFIG_PATH=/mnt/fastdata/bo1lpg/Aristida-v3/A04_results/config

# busco version: BUSCO 3.1.0

mkdir -p ${tempfiles}-${sample}
busco --in ${assemblies} --lineage_path ${lineages}/ --out ${sample} --mode geno --cpu 2 -r \
  --tmp_path ${tempfiles}-${sample}/ --tarzip
