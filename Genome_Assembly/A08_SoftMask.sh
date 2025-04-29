#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=16G
#$ -l rmem=16G
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 2
#$ -v OMP_NUM_THREADS=2
#$ -j y

#####################################################################################
#	Script Name: 	SoftMask.sh
#	Description: 	Soft mask genome using EDTA and RepeatMasker libraries
#	Author:		LPereiraG
#	Last updated:	28/06/2022, LPG
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and files
wd=/mnt/fastdata/bo1lpg/Aristida-v3
input=${wd}/A03_results/Aristida_nuclear.asm.bp.p_ctg.fa
out=${wd}/A07_results/Aristida_nuclear.asm.bp.p_ctg.fa.out

#### Parameters
label="Aristida_nuclear"

#### Step 1:
cd ${wd}
mkdir A08_results
cd A08_results
sed -e '1,3d' ${out} | awk -v OFS='\t' '{print $5, $6-1, $7}' | sort -k1,1 -k2,2 -V > ${label}.bed
bedtools maskfasta -soft -fi ${input} -bed ${label}.bed -fo ${label}-softMasked.fa
