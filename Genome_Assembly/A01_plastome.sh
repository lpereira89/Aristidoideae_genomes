#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=12G
#$ -l rmem=12G
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 4
#$ -v OMP_NUM_THREADS=4
#$ -j y

#####################################################################################
#       Script Name:    hifiasm_assembly.sh
#       Description:    Assembly HiFi reads
#       Author:         LPereiraG
#       Last updated:   22/04/2022
#####################################################################################

module load apps/python/anaconda2-4.2.0

#### Directories and files
wd=/mnt/fastdata/bo1lpg/Aristida-v3
chloroplast=${wd}/data/Aristida-adscenionis-chloroplast.fa
reads=/shared/dunning_lab/User/bo1lpg-backup/Aristida/data/Aristida_reads.fa

#### Parameters


#### Step 1: Fish chloroplast reads
cd ${wd}
mkdir A01_results
cd A01_results
source activate /shared/dunning_lab/Shared/conda_env/blasr
blasr ${reads} ${chloroplast} -m 1 --header --minAlnLength 5000 --minPctSimilarity 99 >> Blast_raw
cat Blast_raw | cut -f 1 -d " " | sort | uniq > reads_to_fish_chloroplast.txt
sed -i 's/\/0_[0-9]\{5\}$//g' reads_to_fish_chloroplast.txt
sed -i 's/\/0_[0-9]\{4\}$//g' reads_to_fish_chloroplast.txt
source /usr/local/extras/Genomics/.bashrc
seqtk subseq ${reads} reads_to_fish_chloroplast.txt >> chloroplast_reads.fa

#### Step 2: Assemble plastome
/shared/dunning_lab/Shared/programs/hifiasm-0.16.1/./hifiasm -o Aristida_chloroplast.asm -t 4 chloroplast_reads.fa
awk '/^S/{print ">"$2;print $3}' Aristida_chloroplast.asm.bp.p_ctg.gfa > Aristida_chloroplast.asm.bp.p_ctg.fa

#### Step 3: Circularize plastome
source activate /shared/dunning_lab/Shared/conda_env/circlator
circlator clean Aristida_chloroplast.asm.bp.p_ctg.fa Aristida_chloroplast_clean
