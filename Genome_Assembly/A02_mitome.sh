#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=12G
#$ -l rmem=12G
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 6
#$ -v OMP_NUM_THREADS=6
#$ -j y

#####################################################################################
#       Script Name:    mitome.sh
#       Description:    Fish mitochondria HiFi reads and assemble them using MitoHiFi
#       Author:         LPereiraG
#       Last updated:   03/08/2022
#####################################################################################

#### Directories and files
wd=/mnt/fastdata/bo1lpg/Aristida-v2
mitochondria=${wd}/data
reads=/shared/dunning_lab/User/bo1lpg-backup/Aristida/data/Aristida_reads.fa

#### Parameters

#### Step 1: Get a similar mitome
mkdir ${wd}/A02_results-v2
cd ${wd}/A02_results-v2
apptainer exec --bind ${wd} /usr/local/packages/singularity/images/mitohifi/mitohifi.sif findMitoReference.py \
  --species "Aristida adscensionis" --outfolder ${wd}/A02_results-v2 -t mitochondrion

#### Step 2: Run the main script to assemble mitome
cp ${reads} ${mitochondria}
apptainer exec --bind ${wd} /usr/local/packages/singularity/images/mitohifi/mitohifi.sif mitohifi.py \
  -r ${mitochondria}/Aristida_reads.fa -f ${wd}/A02_results-v2/*.fasta -g ${wd}/A02_results-v2/*.gb \
  -t 6 -a plant

#### Step 3: Clean hifiasm results to keep only unique sequences
cd reads_mapping_and_assembly
module load apps/python/anaconda2-4.2.0
source activate /shared/dunning_lab/Shared/conda_env/circlator
circlator clean hifiasm.contigs.fasta hifiasm.contigs.clean
