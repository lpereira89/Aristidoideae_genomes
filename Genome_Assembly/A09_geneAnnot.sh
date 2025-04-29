#!/bin/bash

# Request resources:

#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:2
#SBATCH --mem=82G
#SBATCH -t 96:00:00 # time limit for job (format: days-hours:minutes:seconds) max 96hrs

####################################################################################################
#       Script Name:    Helixer.sh
#       Description:    Get helixer annotations using gpu
#       Author:         Noah Bourne
#       Last updated:   22/02/2024
####################################################################################################

assembly_p=/mnt/fastdata/bo1lpg/Aristida-v3/

export HELIXER_SIF=/mnt/parscratch/users/bob21lp/public/Helixer/helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif

apptainer exec ${HELIXER_SIF} Helixer.py --fasta-path ${assembly_p}/Aristida_nuclear-softMasked.fa \
--lineage land_plant --batch-size 512 --gff-output-path ${assembly_p}/Aristida_nuclear-softMasked.gff3
