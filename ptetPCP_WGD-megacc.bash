#!/bin/bash
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH --time=04:00:00
#SBATCH -c 1
#SBATCH -e signal_processing_job.%j.err
#SBATCH -o signal_processing_job.%j.out

module load megacc/10.2.5

in="/data/LynchLabCME/Paramecium/Tim/ptetLOPITv3/ptetPCP_WGD-aligned"
outt="/data/LynchLabCME/Paramecium/Tim/ptetLOPITv3/ptetPCP_WGD-dist"

for fasta in $in/*
do
	outfile="$(basename $fasta)"
	megacc -a ../distance_estimation_pairwise_amino_acid.mao -d $fasta -o $outt/$outfile
done
