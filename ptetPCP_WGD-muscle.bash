#!/bin/bash
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH --time=04:00:00
#SBATCH -c 1
#SBATCH -e signal_processing_job.%j.err
#SBATCH -o signal_processing_job.%j.out

module load muscle/3.8.1551

in="/data/LynchLabCME/Paramecium/Tim/ptetLOPITv3/ptetPCP_WGD"
outt="/data/LynchLabCME/Paramecium/Tim/ptetLOPITv3/ptetPCP_WGD-aligned"

cd $outt

for fasta in $in/*
do
	outfile="$(basename $fasta)"
	muscle -in $fasta -out $outt/$outfile
done
