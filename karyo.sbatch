#!/bin/bash

#SBATCH --job-name=karyo_HL1-29_P23_1
#SBATCH --error=karyo_HL1-29_P23_1.err
#SBATCH --output=karyo_HL1-29_P23_1.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH -p hns,normal,hbfraser,owners

module load R/3.5.1
srun ~/bin/karyotype.R iPSC_P23_HL1-29/iPSC_P23_HL1-29_ase_by_reads_merged.txt ~/bin/positions.filled.txt iPSC_P23_HL1-29_karyo.pdf
