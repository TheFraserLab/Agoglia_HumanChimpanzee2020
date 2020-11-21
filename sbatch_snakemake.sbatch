#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --error=snake.err
#SBATCH --output=snake.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH -p hns,normal,hbfraser

module load conda
source activate fraserconda
srun snakemake --unlock
srun snakemake --jobs 100 --cluster-config cluster_config.json --rerun-incomplete --keep-going --allowed-rules count_ase --use-conda --conda-prefix /home/groups/hbfraser/modules/packages/conda/4.6.14/envs/ --cluster "sbatch --mem={cluster.memory} --partition={cluster.partition} --nodes=1 --ntasks-per-node=1 --cpus-per-task={cluster.cpus} --time={cluster.time}" 
