#!/bin/bash

#SBATCH --partition=test
#SBATCH --job-name=generate_frequencies_eigenvectors
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=00:20:00
#SBATCH --account=pn73da
#SBATCH --mail-user=jonas.grandel@bam.de
#SBATCH --mail-type=ALL
#SBATCH --output=/hppfs/scratch/02/di35jom/PCM_calculations/06/b2/4e/06b24efe-1443-43c8-90c0-554a58725a11_1/queue.out
#SBATCH --error=/hppfs/scratch/02/di35jom/PCM_calculations/06/b2/4e/06b24efe-1443-43c8-90c0-554a58725a11_1/queue.err
#SBATCH --exclusive=exclusive
#SBATCH --get-user-env
cd /hppfs/scratch/02/di35jom/PCM_calculations/06/b2/4e/06b24efe-1443-43c8-90c0-554a58725a11_1
module load slurm_setup
module load vasp/6.2
source activate jobflow_remote

jf execution run /hppfs/scratch/02/di35jom/PCM_calculations/06/b2/4e/06b24efe-1443-43c8-90c0-554a58725a11_1