#!/bin/bash

#SBATCH -J mpi_ofjob
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -A birthright
#SBATCH -p testing
#SBATCH -o ./test-output.txt
#SBATCH -e ./test-error.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yangf@ornl.gov

module purge
module load PE-gnu
source $HOME/OpenFOAM/OpenFOAM-v1906/etc/bashrc
./Allclean
./Allrun
foamToVTK

