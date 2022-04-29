#! /bin/bash
#BATCH -p normal <---- the type of 'computer' you want
#SBATCH -t 15:00:00 <---- The computation time in wall-clock. In this case, it is 15hrs
#SBATCH -n 1 -c 24 <----- How many nodes and how many cores

module load 2020
module load GCC/10.3.0
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export FC=/usr/bin/gfortran
export GFORTRAN_UNBUFFERED_ALL=y

make
