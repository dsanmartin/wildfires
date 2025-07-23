#!/bin/bash
#---------------Script SBATCH - NLHPC ----------------
#SBATCH -J wildfire
#SBATCH -p general
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=dsanmartinreyes@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 10-24:00:00
#SBATCH -o wildfire_%j.out
#SBATCH -e wildfire_%j.err

#-----------------Toolchain---------------------------
# ----------------Modulos----------------------------
ml GCC/13.2.0 
ml FFTW/3.3.10 
# ----------------Comando--------------------------
cd /home/dsanmartin/fire/wildfires
./bin/wildfires data/input/debug.txt