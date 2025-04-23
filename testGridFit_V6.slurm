#!/bin/bash
#SBATCH --job-name=testGridFit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=02:00:00
#SBATCH --account=bgross
#SBATCH --partition=normal
#SBATCH --output=%x-%j.out

cd ~/GridFit
julia -L testGridFit_V6.jl -e 'test_grid_opt('$1','$2','$3')'
