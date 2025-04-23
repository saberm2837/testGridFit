#!/bin/bash
export PATH=$PATH:$HOME/julia/bin

job_count=0
sub_n=$1
points=$2

for rep_n in {1..300}
do
	sbatch testGridFit_V6.slurm $sub_n $points $rep_n
	((job_count++))
done

echo "Total "$job_count" jobs submitted."
