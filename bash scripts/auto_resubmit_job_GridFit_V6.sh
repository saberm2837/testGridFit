#!/bin/bash
export PATH=$PATH:$HOME/julia/bin

# Initialize parameters
sub_n=$1
points=$2

count=0
for rep_n in {1..300}
do
	if ! [ -f "GridFit/sub"$sub_n"_pt"$points"_run"$rep_n"_GridFitOptCoords.csv" ]; then
		sbatch testGridFit_V6.slurm $sub_n $points $rep_n
		((count++))
	fi
done
echo "Total "$count" jobs resubmitted for sub"$sub_n" with "$points" points"
