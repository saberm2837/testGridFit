#!/bin/bash
export PATH=$PATH:$HOME/julia/bin

# Initialize subject number and initial points
sub_n=1
points=6

# Initialize failed repetition numbers
rep_n=(299 300)
 
for i in ${rep_n[@]};
do
	sbatch testGridFit_V6.slurm $sub_n $points $i
done

echo "Total "${#rep_n[@]}" jobs submitted."
