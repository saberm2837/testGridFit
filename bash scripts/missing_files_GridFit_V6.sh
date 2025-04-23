#!/bin/bash

# Initialize parameters
sub_n=$1
points=$2

count=0
for run in {1..300}
do
	if ! [ -f "GridFit/sub"$sub_n"_pt"$points"_run"$run"_GridFitOptCoords.csv" ]; then
		echo "Missing Run "$run
		((count++))
	fi
done
echo "Total "$count" files missing for sub"$sub_n" with "$points" points"

if (($count == 0)); then
	echo "DONE!"
fi
