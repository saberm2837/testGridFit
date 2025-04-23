#!/bin/bash

startJobId=$1
endJobId=$2
for (( jobId=$startJobId; jobId<=$endJobId; jobId++ ))
do
	scancel $jobId

done

echo "Total "$((endJobId-startJobId+1))" jobs cancelled."
