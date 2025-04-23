# GridFit Project: Mapping electrocortical grids to the brain surface

## Overview
This project evaluates the accuracy of the GridFit algorithm in estimating electrode grid locations placed inside the human brain during surgery, using fMRI data.
We systematically test GridFit, detect and resubmit any failed runs, aggregate outputs, and visualize the results.

## Project Structure

 - Raw DatafMRI skull-stripped images + ground truth locations:https://osf.io/k8v9u/
 - Generated outputs (plots, statistics):https://osf.io/gqnef/
 - Code Components
   - `testGridFit_V6.jl`: Main Julia script that runs GridFit and generates output files.
   - `testGridFit_V6.slurm`: SLURM batch script to run a single instance of testGridFit.
   - `submit_job_GridFit_V6.sh`: Shell script to submit 300 SLURM jobs for each subject/grid point configuration.
   - `missing_files_GridFit_V6.sh`: Shell script to check missing outputs.
   - `resubmit_job_GridFit_V6.sh`: Shell script to resubmit failed runs manually.
   - `auto_resubmit_job_GridFit_V6.sh`: Automatically finds and resubmits missing jobs.
   - `cancel_job.sh`: Cancel a range of job IDs if needed.
   - `comb_GridFit_out_files_V6.ipynb`: Jupyter Notebook to combine and organize outputs into analysis-ready format.
   - `visualize_GridFit_V6.ipynb`: Jupyter Notebook to visualize statistical results.
All code uses the Julia programming language (with Julia kernel for Jupyter notebooks).

## Output Files
Each GridFit run produces three main output files:
 - `*_GridFitFixedCoords.csv`: Coordinates from initial fixed fit.
 - `*_GridFitOptCoords.csv`: Coordinates after optimization.
 - `*_GridFitErrbyDist.csv`: Distance error for each grid point.

## How to Reproduce This Analysis

### 1. Set up environment:
 - Install Julia.
 - Ensure access to a SLURM cluster.

### 2. Download the source data from https://osf.io/k8v9u/.

### 3. Run batch jobs:
 - Use `run_testGridFit_V6.sh` to submit 300 jobs per subject/point set.

### 4. Check and resubmit failed runs (if necessary):
 - Run `find_missing_files.sh` and/or `auto_rerun_missing_jobs.sh`.

### 5. Aggregate outputs:
 - Run `comb_GridFit_out_files_V6.ipynb`.

### 6. Visualize results:
 - Run `visualize_GridFit_V6.ipynb`.

### 7. Review final figures or check ready-made plots at https://osf.io/gqnef/.

## Notes
 - Each subject is processed separately.
 - Each configuration (number of grid points) is tested independently.
 - Repetitions: 300 times per configuration.
 - Job scripts are modular to make monitoring, resubmission, and error recovery efficient.
