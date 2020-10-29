#!/bin/sh
#SBATCH --job-name=flash_array   # Job name
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=2gb           # Memory per processor
#SBATCH --time=10-00:00:00          # Time limit days-hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=6-10                 # Array range
#SBATCH --account=lsa1
#SBATCH --partition=standard


# Run the loop of runs for this task.
python3.6 ./notebooks/flash_scan.py 10 $SLURM_ARRAY_TASK_ID
