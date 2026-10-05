#!/bin/bash
#SBATCH --job-name=atlas_struct_src
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --array=0-23
#SBATCH --output=<scratch>/atlas_structure_ca/logs/%x_%A_%a.out
#SBATCH --error=<scratch>/atlas_structure_ca/logs/%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>

# TODO 34a. 24 shards x ~1,700 models; the pydssp job this mirrors took ~2 min per shard.
# Output lives OUTSIDE the repo: <scratch>/atlas_structure_ca/ca_NN.jsonl.gz
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
module load Miniforge3/24.11.3-0
source activate <home>/envs/esmfold_cor

python <scratch>/poplar-717-atlas/scripts/make_structure_source.py \
  <scratch>/atlas_structure_ca "${SLURM_ARRAY_TASK_ID}" 24
