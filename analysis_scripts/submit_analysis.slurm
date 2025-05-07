#!/usr/bin/env bash
#
# SLURM submission script for processing VR data workflows
#
# Usage: sbatch this_script.sh
#
# This script runs a series of Python analysis steps under SLURM,
# auto-detects the “adiabatic” mode based on its directory name,
# and logs each step with timestamps and error trapping.

#SBATCH --job-name=processing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=cosma8
#SBATCH --account=username
#SBATCH --exclusive
#SBATCH --time=20:00:00
#SBATCH --output=../%j.analysis.out
#SBATCH --error=../%j.analysis.err

#### Strict error handling & tracing ####
set -euo pipefail
# Uncomment for debugging:
# set -x

#### Environment setup ####
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#### Helper functions ####
log() {
    local ts; ts=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[$ts] $*"
}

on_error() {
    local exit_code=$?
    log "ERROR: Script exited with code ${exit_code}."
    sacct -j "$SLURM_JOBID" --format=JobID,JobName,Partition,AveRSS,MaxRSS,AveVMSize,MaxVMSize,Elapsed,ExitCode
    exit $exit_code
}
trap on_error ERR

#### Detect adiabatic mode ####
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $script_dir == *adiabatic* ]]; then
    adiabatic_flag="--adiabatic"
    log "Adiabatic mode detected."
else
    adiabatic_flag=""
    log "Standard mode."
fi

#### List of analysis steps ####
python_steps=(
    "rename_vr_mpi.py"
    "clean_vr_submit.py"
    "final_halo_centre.py"
    "final_gas_ids.py"
    "plot_runtime.py"
    "phase_space.py"
    "basic_properties.py ${adiabatic_flag}"
    "particle_history.py"
    "vis_temperatures.py"
)

#### Main execution loop ####
log "Starting analysis pipeline."
for step in "${python_steps[@]}"; do
    log "Running: python3 ${step}"
    python3 ${step}
    log "Completed: ${step}"
done

#### Final reporting ####
log "Job done. Gathering SLURM accounting information:"
sacct -j "$SLURM_JOBID" --format=JobID,JobName,Partition,AveRSS,MaxRSS,AveVMSize,MaxVMSize,Elapsed,ExitCode
