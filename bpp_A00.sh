#!/bin/bash
#
#SBATCH --job-name=bpp_A00
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --gres=localtmp:100G
#SBATCH --mem=96G
#SBATCH -t 5-00:00:00
#SBATCH -o /home/dmendez/outdir/%j_%x
#SBATCH -e /home/dmendez/errdir/%j_%x

set -eo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PROC_BIND=false
export OMP_PLACES=threads

cd /fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/OrthoFinder/Results_Jan09/bpp_from_SCO/A00

/fast/AG_Lewin/dmendez/tools/bpp-4.8.7-linux-x86_64/bin/bpp --no-pin --cfile /fast/AG_Lewin/dmendez/orthofinder/20251219/primary/input_clean/OrthoFinder/Results_Jan09/bpp_from_SCO/A00/A00.ctl