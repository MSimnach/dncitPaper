#!/bin/bash -eux
#SBATCH --job-name=CITs_benchmarking
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.simnacher@guest.hpi.de
#SBATCH --partition=cpu # -p
#SBATCH --cpus-per-task=1 # -c
#SBATCH --mem=128gb
#SBATCH --gpus=0
#SBATCH --time=02:00:00 
#SBATCH --output=logs/Y_gen_%j.log # %j is job id

set +u
eval "$(conda shell.bash hook)"
conda activate /dhc/home/marco.simnacher/conda3/envs/CIT_R

#### args =['/CI/ or /No_CI/', post_non_lin, eps_sigmaX, eps_sigmaY, eps_sigmaZ, embedding_orig, embedding_obs, confounder, response, CIT, CIT_params]
## same embedding
Rscript src/Y_generation.R /CI/ 1 0 1 0 fastsurfer fastsurfer AS simulated 
Rscript src/Y_generation.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /CI/ 1 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /No_CI/ 1 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /CI/ 4 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /No_CI/ 4 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /CI/ 4 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /No_CI/ 4 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /CI/ 5 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /No_CI/ 5 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /CI/ 5 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/Y_generation.R /No_CI/ 5 0 50 0 fastsurfer fastsurfer AS simulated 
set -u