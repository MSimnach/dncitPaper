#!/bin/bash -eux
#SBATCH --job-name=CITs_benchmarking
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.simnacher@guest.hpi.de
#SBATCH --partition=gpua100 # -p
#SBATCH --cpus-per-task=10 # -c
#SBATCH --mem=128gb
#SBATCH --gpus=0
#SBATCH --time=72:00:00 
#SBATCH --output=logs/kpc_graph_3_1_%j.log # %j is job id

set +u
eval "$(conda shell.bash hook)"
conda activate /dhc/home/marco.simnacher/conda3/envs/CIT_R

# args =['/CI/ or /No_CI/', post_non_lin, eps_sigmaX, eps_sigmaY, eps_sigmaZ, embedding_orig, embedding_obs, confounder, response, CIT, CIT_params]
## same embedding
Rscript src/R_CITs.R /CI/ 1 0 1 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 0 50 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 50 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 1 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 1 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 50 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 50 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 1 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 1 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 50 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 50 0 fastsurfer fastsurfer AS simulated kpc_graph 3 1

## fastsurfer and condVAE
Rscript src/R_CITs.R /CI/ 1 0 1 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 1 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 0 50 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 50 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 1 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 1 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 50 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 50 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 1 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 1 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 50 0 fastsurfer condVAE AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 50 0 fastsurfer condVAE AS simulated kpc_graph 3 1

## noisy embedding with increasing noise
# 1 funct.
Rscript src/R_CITs.R /CI/ 1 10 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 10 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 50 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 50 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 500 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 500 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 10 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 10 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 50 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 50 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 500 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 500 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
# 4 funct.
Rscript src/R_CITs.R /CI/ 4 10 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 10 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 50 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 50 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 500 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 500 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 10 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 10 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 50 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 50 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 500 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 500 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
# 5 funct.
Rscript src/R_CITs.R /CI/ 5 10 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 10 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 50 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 50 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 500 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 500 1 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 10 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 10 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 50 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 50 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 500 50 0 fastsurfer noisy AS simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 500 50 0 fastsurfer noisy AS simulated kpc_graph 3 1

### All for confounder with 10 PCs
## same embedding
Rscript src/R_CITs.R /CI/ 1 0 1 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 0 50 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 50 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 1 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 1 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 50 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 50 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 1 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 1 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 50 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 50 0 fastsurfer fastsurfer genes10 simulated kpc_graph 3 1

## fastsurfer and condVAE
Rscript src/R_CITs.R /CI/ 1 0 1 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 1 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 0 50 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 0 50 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 1 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 1 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 0 50 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 0 50 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 1 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 1 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 0 50 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 0 50 0 fastsurfer condVAE genes10 simulated kpc_graph 3 1

# ## noisy embedding with increasing noise
# # 1 funct.
Rscript src/R_CITs.R /CI/ 1 10 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 10 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 50 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 50 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 500 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 500 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 10 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 10 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 50 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 50 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 1 500 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 1 500 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
# 4 funct.
Rscript src/R_CITs.R /CI/ 4 10 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 10 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 50 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 50 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 500 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 500 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 10 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 10 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 50 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 50 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 4 500 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 4 500 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
# # 5 funct.
Rscript src/R_CITs.R /CI/ 5 10 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 10 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 50 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 50 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 500 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 500 1 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 10 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 10 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 50 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 50 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /CI/ 5 500 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
Rscript src/R_CITs.R /No_CI/ 5 500 50 0 fastsurfer noisy genes10 simulated kpc_graph 3 1
set -u
