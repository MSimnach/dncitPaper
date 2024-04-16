#!/bin/bash -eux
#SBATCH --job-name=CITs_benchmarking
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.simnacher@guest.hpi.de
#SBATCH --partition=cpu # -p
#SBATCH --cpus-per-task=8 # -c
#SBATCH --mem=256gb
#SBATCH --gpus=0
#SBATCH --time=72:00:00 
#SBATCH --output=logs/cond_dep_embeddings_%j.log # %j is job id

set +u
eval "$(conda shell.bash hook)"
conda activate /dhc/home/marco.simnacher/conda3/envs/CIT_R

# args =['/CI/ or /No_CI/', post_non_lin, eps_sigmaX, eps_sigmaY, eps_sigmaZ, embedding_orig, embedding_obs, confounder, response, cond_dep_measure, cond_dep_measure_params]
## same embedding
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 1 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 50 0 fastsurfer fastsurfer AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 50 0 fastsurfer fastsurfer AS simulated 

# ## fastsurfer and condVAE
Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 1 0 fastsurfer condVAE AS simulated 
Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 1 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 50 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 50 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 1 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 1 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 50 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 50 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 1 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 1 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 50 0 fastsurfer condVAE AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 50 0 fastsurfer condVAE AS simulated 

# ## noisy embedding with increasing noise
# # 1 funct.
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 10 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 10 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 50 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 50 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 500 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 500 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 10 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 10 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 50 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 50 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 500 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 500 50 0 fastsurfer noisy AS simulated 
# # 4 funct.
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 10 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 10 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 50 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 50 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 500 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 500 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 10 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 10 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 50 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 50 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 500 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 500 50 0 fastsurfer noisy AS simulated 
# # 5 funct.
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 10 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 10 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 50 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 50 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 500 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 500 1 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 10 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 10 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 50 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 50 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 500 50 0 fastsurfer noisy AS simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 500 50 0 fastsurfer noisy AS simulated 

# ### All for confounder with 10 PCs
# ## same embedding
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 1 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 50 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 50 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 1 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 1 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 50 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 50 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 1 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 1 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 50 0 fastsurfer fastsurfer genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 50 0 fastsurfer fastsurfer genes10 simulated 

# ## fastsurfer and condVAE
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 1 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 1 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 0 50 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 0 50 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 1 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 1 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 0 50 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 0 50 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 1 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 1 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 0 50 0 fastsurfer condVAE genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 0 50 0 fastsurfer condVAE genes10 simulated 

# # ## noisy embedding with increasing noise
# # # 1 funct.
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 10 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 10 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 50 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 50 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 500 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 500 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 10 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 10 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 50 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 50 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 1 500 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 1 500 50 0 fastsurfer noisy genes10 simulated 
# # 4 funct.
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 10 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 10 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 50 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 50 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 500 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 500 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 10 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 10 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 50 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 50 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 4 500 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 4 500 50 0 fastsurfer noisy genes10 simulated 
# # # 5 funct.
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 10 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 10 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 50 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 50 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 500 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 500 1 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 10 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 10 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 50 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 50 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /CI/ 5 500 50 0 fastsurfer noisy genes10 simulated 
# Rscript src/cond_dep/cond_dep_embeddings.R /No_CI/ 5 500 50 0 fastsurfer noisy genes10 simulated 

set -u