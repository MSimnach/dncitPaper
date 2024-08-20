#!/bin/bash -eux
#SBATCH --job-name=CITs_benchmarking
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.simnacher@guest.hpi.de
#SBATCH --partition=gpu,gpupro,cpu,gpua100 # -p
#SBATCH --cpus-per-task=50 # -c
#SBATCH --mem=128gb
#SBATCH --gpus=0
#SBATCH --time=72:00:00 
#SBATCH --output=logs/kpc_graph_2_10_%j.log # %j is job id

set +u
eval "$(conda shell.bash hook)"
conda activate /dhc/home/marco.simnacher/conda3/envs/install-dncit

# args = cit script dependence fct_relation eps_sigmaX eps_sigmaY eps_sigmaZ embedding_orig embedding_obs confounder response CIT
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 linear kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 realistic kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 breakpoint3 kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z4 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z4 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z4 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z4 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z6 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z10 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z10 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z10 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z10 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z15 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z15 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z15 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z15 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z1 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 squared kpc_graph 2 10
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z4 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z4 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z4 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z4 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z6 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z10 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z10 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z10 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z10 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z15 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z15 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z15 squared kpc_graph 2 10
Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z15 squared kpc_graph 2 10
set -u
