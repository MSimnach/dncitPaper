#!/bin/bash -eux
#SBATCH --job-name=CITs_benchmarking
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.simnacher@guest.hpi.de
#SBATCH --partition=cpu,gpu,gpupro
 # -p
#SBATCH --cpus-per-task=50 # -c
#SBATCH --mem=128gb
#SBATCH --gpus=0
#SBATCH --time=72:00:00 
#SBATCH --output=logs/comets_pcm_%j.log # %j is job id

set +u
eval "$(conda shell.bash hook)"
conda activate /dhc/home/marco.simnacher/conda3/envs/install-dncit

# args = [dependence fct_relation eps_sigmaX eps_sigmaY eps_sigmaZ embedding_orig embedding_obs confounder response cit cit params]
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 linear comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 realistic comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z1 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z2 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z4 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z4 squared comets_pcm
#Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z4 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z4 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z10 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z10 squared comets_pcm
#Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z10 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z10 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z15 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 3 1 0 fastsurfer noisy ukb_z15 squared comets_pcm
#Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer freesurfer ukb_z15 squared comets_pcm
Rscript raw_code/sim_ukb_brainmri.R /CI/ 1 0 1 0 fastsurfer condVAE ukb_z15 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z1 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z2 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z2 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z2 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z2 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z4 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z4 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z4 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z4 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z6 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z10 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z10 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z10 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z10 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer fastsurfer ukb_z15 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 3 1 0 fastsurfer noisy ukb_z15 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer freesurfer ukb_z15 squared comets_pcm
# Rscript raw_code/sim_ukb_brainmri.R /No_CI/ 1 0 1 0 fastsurfer condVAE ukb_z15 squared comets_pcm
set -u
