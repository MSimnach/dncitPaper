# load T1 paths for diagnostic
library(data.table)
library(dplyr)
library(devtools)
load_all()
library(reticulate)
np <- import("numpy")

# load UKB image_path.csv 
readRenviron(".Renviron")
Sys.getenv(c("UKB_PATH","IXI_PATH"))
ukb_path <- Sys.getenv("UKB_PATH", unset = NA)
t1_paths_diag <- read.csv(paste0(ukb_path, "/t1_paths_diag.csv"))

### Load embeddings
# load fastsurfer, condVAE, freesurfer embeddings
fastsurfer_embedding <- data.table::fread(paste0(ukb_path, "/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
colnames(fastsurfer_embedding)[1] <- "id"
condVAE_embedding <- data.table::fread(paste0(ukb_path, "/ukb_condVAE.csv"), header=FALSE, nThread = 1)
colnames(condVAE_embedding)[1] <- "id"
freesurfer_embedding <- data.table::fread(paste0(ukb_path, "/ukb_freesurfer.csv"), header=TRUE, nThread = 1)
colnames(freesurfer_embedding)[1] <- "id"

# Extract pretrained MedicalNet embeddings
out_dir <- paste0(ukb_path, "/medicalnet_embeddings")
system2("python", c("inst/learn_embedding/extract_embeddings.py",
                    "--pretrained", "--model_name", "resnet18",
                    "--input_csv", paste0(ukb_path, "/t1_paths_diag.csv"),
                    "--output_dir", out_dir))
medicalnet_emb <- np$load(paste0(out_dir, "/embeddings.npy"))
medicalnet_idx <- fread(paste0(out_dir, "/embeddings_index.csv"))
medicalnet_embedding <- cbind(id = medicalnet_idx$subject_id, as.data.frame(medicalnet_emb))

# Extract embeddings from trained checkpoint
checkpoint_path <- paste0(ukb_path, "/t1/embeddings/550_scratch_ukb_z4_linear/1/best_model.ckpt")
out_dir_trained <- paste0(ukb_path, "/trained_model_embeddings")
system2("python", c("inst/learn_embedding/extract_embeddings.py",
                    "--checkpoint", checkpoint_path,
                    "--input_csv", paste0(ukb_path, "/t1_paths_diag.csv"),
                    "--output_dir", out_dir_trained,
                    "--config_name", "trained", 
                    "--batch_size", "16",
                    "--num_workers", "4"))
trained_emb <- np$load(paste0(out_dir_trained, "/trained/embeddings.npy"))
trained_idx <- fread(paste0(out_dir_trained, "/trained/embeddings_index.csv"))
trained_embedding <- cbind(id = trained_idx$subject_id, as.data.frame(trained_emb))

Y <- read.csv(paste0(ukb_path, "/ukb_z1_age.csv"))

# Find common IDs across all datasets
# Assuming Y is stored in a dataframe/data.table with an 'id' column
common_ids <- Reduce(intersect, list(
  t1_paths_diag$id,
  fastsurfer_embedding$id,  # or $id depending on column name
  #condVAE_embedding$id,
  freesurfer_embedding$id,
  medicalnet_embedding$id,
  trained_embedding$id,
  Y$id  # adjust column name as needed
))

cat("Total common IDs across all datasets:", length(common_ids), "\n")

# Filter and sort all datasets by common IDs
# Sort common_ids first for consistent ordering
common_ids <- sort(common_ids)

# Filter and reorder each dataset
t1_paths_diag <- t1_paths_diag[match(common_ids, t1_paths_diag$id), ]
fastsurfer_emb <- fastsurfer_embedding[match(common_ids, fastsurfer_embedding$eid), ]  # adjust column name
condVAE_emb <- condVAE_embedding[match(common_ids, condVAE_embedding$id), ]
freesurfer_emb <- freesurfer_embedding[match(common_ids, freesurfer_embedding$id), ]
medicalnet_emb <- medicalnet_embedding[match(common_ids, medicalnet_embedding$id), ]
trained_emb <- trained_embedding[match(common_ids, trained_embedding$id), ]
Y_aligned <- Y[match(common_ids, Y$id), ]

# Verify alignment (all should return TRUE)
cat("Checking alignment:\n")
cat("  t1_paths_diag IDs match:", all(t1_paths_diag$id == common_ids), "\n")
cat("  fastsurfer IDs match:", all(fastsurfer_emb$eid == common_ids), "\n")
#cat("  condVAE IDs match:", all(condVAE_emb$id == common_ids), "\n")
cat("  freesurfer IDs match:", all(freesurfer_emb$id == common_ids), "\n")
cat("  medicalnet IDs match:", all(medicalnet_emb$id == common_ids), "\n")
cat("  trained IDs match:", all(trained_emb$id == common_ids), "\n")
cat("  Y IDs match:", all(Y_aligned$id == common_ids), "\n")

