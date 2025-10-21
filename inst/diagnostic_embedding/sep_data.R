#### prepare separate data for diagnostic of embeddings
library(data.table)
library(dplyr)
library(devtools)
load_all()

# load UKB image_path.csv 
readRenviron(".Renviron")
Sys.getenv(c("UKB_PATH","IXI_PATH"))
ukb_path <- Sys.getenv("UKB_PATH", unset = NA)
ukb_data <- read.csv(paste0(ukb_path, "/t1_paths.csv"))

n_diag <- 1000
idx <- sample(1:nrow(ukb_data), n_diag)
ukb_data_diag <- ukb_data[idx,]

# save ukb_data_diag
write.csv(ukb_data_diag, file = paste0(ukb_path, "/t1_paths_diag.csv"), row.names = FALSE)

idx_cit <- !(1:nrow(ukb_data) %in% idx)
ukb_data_for_cit <- ukb_data[idx_cit,]
write.csv(ukb_data_for_cit, file = paste0(ukb_path, "/t1_paths_cit.csv"), row.names = FALSE)
