### Parallelisation of the simulation study for the UK Biobank brain MRI data
library(doParallel)
library(foreach)
library(DNCIT)
library('dncitPaper')
data.table::setDTthreads(1)
args = commandArgs(trailingOnly=TRUE)
cat(args[10])

n_cits <- 1
cit <- c(args[10])
print(args)

####### In parallel #######
if(cit == 'KCIT' || tail(args,1)=='20' || cit=='CMIknn'){
  n_sample = list(145, 256, 350, 460, 825, 1100)
}else if(cit=='WALD'){
  n_sample = list(350, 460, 825, 1100, 1475, 1964, 5000, 10000)
}else if(cit=='pred_cit'){
  n_sample = list(460, 825, 1100)
}else{
  n_sample = list(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)
}
n_seeds = 1:200

#vary also effect sizes -> change also how to save results in for loops
#beta2s_all <- list()
#if(args[7]=="condVAE"){
#  for (k in -2:2){
#    beta2s_all <- append(beta2s_all, c(10,7.5,5, 2.5)*10^(-k))
#  }
#}else{
#  for (k in -1:3){
#    beta2s_all <- append(beta2s_all, c(10,7.5,5, 2.5)*10^(-k))
#  }
#}
#beta2s <- beta2s_all

beta2s <- list(1)

if (cit %in% c('cpi', 'pred_cit')) {
   pkgs_for_each <- c('DNCIT', 'dncitPaper', 'mlr3', 'mlr3learners')
} else {
   pkgs_for_each <- c('DNCIT', 'dncitPaper')
}

#if(cit %in% c("FCIT")){
#  Sys.setenv(OMP_NUM_THREADS = "50")
#  cl <- parallel::makeCluster(50, outfile="")
#  doParallel::registerDoParallel(cl)
#}else if(cit %in% c("CMIknn")){
#  Sys.setenv(OMP_NUM_THREADS = "50")
#  cl <- parallel::makeCluster(50, outfile="")
#  doParallel::registerDoParallel(cl)
#}else{
Sys.setenv(OMP_NUM_THREADS = "50")
cl <- parallel::makeCluster(50, outfile="")
doParallel::registerDoParallel(cl)
#}

#cl <- parallel::makeCluster(10, outfile="")
#doParallel::registerDoParallel(cl)

res_time <- foreach::foreach (i= n_seeds, .packages = pkgs_for_each) %dopar% {
  if(cit %in% c('pred_cit')){
    # Set the logging level for mlr3
   lgr::get_logger("mlr3")$set_threshold("fatal")
   # Set the logging level for bbotk
   lgr::get_logger("bbotk")$set_threshold("fatal")
  }else if (cit %in% c('cpi')){
    lgr::get_logger("cpi")$set_threshold("fatal")
  }
                                                 if (grepl('/CI',args[1],fixed=TRUE)){
                                                   #cat('CIT:', cit)
                                                   res <- rep(0,length(n_sample))
                                                   runtime <- rep(0,length(n_sample))
                                                   for (idx_sample in seq_along(n_sample)){
                                                     if(i %% 50 == 1) {
                                                      cat(paste("Iteration",i,"for sample size", n_sample[[idx_sample]], "\n"))
                                                     }
                                                     XYZ_list <- dncitPaper::data_gen(seed=i, idx_sample=idx_sample, n_sample=n_sample, idx_beta2=NULL, beta2s=NULL,
                                                                          post_non_lin=as.numeric(args[2]), eps_sigmaX=as.numeric(args[3]), eps_sigmaY=as.numeric(args[4]),
                                                                          eps_sigmaZ=as.numeric(args[5]), embedding_orig=args[6], embedding_obs=args[7],
                                                                          confounder=args[8], g_z=args[9])
                                                     X <- as.matrix(XYZ_list[[1]])
                                                     Y <- as.matrix(XYZ_list[[2]])
                                                     Z <- as.matrix(XYZ_list[[3]])

                                                     if (args[10] == 'RCOT'){
                                                       cit_params <- list(cit='RCOT', params_cit=list(seed=as.numeric(args[11])))
                                                     }else if(args[10] == 'CMIknn'){
                                                       cit_params <- list(cit='cmiknn', params_cit=list())
                                                     }else if(args[10] == 'kpc_graph'){
                                                      if (args[11]=='1'){
                                                          k = kernlab::vanilladot()
                                                      }else if (args[11]=='2') {
                                                          k = kernlab::rbfdot(1/(2*stats::median(stats::dist(X))^2))
                                                      }else if (args[11]=='3') {
                                                          k = kernlab::laplacedot(1/(2*stats::median(stats::dist(X))^2))
                                                      }else if (args[11]=='4') {
                                                          k = kernlab::tanhdot()
                                                      }
                                                      if(args[8]=='ukb_z1'){
                                                        model_formula_YZ <- "V1~1+s(V2)"
                                                      }else if(args[8]=='ukb_z2'){
                                                        model_formula_YZ <- 'V1~1+s(V2)+s(V3)'
                                                      }else if(args[8]=='ukb_z4'){
                                                        n_covs <- ncol(Z)
                                                        lin_covs <- paste0("V", seq(4, n_covs+1))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste0('V1~1+s(V2)+s(V3)+', lin_covs_string)
                                                      }else if(args[8]=='ukb_z6'){
                                                        n_covs <- ncol(Z)
                                                        date_diff <- as.character(n_covs)
                                                        qc <- as.character(n_covs+1)
                                                        lin_covs <- paste0("V", seq(4, n_covs-1))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+', lin_covs_string, sep="")
                                                      }else if(args[8]=='ukb_z10'){
                                                        n_covs <- ncol(Z)
                                                        date_diff <- as.character(n_covs-4)
                                                        qc <- as.character(n_covs-3)
                                                        head_loc_1 <- as.character(n_covs-2)
                                                        head_loc_2 <- as.character(n_covs-1)
                                                        head_loc_3 <- as.character(n_covs)
                                                        head_loc_4 <- as.character(n_covs+1)
                                                        lin_covs <- paste0("V", seq(4, n_covs-5))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+s(V', head_loc_1, ', k=3)+s(V', head_loc_2, ', k=3)+s(V', head_loc_3, ', k=3)+s(V', head_loc_4, ', k=3)+', lin_covs_string, sep="")
                                                      }else if(args[8]=='ukb_z15'){
                                                        n_covs <- ncol(Z)
                                                        date_diff <- as.character(n_covs-9)
                                                        qc <- as.character(n_covs-8)
                                                        head_loc_1 <- as.character(n_covs-7)
                                                        head_loc_2 <- as.character(n_covs-6)
                                                        head_loc_3 <- as.character(n_covs-5)
                                                        head_loc_4 <- as.character(n_covs-4)
                                                        gene_1 <- as.character(n_covs-3)
                                                        gene_2 <- as.character(n_covs-2)
                                                        gene_3 <- as.character(n_covs-1)
                                                        gene_4 <- as.character(n_covs)
                                                        gene_5 <- as.character(n_covs+1)
                                                        lin_covs <- paste0("V", seq(4, n_covs-10))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+s(V', head_loc_1, ', k=3)+s(V', head_loc_2, ', k=3)+s(V', head_loc_3, ', k=3)+s(V', head_loc_4, ', k=3)+s(V', gene_1, ', k=3)+s(V', gene_2, ', k=3)+s(V', gene_3, ', k=3)+s(V', gene_4, ', k=3)+s(V', gene_5, ', k=3)+', lin_covs_string, sep="")
                                                      }
                                                      cit_params <- list(cit='cpt_kpc', params_cit=list(k=k, Knn = as.numeric(args[12]), model.formula.YZ=model_formula_YZ))
                                                     }else if(args[10]=='FCIT'){
                                                        cit_params <- list(cit='fcit')
                                                     }else if(args[10]=='cpi'){
                                                        cit_params <- list(cit='cpi')
                                                     }else if(args[10]=='comets_pcm'){
                                                       cit_params <- list(cit='comets')
                                                     }else if(args[10]=='comets_gcm'){
                                                       cit_params <- list(cit='comets', params_cit=list(method='gcm', alternative='less'))
                                                     }else if(args[10]=='pred_cit'){
                                                        min_samples <- min(unlist(n_sample))
                                                        max_samples <- max(unlist(n_sample))
                                                        current_sample <- n_sample[[idx_sample]]
                                                        term_time <- round(exp(1.5)+(current_sample-min_samples)/(max_samples-min_samples)*(exp(2.25)-exp(1.5))/3)
                                                        cit_params <- list(cit='pred_cit', params_cit=list(term_time = term_time))
                                                     }else if(args[10]=='WALD'){
                                                        cit_params <- list(cit='wald')
                                                     }

                                                     tmp <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                                           cit_with_parameters = cit_params)
                                                     res[idx_sample] <- tmp$p
                                                     runtime[idx_sample] <- tmp$runtime
                                                   }
                                                   p_time <- cbind(res, runtime)
                                                 }else if(grepl('/No_CI',args[1],fixed=TRUE)){
                                                   #cat('CIT:', cit)
                                                   res <- rep(0,length(n_sample))
                                                   runtime <- rep(0,length(n_sample))
                                                   for (idx_sample in seq_along(n_sample)){
                                                     if(i %% 50 == 1) {
                                                      cat(paste("Iteration",i,"for sample size", n_sample[[idx_sample]], "\n"))
                                                     }
                                                     for (idx_beta2 in seq_along(beta2s)){
                                                       #cat(paste("Iteration",i, "for beta2", beta2s[[idx_beta2]], "\n"))
                                                       XYZ_list <- dncitPaper::data_gen(seed=i, idx_sample=idx_sample, n_sample=n_sample,idx_beta2=idx_beta2, beta2s=beta2s,
                                                                            post_non_lin=as.numeric(args[2]), eps_sigmaX=as.numeric(args[3]), eps_sigmaY=as.numeric(args[4]),
                                                                            eps_sigmaZ=as.numeric(args[5]), embedding_orig=args[6], embedding_obs=args[7], confounder=args[8],
                                                                            g_z=args[9])
                                                       X <- as.matrix(XYZ_list[[1]])
                                                       Y <- as.matrix(XYZ_list[[2]])
                                                       Z <- as.matrix(XYZ_list[[3]])
                                                       #XYZ <- data.frame(XYZ_list)

                                                       if (args[10] == 'RCOT'){
                                                         cit_params <- list(cit='RCOT', params_cit=list(seed=as.numeric(args[11])))
                                                       }else if(args[10] == 'CMIknn'){
                                                         cit_params <- list(cit='cmiknn', params_cit=list())
                                                       }else if(args[10] == 'kpc_graph'){
                                                          if (args[11]=='1'){
                                                              k = kernlab::vanilladot()
                                                          }else if (args[11]=='2') {
                                                              k = kernlab::rbfdot(1/(2*stats::median(stats::dist(X))^2))
                                                          }else if (args[11]=='3') {
                                                              k = kernlab::laplacedot(1/(2*stats::median(stats::dist(X))^2))
                                                          }else if (args[11]=='4') {
                                                              k = kernlab::tanhdot()
                                                          }
                                                          if(args[8]=='ukb_z1'){
                                                            model_formula_YZ <- "V1~1+s(V2)"
                                                          }else if(args[8]=='ukb_z2'){
                                                            model_formula_YZ <- 'V1~1+s(V2)+s(V3)'
                                                          }else if(args[8]=='ukb_z4'){
                                                            n_covs <- ncol(Z)
                                                            lin_covs <- paste0("V", seq(4, n_covs+1))
                                                            lin_covs_string <- paste(lin_covs, collapse = "+")
                                                            model_formula_YZ <- paste0('V1~1+s(V2)+s(V3)+', lin_covs_string)
                                                          }else if(args[8]=='ukb_z6'){
                                                            n_covs <- ncol(Z)
                                                            date_diff <- as.character(n_covs)
                                                            qc <- as.character(n_covs+1)
                                                            lin_covs <- paste0("V", seq(4, n_covs-1))
                                                            lin_covs_string <- paste(lin_covs, collapse = "+")
                                                            model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+', lin_covs_string, sep="")
                                                          }else if(args[8]=='ukb_z10'){
                                                            n_covs <- ncol(Z)
                                                            date_diff <- as.character(n_covs-4)
                                                            qc <- as.character(n_covs-3)
                                                            head_loc_1 <- as.character(n_covs-2)
                                                            head_loc_2 <- as.character(n_covs-1)
                                                            head_loc_3 <- as.character(n_covs)
                                                            head_loc_4 <- as.character(n_covs+1)
                                                            lin_covs <- paste0("V", seq(4, n_covs-5))
                                                            lin_covs_string <- paste(lin_covs, collapse = "+")
                                                            model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+s(V', head_loc_1, ', k=3)+s(V', head_loc_2, ', k=3)+s(V', head_loc_3, ', k=3)+s(V', head_loc_4, ', k=3)+', lin_covs_string, sep="")
                                                          }else if(args[8]=='ukb_z15'){
                                                            n_covs <- ncol(Z)
                                                            date_diff <- as.character(n_covs-9)
                                                            qc <- as.character(n_covs-8)
                                                            head_loc_1 <- as.character(n_covs-7)
                                                            head_loc_2 <- as.character(n_covs-6)
                                                            head_loc_3 <- as.character(n_covs-5)
                                                            head_loc_4 <- as.character(n_covs-4)
                                                            gene_1 <- as.character(n_covs-3)
                                                            gene_2 <- as.character(n_covs-2)
                                                            gene_3 <- as.character(n_covs-1)
                                                            gene_4 <- as.character(n_covs)
                                                            gene_5 <- as.character(n_covs+1)
                                                            lin_covs <- paste0("V", seq(4, n_covs-10))
                                                            lin_covs_string <- paste(lin_covs, collapse = "+")
                                                            model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+s(V', head_loc_1, ', k=3)+s(V', head_loc_2, ', k=3)+s(V', head_loc_3, ', k=3)+s(V', head_loc_4, ', k=3)+s(V', gene_1, ', k=3)+s(V', gene_2, ', k=3)+s(V', gene_3, ', k=3)+s(V', gene_4, ', k=3)+s(V', gene_5, ', k=3)+', lin_covs_string, sep="")
                                                          }
                                                          cit_params <- list(cit='cpt_kpc', params_cit=list(k=k, Knn = as.numeric(args[12]), model.formula.YZ=model_formula_YZ))
                                                       }else if(args[10]=='FCIT'){
                                                          cit_params <- list(cit='fcit')
                                                       }else if(args[10]=='cpi'){
                                                         cit_params <- list(cit='cpi')
                                                       }else if(args[10]=='comets_pcm'){
                                                         cit_params <- list(cit='comets')
                                                       }else if(args[10]=='comets_gcm'){
                                                         cit_params <- list(cit='comets', params_cit=list(method='gcm', alternative='less'))
                                                       }else if(args[10]=='pred_cit'){
                                                        min_samples <- min(unlist(n_sample))
                                                        max_samples <- max(unlist(n_sample))
                                                        current_sample <- n_sample[[idx_sample]]
                                                        term_time <- round(exp(1.5)+(current_sample-min_samples)/(max_samples-min_samples)*(exp(2.25)-exp(1.5))/3)
                                                        cit_params <- list(cit='pred_cit', params_cit=list(term_time = term_time))
                                                       }else if(args[10]=='WALD'){
                                                          cit_params <- list(cit='wald', params_cit=NULL)
                                                       }

                                                       tmp <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                                                    cit_with_parameters = cit_params)
                                                       res[idx_sample] <- tmp$p
                                                       runtime[idx_sample] <- tmp$runtime
                                                     }
                                                   }
                                                   p_time <- cbind(res, runtime)
                                                 }
                                                 p_time
                                               }
parallel::stopCluster(cl)


#save results
p_res <- matrix(nrow=length(n_seeds), ncol=length(n_sample))
runtime <- matrix(nrow=length(n_seeds), ncol=length(n_sample))
for (i in n_seeds){
  for (idx_sample in seq_along(n_sample)){
    p_res[i, idx_sample] <- res_time[[i]][idx_sample,1]
    runtime[i, idx_sample] <- res_time[[i]][idx_sample,2]
  }
}
write.csv(p_res,file=paste0("Results", args[1], "p-values/", paste(args[-1], collapse="_"), ".csv", collapse=''))
write.csv(runtime,file=paste0("Results", args[1], "Runtime/", paste(args[-1], collapse="_"), ".csv", collapse=''))

#rejection rates
n_seeds_ <- length(n_seeds)
for (idx_sample in seq_along(n_sample)){
  rejected <- colSums(p_res < 0.05, na.rm=TRUE) / n_seeds_
}
rejected <- data.frame(rejected, row.names = n_sample)
colnames(rejected) <- cit
#print(xtable(rejected))
write.csv(rejected,file=paste0("Results", args[1], "rejection_rates/", paste(args[-1], collapse="_"), ".csv", collapse=''))


if (FALSE){
  if(grepl("/CI",args[1],fixed=TRUE)){
    p_res <- matrix(nrow=length(n_seeds), ncol=length(n_sample))
    runtime <- matrix(nrow=length(n_seeds), ncol=length(n_sample))
    for (i in n_seeds){
      for (idx_sample in seq_along(n_sample)){
        p_res[i, idx_sample] <- res_time[[i]][idx_sample,1]
        runtime[i, idx_sample] <- res_time[[i]][idx_sample,2]
      }
    }
    write.csv(p_res,file=paste0("Results", args[1], "p-values/", paste(args[-1], collapse="_"), ".csv", collapse=''))
    write.csv(runtime,file=paste0("Results", args[1], "Runtime/", paste(args[-1], collapse="_"), ".csv", collapse=''))

    #rejection rates
    n_seeds_ <- length(n_seeds)
    for (idx_sample in seq_along(n_sample)){
      rejected <- colSums(p_res < 0.05, na.rm=TRUE) / n_seeds_
    }
    rejected <- data.frame(rejected, row.names = n_sample)
    colnames(rejected) <- cit
    #print(xtable(rejected))
    write.csv(rejected,file=paste0("Results", args[1], "rejection_rates/", paste(args[-1], collapse="_"), ".csv", collapse=''))

    # #ks statistics
    # ks <- matrix(nrow = length(n_sample), ncol=n_cits)
    # for (idx_sample in seq_along(n_sample)){
    #   ks[idx_sample,1] <- ks.test(p_res[,idx_sample], punif)$statistic
    # }
    # ks_df <- data.frame(ks)
    # colnames(ks_df)<- cit
    # row.names(ks_df) <- n_sample
    # write.csv(ks_df,file=paste0("Results", args[1], "Kolmogorov-Smirnov-Statistics/", paste(args[-1], collapse="_"), ".csv", collapse=''))
  }else if(grepl("/No_CI",args[1],fixed=TRUE)){
    #save p-values and runtimes
    p_res <- matrix(nrow=length(n_seeds), ncol=length(beta2s))
    runtime <- matrix(nrow=length(n_seeds), ncol=length(beta2s))
    for (i in n_seeds){
      for (idx_beta2 in seq_along(beta2s)){
        p_res[i, idx_beta2] <- res_time[[i]][idx_beta2,1]
        runtime[i, idx_beta2] <- res_time[[i]][idx_beta2,2]
      }
    }
    write.csv(p_res,file=paste0("Results", args[1], "p-values/", paste(args[-1], collapse="_"), ".csv", collapse=''))
    write.csv(runtime,file=paste0("Results", args[1], "Runtime/", paste(args[-1], collapse="_"), ".csv", collapse=''))

    #rejection rates
    n_seeds_ <- length(n_seeds)
    rejected_ <- colSums(p_res < 0.05, na.rm=TRUE) / n_seeds_
    rejected <- data.frame(rejected_, row.names = beta2s)
    colnames(rejected) <- cit
    #print(xtable(rejected))
    write.csv(rejected,file=paste0("Results", args[1], "rejection_rates/", paste(args[-1], collapse="_"), ".csv", collapse=''))

  }
}
