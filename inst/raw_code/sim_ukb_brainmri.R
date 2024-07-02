library(doParallel)
library(foreach)
#devtools::load_all('/dhc/home/marco.simnacher/DNCIT')
#devtools::install('/home/RDC/simnacma/Coding/DNCIT')
library(DNCIT)
#devtools::load_all('/dhc/home/marco.simnacher/dncitPaper')
#devtools::install('/home/RDC/simnacma/Coding/dncitPaper')
library('dncitPaper')
data.table::setDTthreads(1)
args = commandArgs(trailingOnly=TRUE)
cat(args[10])

n_cits <- 1
cit <- c(args[10])
print(args)

####### In parallel #######
if(cit == 'KCIT' || tail(args,1)=='10' || tail(args,1)=='20' || cit=='CMIknn'){
  n_sample = list(50, 100, 500, 1000)
}else if(cit=='WALD' || cit=='RCOT'){
  n_sample = list(500, 1000, 5000, 10000)
}else{
  n_sample = list(50, 100, 500, 1000, 5000, 10000)
}
n_seeds = 1:200
n <- 1000
beta2s_all <- list()
if(args[7]=='condVAE'){
  for (k in -2:2){
    beta2s_all <- append(beta2s_all, c(10,7.5,5, 2.5)*10^(-k))
  }
}else{
  for (k in -1:3){
    beta2s_all <- append(beta2s_all, c(10,7.5,5, 2.5)*10^(-k))
  }
}
beta2s <- beta2s_all

if(cit %in% c('CMIknn', 'FCIT')){
  Sys.setenv(OMP_NUM_THREADS = "50")
  cl <- parallel::makeCluster(2, outfile="")
  doParallel::registerDoParallel(cl)
}else{
  Sys.setenv(OMP_NUM_THREADS = "50")
  cl <- parallel::makeCluster(50, outfile="")
  doParallel::registerDoParallel(cl)
}

res_time <- foreach::foreach (i= n_seeds, .packages = c('DNCIT', 'dncitPaper')) %dopar% {
                                                 if (grepl('/CI',args[1],fixed=TRUE)){
                                                   cat('CIT:', cit)
                                                   res <- rep(0,length(n_sample))
                                                   runtime <- rep(0,length(n_sample))
                                                   for (idx_sample in seq_along(n_sample)){
                                                     cat(paste("Iteration",i,"for sample size", n_sample[[idx_sample]], "\n"))
                                                     XYZ_list <- dncitPaper::data_gen(seed=i, idx_sample=idx_sample, n_sample=n_sample, idx_beta2=NULL, beta2s=NULL, n=NULL,
                                                                          post_non_lin=as.numeric(args[2]), eps_sigmaX=as.numeric(args[3]), eps_sigmaY=as.numeric(args[4]),
                                                                          eps_sigmaZ=as.numeric(args[5]), embedding_orig=args[6], embedding_obs=args[7],
                                                                          confounder=args[8], response=args[9])
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
                                                      if(args[8]=='AS'){
                                                        model_formula_YZ <- "V1~1+s(V3, by=V2)"
                                                      }else if(args[8]=='genes10'){
                                                        add_confounders <- paste('+',paste('s(V', 4:(ncol(Y)+ncol(Z)), ')', collapse='+', sep=""), sep="")
                                                        model_formula_YZ <- paste('V1~1+s(V3, by=V2)', add_confounders, sep="")
                                                      }
                                                      cit_params <- list(cit='cpt_kpc', params_cit=list(k=k, Knn = as.numeric(args[12]), model.formula.YZ=model_formula_YZ))
                                                     }else if(args[10]=='FCIT'){
                                                        cit_params <- list(cit='fcit')
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
                                                   cat('CIT:', cit)
                                                   res <- rep(0,length(beta2s))
                                                   runtime <- rep(0,length(beta2s))
                                                   for (idx_beta2 in seq_along(beta2s)){
                                                     cat(paste("Iteration",i, "for beta2", beta2s[[idx_beta2]], "\n"))
                                                     XYZ_list <- dncitPaper::data_gen(seed=i, idx_sample=NULL, n_sample=NULL,idx_beta2=idx_beta2, beta2s=beta2s, n=n,
                                                                          post_non_lin=as.numeric(args[2]), eps_sigmaX=as.numeric(args[3]), eps_sigmaY=as.numeric(args[4]),
                                                                          eps_sigmaZ=as.numeric(args[5]), embedding_orig=args[6], embedding_obs=args[7], confounder=args[8],
                                                                          response=args[9])
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
                                                        if(args[8]=='AS'){
                                                          model_formula_YZ <- "V1~1+s(V3, by=V2)"
                                                        }else if(args[8]=='genes10'){
                                                          add_confounders <- paste('+',paste('s(V', 4:(ncol(Y)+ncol(Z)), ')', collapse='+', sep=""), sep="")
                                                          model_formula_YZ <- paste('V1~1+s(V3, by=V2)', add_confounders, sep="")
                                                        }
                                                        cit_params <- list(cit='cpt_kpc', params_cit=list(k=k, Knn = as.numeric(args[12]), model.formula.YZ=model_formula_YZ))
                                                     }else if(args[10]=='FCIT'){
                                                        cit_params <- list(cit='fcit')
                                                     }else if(args[10]=='WALD'){
                                                        cit_params <- list(cit='wald', params_cit=NULL)
                                                     }

                                                     tmp <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                                                  cit_with_parameters = cit_params)
                                                     res[idx_beta2] <- tmp$p
                                                     runtime[idx_beta2] <- tmp$runtime
                                                   }
                                                   p_time <- cbind(res, runtime)
                                                 }
                                                 p_time
                                               }
parallel::stopCluster(cl)

if(grepl('/CI',args[1],fixed=TRUE)){
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
}else if(grepl('/No_CI',args[1],fixed=TRUE)){
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
