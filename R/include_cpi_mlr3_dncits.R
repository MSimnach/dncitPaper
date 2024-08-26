library(mlr3)
library(mlr3learners)
library(data.table)
library(cpi)
yxz <- do.call(cbind, list(Y,X,Z))
colnames(yxz) <- paste0("V", 1:ncol(yxz))
tsk_yxz <-  as_task_regr(yxz, target = "V1")
lrn_ranger <- lrn("regr.ranger")

measure = msrs("regr.mse")
resampling = rsmp("cv", folds = 5)
#resample_result <- resample(tsk_yxz, lrn_ranger, cv3)
#performance <- resample_result$aggregate(measure)
params_cit <- list(task = tsk_yxz, learner = lrn_ranger,
                   resampling = resampling, groups = list(X = 2:(ncol(X)+1)))
updated_parameters <- update_params_cits(cpi::cpi, X,Y,Z,params_cit)[1:13]
do.call(cpi::cpi, updated_parameters)

