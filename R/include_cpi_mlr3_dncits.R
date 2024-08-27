##### Implementation tests for CPI based CIT and mlr3 prediction based CIT
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3tuningspaces)
library(data.table)
library(cpi)

#### 1) CPI
yxz <- do.call(cbind, list(Y,X,Z))
colnames(yxz) <- paste0("V", 1:ncol(yxz))
tsk_yxz <-  as_task_regr(yxz, target = "V1")
lrn_ranger <- lrn("regr.ranger")

measure = msr("regr.mse")
resampling = rsmp("cv", folds = 5)
#resample_result <- resample(tsk_yxz, lrn_ranger, cv3)
#performance <- resample_result$aggregate(measure)
params_cit <- list(task = tsk_yxz, learner = lrn_ranger,
                   resampling = resampling, groups = list(X = 2:(ncol(X)+1)))
updated_parameters <- update_params_cits(cpi::cpi, X,Y,Z,params_cit)[1:13]
cpitest <- do.call(cpi::cpi, updated_parameters)

#### Prediction based CIT with mlr3 and random forest
## Auto tuner for hyperparameter optimization
## Nested resampling for performance estimation
# termination time for auto tuner
term_time_at <- round(exp(4.5)/3)

# y on x,z
tsk_yxz <-  as_task_regr(yxz, target = "V1")
measure <- msr("regr.mse")
rsmp_holdout <- rsmp("holdout")
rsmp_cv3 <- rsmp("cv", folds = 3)
learner <- mlr3tuningspaces::lts(lrn("regr.ranger"))

at_yxz <- auto_tuner(
  tuner = tnr("random_search"),
  learner = learner,
  resampling = rsmp_holdout,
  measure = measure,
  term_time = term_time_at
)

rr_yxz <- resample(tsk_yxz, at_yxz, rsmp_cv3)


gen_error_estimate_yxz <- rr_yxz$aggregate()

# y on z
yz <- do.call(cbind, list(Y,Z))
colnames(yz) <- paste0("V", 1:ncol(yz))
tsk_yz <-  as_task_regr(yz, target = "V1")

at_yz <- auto_tuner(
  tuner = tnr("random_search"),
  learner = learner,
  resampling = rsmp_holdout,
  measure = measure,
  term_time = term_time_at
)

rr_yz <- resample(tsk_yxz, at_yz, rsmp_cv3)

gen_error_estimate_yz <- rr_yz$aggregate()

## compare residuals
preds_yxz <- rr_yxz$prediction()
res_yxz <- (preds_yxz$truth - preds_yxz$response)^2
preds_yz <- rr_yz$prediction()
res_yz <- (preds_yz$truth - preds_yz$response)^2
t_test <- stats::t.test(res_yxz, res_yz, paired = TRUE, alternative = "less")
wilcox_test <- stats::wilcox.test(res_yxz, res_yz, paired = TRUE, alternative = "less")
t_test$p.value
wilcox_test$p.value

