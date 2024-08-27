library(ranger)
library(tidymodels)

yxz <- as.data.frame(do.call(cbind, list(Y,X,Z)))
colnames(yxz) <- paste0("V", 1:ncol(yxz))

data_split <- initial_split(yxz, prop = 0.8, strata = "V1")
yxz_train <- training(data_split)
yxz_test <- testing(data_split)

