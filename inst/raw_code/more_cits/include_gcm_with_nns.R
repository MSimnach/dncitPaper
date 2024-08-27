library(torch)
library(magrittr)
library(luz)
library(rsample)

#### Data preparation
## Define split indices
train_ratio <- 0.7
val_ratio <- 0.15
test_ratio <- 0.15

n <- nrow(Z)
batch_size <- n

n_train <- floor(train_ratio * n)
n_val <- floor(val_ratio * n)
n_test <- n - n_train - n_val

indices <- sample(seq_len(n))

train_indices <- indices[1:n_train]
val_indices <- indices[(n_train + 1):(n_train + n_val)]
test_indices <- indices[(n_train + n_val + 1):n]

# Split the data
Z_train <- torch_tensor(Z[train_indices, ,drop=FALSE])
Y_train <- torch_tensor(Y[train_indices,,drop=FALSE])
X_train <- torch_tensor(X[train_indices,])

Z_val <- torch_tensor(Z[val_indices,,drop=FALSE])
Y_val <- torch_tensor(Y[val_indices,,drop=FALSE])
X_val <- torch_tensor(X[val_indices,])

Z_test <- torch_tensor(Z[-c(train_indices, val_indices),,drop=FALSE])
Y_test <- torch_tensor(Y[-c(train_indices, val_indices),,drop=FALSE])
X_test <- torch_tensor(X[-c(train_indices, val_indices),])


#### Hyperparameters
grid <- expand.grid(
  lr = c(0.001, 0.003, 0.01),    # Learning rates to try
  dropout_rate = c(0.2, 0.3, 0.4) # Dropout rates to try
)
dropout_rate <- 0.3
lr = 0.003
num_epochs <- 100


##### Regress Y on Z
# Define ZY dataset
dataset_ZY <- dataset(
  initialize = function(Z, Y) {
    self$Z <- Z
    self$Y <- Y
  },
  .getitem = function(i) {
    list(Z = self$Z[i, ], Y = self$Y[i])
  },
  .length = function() {
    self$Z$size(1)
  }
)

train_ds_ZY <- dataset_ZY(Z_train, Y_train)
train_dl_ZY <- dataloader(train_ds_ZY, batch_size = batch_size, shuffle = TRUE)
val_ds_ZY <- dataset_ZY(Z_val, Y_val)
val_dl_ZY <- dataloader(val_ds_ZY, batch_size = batch_size)
test_ds_ZY <- dataset_ZY(Z_test, Y_test)
test_dl_ZY <- dataloader(test_ds_ZY, batch_size = batch_size)

### Model predicting Y from Z
net_Y <- nn_module(
  initialize = function() {
    self$linear1 <- nn_linear(ncol(Z), 10) # Adjust input/output sizes as needed
    self$linear2 <- nn_linear(10, 1)

    self$bn1 <- nn_batch_norm1d(num_features = 10)
  },
  forward = function(Z) {
    Z %>%
      self$linear1() %>%
      nnf_relu() %>%
      self$bn1() %>%
      nnf_dropout(p = dropout_rate) %>%

      self$linear2()
  }
)

fitted_Y <- net_Y %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true),
    optimizer = optim_adam,
    metrics = list(luz_metric_mse(), luz_metric_mae())
  ) %>%
  set_opt_hparams(lr = lr) %>%
  fit(train_dl_ZY,
      epochs=num_epochs,
      valid_data  = val_dl_ZY,
      callbacks = list(
        luz_callback_early_stopping(patience = 10)
      )
  )

predictions_Y <- predict(fitted_Y, test_dl_ZY)

##### Regress X on Z
### Define ZX dataset
dataset_ZX <- dataset(
  initialize = function(Z, X) {
    self$Z <- Z
    self$X <- X
  },
  .getitem = function(i) {
    list(Z = self$Z[i, ], X = self$X[i, ])
  },
  .length = function() {
    self$Z$size(1)
  }
)

train_ds_ZX <- dataset_ZX(Z_train, X_train)
train_dl_ZX <- dataloader(train_ds_ZX, batch_size = batch_size, shuffle = TRUE)
val_ds_ZX <- dataset_ZY(Z_val, X_val)
val_dl_ZX <- dataloader(val_ds_ZX, batch_size = batch_size)
test_ds_ZX <- dataset_ZX(Z_test, X_test)
test_dl_ZX <- dataloader(test_ds_ZX, batch_size = batch_size)

### Model predicting X from Z
net_X <- nn_module(
  initialize = function() {
    self$linear1 <- nn_linear(in_features = ncol(Z), out_features = 10)
    self$linear2 <- nn_linear(in_features = 10, out_features =100)
    self$linear3 <- nn_linear(in_features = 100, out_features =139)

    self$bn1 <- nn_batch_norm1d(num_features = 10)
    self$bn2 <- nn_batch_norm1d(num_features = 100)
  },
  forward = function(Z) {
    Z %>%
      self$linear1() %>%
      nnf_relu() %>%
      self$bn1() %>%
      nnf_dropout(p = dropout_rate) %>%

      self$linear2() %>%
      nnf_relu() %>%
      self$bn2() %>%
      nnf_dropout(p = dropout_rate) %>%

      self$linear3()
  }
)

### model fit
#find lr
model_X <- net_X %>% setup(
  loss = function(x_hat, x_true) nnf_mse_loss(x_hat, x_true),
  optimizer = torch::optim_adam
)

records_X <- lr_finder(
  object = model_X,
  data = train_ds_ZX,
  verbose = FALSE,
  dataloader_options = list(batch_size = batch_size),
  start_lr = 1e-3, # the smallest value that will be tried
  end_lr = 1 # the largest value to be experimented with
)

fitted_X <- model_X %>%
  setup(
    loss = function(x_hat, x_true) nnf_mse_loss(x_hat, x_true),
    optimizer = optim_adam,
    metrics = list(luz_metric_mse(), luz_metric_mae())
  ) %>%
  set_opt_hparams(lr = lr) %>%
  fit(train_dl_ZX,
      epochs=num_epochs,
      valid_data  = val_dl_ZX,
      callbacks = list(
        luz_callback_early_stopping(patience = 10)
      )
  )


#### GCM
# Residuals
predictions_Y <- as_array(predict(fitted_Y, test_dl_ZY))
true_Y <- as.array(test_dl_ZY$dataset$Y)
res_Y <- true_Y - predictions_Y

predictions_X <- as_array(predict(fitted_X, test_dl_ZX))
true_X <- as.array(test_dl_ZX$dataset$X)
res_X <- true_X - predictions_X
library(GeneralisedCovarianceMeasure)
gcm.test(X,Y,Z, resid.YonZ = res_Y, resid.XonZ = res_X)

gcm.test(X,Y,Z)

#compare to lm
lm_YZ <- lm(as.matrix(Y)~as.matrix(Z))
mse(lm_YZ$fitted.values,Y)
sum(lm_YZ$residuals^2)
