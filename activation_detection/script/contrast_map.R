
## ----------------------------------------------------------------------------------------------
# load package
library(oro.nifti)
library(neurobase)
library(EBARS)
library(mgcv)

# file path
data_path = "/home/hejunhui/R/ABMRS"

set.seed(1234)

slice_idxs = c(25, 35, 50, 65)

# test mse
test_mse = array(NA, dim = c(3,4))

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No SLURM_ARRAY_TASK_ID provided. Please pass it as an argument.")
}
task_id <- as.numeric(args[1])
cat("The SLURM_ARRAY_TASK_ID is:", task_id, "\n")

# run simulations
for (i in 1:4) {
  idx = slice_idxs[i]
  load(file = file.path(data_path, "data", paste0("slice_", idx, ".RData")))
  batch = round(nrow(x) / 10)
  if(task_id<10) {test_idx = c(((task_id-1)*batch+1):(task_id*batch))}
  else {test_idx = c(((task_id-1)*batch+1):(nrow(x)))}
  train_idx = c(1:nrow(x))[-test_idx]

  ## ----------------------------------------------------------------------------------------------
  # apply EBARS
  model_bebars = binebars(x[train_idx,], y[train_idx], k_1=8, k_2=8, n_1=1000, n_2=1000)

  time_start = Sys.time()
  model_bebars$mcmc(burns=1000, steps=1000)
  time_end = Sys.time()
  print(time_end - time_start)

  ## ----------------------------------------------------------------------------------------------
  # fit a spline model
  knots = model_bebars$knots()
  z = tensor_spline(x, xi_1=knots$xi_1, xi_2=knots$xi_2, intercept_1 = TRUE, intercept_2 = TRUE)
  reg_df = as.data.frame(z)
  reg_df$y = y
  lm_model = lm(y ~ ., data=reg_df, subset = train_idx)

  # make predictions
  y_hat = predict(lm_model, reg_df)[-train_idx]

  # compute test mse
  test_mse[1,i] = mean((y_hat-y[-train_idx])^2)
  # Thin Plate Spline
  x_train_df = data.frame(x1 = x[train_idx,1], x2 = x[train_idx,2])
  x_new_df = data.frame(x1 = x[-train_idx,1], x2 = x[-train_idx,2])
  model_tps = gam(y[train_idx]~s(x1, x2, k = 144), data = x_train_df)
  y_tps = predict(model_tps,x_new_df)
  test_mse[2,i] = mean((y[-train_idx]-y_tps)^2)
  # Tensor Product Spline
  model_tep = gam(y[train_idx]~te(x1, x2, k = 12), data = x_train_df)
  y_tep = predict(model_tep, x_new_df)
  test_mse[3,i] = mean((y[-train_idx]-y_tep)^2)
}



## ----------------------------------------------------------------------------------------------
save(test_mse, file = file.path(data_path, "output", paste0("test_mse_", task_id, ".RData")))

