
## --------------------------------------------------------------------------------------
# load package
library(oro.nifti)
library(neurobase)
library(EBARS)

# file path
data_path = "/home/hejunhui/R/ABMRS/fmri"

set.seed(1234)

aal_codes = c(6201, 6101, 6102, 2001, 2201)

load(file = file.path(data_path, "data/2bk-0bk_contrast.RData"))

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No SLURM_ARRAY_TASK_ID provided. Please pass it as an argument.")
}
task_id <- as.numeric(args[1])
cat("The SLURM_ARRAY_TASK_ID is:", task_id, "\n")

# run simulations
idx = contrast_map_array[,5] == aal_codes[task_id]
x = contrast_map_array[idx, 1:3]
# normalize locations
for (i in 1:3) {
  x[,i] = (x[,i]- min(x[,i])) / (max(x[,i]) - min(x[,i]))
}
y = contrast_map_array[idx, 4]

## --------------------------------------------------------------------------------------
# apply EBARS
model_ebars = triebars(x, y, k_1 = 10, k_2 = 10, k_3 = 10, n_1 = 1000, n_2 = 1000, n_3 = 1000, degree_1 = 1, degree_2 = 1, degree_3 = 1)
time_start = Sys.time()
model_ebars$mcmc(burns=500, steps=500)
time_end = Sys.time()
print(time_end - time_start)


## ----------------------------------------------------------------------------------------------
# fit a spline model
knots = model_ebars$knots()
z = tri_tensor_spline(x, xi_1=knots$xi_1, xi_2=knots$xi_2, xi_3=knots$xi_3, degree_1 = 1, degree_2 = 1, degree_3 = 1, intercept_1 = TRUE, intercept_2 = TRUE, intercept_3 = TRUE)
reg_df = as.data.frame(z)
reg_df$y = y
lm_model = lm(y ~ ., data=reg_df)

peak_df = data.frame(expand.grid(knots$xi_1, knots$xi_2, knots$xi_3))
colnames(peak_df) = c("x1", "x2", "x3")
peak_df$y_hat = rep(NA, nrow(peak_df))
peak_df$max = rep(FALSE, nrow(peak_df))
peak_df$min = rep(FALSE, nrow(peak_df))

# determine the peak points
for (i in 1:nrow(peak_df)) {
  x_cand = unlist(peak_df[i,1:3])
  x_nbs = expand.grid(seq(x_cand[1]-0.02, x_cand[1]+0.02, 0.01), seq(x_cand[2]-0.02, x_cand[2]+0.02, 0.01), seq(x_cand[3]-0.02, x_cand[3]+0.02, 0.01))
  z = tri_tensor_spline(as.matrix(rbind(x_cand, x_nbs)), xi_1=knots$xi_1, xi_2=knots$xi_2, xi_3=knots$xi_3, degree_1 = 1, degree_2 = 1, degree_3 = 1, intercept_1 = TRUE, intercept_2 = TRUE, intercept_3 = TRUE)
  z = as.data.frame(z)

  # make predictions
  y_hat = predict(lm_model, z)

  peak_df$y_hat[i] = y_hat[1]
  peak_df$max[i] = y_hat[1] == max(y_hat)
  peak_df$min[i] = y_hat[1] == min(y_hat)
}

for (j in 1:3) {
  peak_df[,j] = peak_df[,j] * (max(contrast_map_array[idx,j]) - min(contrast_map_array[idx,j])) + min(contrast_map_array[idx,j])
}

save(peak_df, file = file.path(data_path, "output", paste0("region_", aal_codes[task_id], ".RData")))
