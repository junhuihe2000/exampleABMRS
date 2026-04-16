
## --------------------------------------------------------------------------------------
# load package
library(oro.nifti)
library(neurobase)
library(EBARS)

# file path
data_path = "/scratch/jiankang_root/jiankang0/hejunhui/WM_contrasts"
output_path = "/home/hejunhui/R/ABMRS/fmri"

set.seed(1234)

aal_codes = c(6201, 2211, 6101, 6102, 2212, 2301, 2001, 2201, 6202, 2111)

n_regions = 10

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No SLURM_ARRAY_TASK_ID provided. Please pass it as an argument.")
}
task_id <- as.numeric(args[1])
cat("The SLURM_ARRAY_TASK_ID is:", task_id, "\n")

# read data
index_subjects = read.csv(file.path(data_path, "index2bk-0bk.csv"))
subject_name = index_subjects$subject[task_id]

# read nii
contrast_map = readnii(file.path(data_path, "Contrasts", paste0(subject_name, "_2bk-0bk.nii")))
aal_mni = readnii(file.path(data_path, "AAL_MNI_2mm.nii"))

# preprocess nii
contrast_map_array = array(NA, c(length(contrast_map), 5))
colnames(contrast_map_array) = c("x1", "x2", "x3", "y", "code")
idx = 0
for (i in 1:dim(contrast_map)[1]) {
  for (j in 1:dim(contrast_map)[2]) {
    for (k in 1:dim(contrast_map)[3]) {
      idx = idx + 1
      contrast_map_array[idx,] = c(i, j, k, contrast_map[i, j, k], aal_mni[i, j, k])
    }
  }
}

contrast_map_array = contrast_map_array[!is.na(contrast_map_array[,4]),]

peak_results = list()

n_regions = 10

# run simulations
for (id in 1:n_regions) {
idx = contrast_map_array[,5] == aal_codes[id]
x = contrast_map_array[idx, 1:3]
# normalize locations
for (i in 1:3) {
  x[,i] = (x[,i]- min(x[,i])) / (max(x[,i]) - min(x[,i]))
}
y = contrast_map_array[idx, 4]

## --------------------------------------------------------------------------------------
# apply EBARS
model_ebars = triebars(x, y, k_1 = 5, k_2 = 5, k_3 = 5, n_1 = 100, n_2 = 100, n_3 = 100, degree_1 = 1, degree_2 = 1, degree_3 = 1)
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

  peak_df[i,4] = all(y_hat[1]>=y_hat[-1])
  peak_df[i,5] = all(y_hat[1]<=y_hat[-1])
}

for (j in 1:3) {
  peak_df[,j] = peak_df[,j] * (max(contrast_map_array[idx,j]) - min(contrast_map_array[idx,j])) + min(contrast_map_array[idx,j])
}

peak_results[[id]] = peak_df
}

save(peak_results, file = file.path(output_path, "output", paste0("peak_", subject_name, ".RData")))
