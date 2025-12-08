
## --------------------------------------------------------------------------------------
# load package
library(oro.nifti)
library(neurobase)
library(EBARS)
library(mgcv)

# file path
data_path = "/home/hejunhui/R/ABMRS"

set.seed(1234)

# read nii
contrast_map = readnii(file.path(data_path, "Contrasts", "168745_2bk-0bk.nii"))
mni152 = readnii(file.path(data_path, "Contrasts", "MNI152_T1_2mm.nii"))

slice_idxs = c(25, 35, 50, 65)

activation_ebars = mni152@.Data
activation_thres = mni152@.Data

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No SLURM_ARRAY_TASK_ID provided. Please pass it as an argument.")
}
task_id <- as.numeric(args[1])
cat("The SLURM_ARRAY_TASK_ID is:", task_id, "\n")

# run simulations
idx = slice_idxs[task_id]

contrast_slice = contrast_map[,,idx]

slice_array = cbind(c(row(contrast_slice)), c(col(contrast_slice)), c(contrast_slice))
slice_array = slice_array[!is.na(slice_array[,3]),]

# normalize locations
x = as.matrix(slice_array[,1:2])
for (i in 1:2) {
  x[,i] = (x[,i] - min(x[,i])) / (max(x[,i]) - min(x[,i]))
}
y = slice_array[,3]

## --------------------------------------------------------------------------------------
# apply EBARS
model_bebars = binebars(x, y, k_1=8, k_2=8, n_1=1000, n_2=1000)

time_start = Sys.time()
model_bebars$mcmc(burns=1000, steps=1000)
time_end = Sys.time()
print(time_end - time_start)


## --------------------------------------------------------------------------------------
# fit a spline model
knots = model_bebars$knots()
z = tensor_spline(x, xi_1=knots$xi_1, xi_2=knots$xi_2, intercept_1 = TRUE, intercept_2 = TRUE)
reg_df = as.data.frame(z)
reg_df$y = y
lm_model = lm(y ~ ., data=reg_df)

# compute predictive intervals
pred_int = as.data.frame(predict(lm_model, reg_df, interval="prediction", level=0.999))

# calculate the activation region
region_pred = (pred_int$lwr>0) | (pred_int$upr<0)
region_act = array(FALSE, dim = dim(mni152)[c(1,2)])
for (i in 1:length(region_pred)) {
  region_act[slice_array[i,1], slice_array[i,2]] = region_pred[i]
}

activation_ebars[,,idx] = region_act
pdf(file.path(data_path, "figures", paste0("activation_ebars_", idx, ".pdf")), width = 4, height = 4)
slice_overlay(mni152, activation_ebars, z = idx, plane = "axial", bg = "transparent", col.y = "red")
dev.off()

region_thres = abs(contrast_slice) > quantile(abs(y), 1 - sum(region_pred)/length(y))
region_thres[is.na(region_thres)] = FALSE
activation_thres[,,idx] = region_thres
pdf(file.path(data_path, "figures", paste0("activation_thres_", idx, ".pdf")), width = 4, height = 4)
slice_overlay(mni152, activation_thres, z = idx, plane = "axial", bg = "transparent", col.y = "red")
dev.off()
