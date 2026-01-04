## Surface spline simulation in appendix for fine-tuned BMARS
## Designed for SLURM array jobs:
## one array task runs one or more experiment replicates

suppressPackageStartupMessages({
  library(ABMRS)
  library(BASS)
  library(GauPro)
  library(mgcv)
})

# ---------- helpers ----------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  return(getwd())
}

get_job_id <- function(args) {
  env_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
  if (!is.na(env_id) && nzchar(env_id)) return(as.integer(env_id))
  if (length(args) >= 1) return(as.integer(args[1]))
  return(1L)
}

get_rep_count <- function(args) {
  if (length(args) >= 2) return(max(1L, as.integer(args[2])))
  return(1L)
}

# ---------- config ----------
script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
job_id <- get_job_id(args)
rep_count <- get_rep_count(args)

set.seed(1234 + job_id)

out_dir <- file.path(repo_root, "output", "surface_fitting_appendix_bmars")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

job_tag <- sprintf("job-%04d", job_id)
raw_mse_out <- file.path(out_dir, sprintf("surface_mse_raw_%s.rds", job_tag))
raw_censor_mse_out <- file.path(out_dir, sprintf("surface_mse_censor_raw_%s.rds", job_tag))
raw_cov_out <- file.path(out_dir, sprintf("surface_coverage_raw_%s.rds", job_tag))
raw_bandwidth_out <- file.path(out_dir, sprintf("surface_bandwidth_raw_%s.rds", job_tag))

# ---------- functions ----------
g1 = function(x) {
  10.391 * ((x[,1]-0.4)*(x[,2]-0.6)) + 0.36
}

g2 = function(x) {
  r2 = (x[,1]-0.5)^2 + (x[,2]-0.5)^2
  24.342 * r2 * (0.75 - r2)
}

g3 = function(x) {
  x1_hat = x[,1] - 0.5
  x2_hat = x[,2] - 0.5
  42.659 * (0.1 + x1_hat*(0.05+x1_hat^4-10*x1_hat^2*x2_hat^2+5*x2_hat^4))
}

g4 = function(x) {
  1.3356 * (1.5*(1-x[,1]) + exp(2*x[,1]-1)*sin(3*pi*(x[,1]-0.6)^2) + exp(3*(x[,2]-0.5))*sin(4*pi*(x[,2]-0.9)^2))
}

g5 = function(x) {
  1.9 * (1.35 + exp(x[,1])*sin(13*(x[,1]-0.6)^2)*exp(-x[,2])*sin(7*x[,2]))
}

gs = list(g1,g2,g3,g4,g5)
noises_bass = c(0.25,0.22,0.22,0.25,0.25)

# ---------- simulation setup ----------
m_trains <- c(250, 1000)
m_tests <- c(250, 1000)
p <- 0.025

if (length(m_tests) != length(m_trains)) {
  stop("m_tests must have same length as m_trains")
}

level <- 0.95
upper_level <- (1 + level) / 2
lower_level <- (1 - level) / 2

cases <- c("Case 1", "Case 2", "Case 3", "Case 4", "Case 5")
methods <- c(
  "BASS (degree=1)",
  "BASS (degree=2)",
  "BASS (degree=3)"
)

row_labels <- paste(rep(cases, times = length(m_trains)), rep(m_trains, each = length(cases)), sep = "-")

mse_array <- array(
  NA_real_,
  dim = c(length(row_labels), length(methods), rep_count),
  dimnames = list(row_labels, methods, NULL)
)
censor_mse_array <- array(
  NA_real_,
  dim = c(length(row_labels), length(methods), rep_count),
  dimnames = list(row_labels, methods, NULL)
)
cov_array <- array(
  NA_real_,
  dim = c(length(row_labels), length(methods), rep_count),
  dimnames = list(row_labels, methods, NULL)
)
bandwidth_array <- array(
  NA_real_,
  dim = c(length(row_labels), length(methods), rep_count),
  dimnames = list(row_labels, methods, NULL)
)

time_start <- Sys.time()
message("Starting job ", job_id, " with ", rep_count, " replicate(s)")

for (k in seq_len(rep_count)) {
  message("Replicate ", k)

  for (s in seq_along(m_trains)) {
    m_train <- m_trains[s]
    m_test <- m_tests[s]
    row_offset <- (s - 1L) * length(cases)
    message("Train size ", m_train, ", test size ", m_test)

    x_train <- cbind(runif(m_train, 0, 1), runif(m_train, 0, 1))
    x_test <- cbind(runif(m_test, 0, 1), runif(m_test, 0, 1))
    xmin <- c(0, 0)
    xmax <- c(1, 1)

    for (i in seq_along(gs)) {
      row_idx <- row_offset + i
      message("  Running ", cases[i])
      f <- gs[[i]]
      noise <- noises_bass[[i]]
      
      y_train <- f(x_train)
      y_obs <- y_train + rnorm(m_train, 0, noise)
      y_test <- f(x_test)

      censor_lo <- max(1L, round(p * m_test))
      censor_hi <- min(m_test, round((1 - p) * m_test))


      # BASS (degree=1)
      bass_model <- bass(xx = x_train, y = y_obs, degree = 1, nmcmc = 10000, nburn = 9000, verbose = FALSE)
      pred_bass <- predict(bass_model, x_test)
      y_bass <- colMeans(pred_bass)
      mse_array[row_idx, "BASS (degree=1)", k] <- mean((y_bass - y_test)^2)
      censor_mse_array[row_idx, "BASS (degree=1)", k] <- mean(sort((y_bass - y_test)^2)[censor_lo:censor_hi])
      lower_bound_bass <- apply(pred_bass, 2, quantile, probs = lower_level)
      upper_bound_bass <- apply(pred_bass, 2, quantile, probs = upper_level)
      cov_array[row_idx, "BASS (degree=1)", k] <- mean((y_test >= lower_bound_bass) & (y_test <= upper_bound_bass))
      bandwidth_array[row_idx, "BASS (degree=1)", k] <- mean(upper_bound_bass - lower_bound_bass)
      
      # BASS (degree=2)
      bass_model <- bass(xx = x_train, y = y_obs, degree = 2, nmcmc = 10000, nburn = 9000, verbose = FALSE)
      pred_bass <- predict(bass_model, x_test)
      y_bass <- colMeans(pred_bass)
      mse_array[row_idx, "BASS (degree=2)", k] <- mean((y_bass - y_test)^2)
      censor_mse_array[row_idx, "BASS (degree=2)", k] <- mean(sort((y_bass - y_test)^2)[censor_lo:censor_hi])
      lower_bound_bass <- apply(pred_bass, 2, quantile, probs = lower_level)
      upper_bound_bass <- apply(pred_bass, 2, quantile, probs = upper_level)
      cov_array[row_idx, "BASS (degree=2)", k] <- mean((y_test >= lower_bound_bass) & (y_test <= upper_bound_bass))
      bandwidth_array[row_idx, "BASS (degree=2)", k] <- mean(upper_bound_bass - lower_bound_bass)
      
      # BASS (degree=3)
      bass_model <- bass(xx = x_train, y = y_obs, degree = 3, nmcmc = 10000, nburn = 9000, verbose = FALSE)
      pred_bass <- predict(bass_model, x_test)
      y_bass <- colMeans(pred_bass)
      mse_array[row_idx, "BASS (degree=3)", k] <- mean((y_bass - y_test)^2)
      censor_mse_array[row_idx, "BASS (degree=3)", k] <- mean(sort((y_bass - y_test)^2)[censor_lo:censor_hi])
      lower_bound_bass <- apply(pred_bass, 2, quantile, probs = lower_level)
      upper_bound_bass <- apply(pred_bass, 2, quantile, probs = upper_level)
      cov_array[row_idx, "BASS (degree=3)", k] <- mean((y_test >= lower_bound_bass) & (y_test <= upper_bound_bass))
      bandwidth_array[row_idx, "BASS (degree=3)", k] <- mean(upper_bound_bass - lower_bound_bass)
    }
  }
}

time_end <- Sys.time()
message("Elapsed: ", format(time_end - time_start))

saveRDS(mse_array, file = raw_mse_out)
saveRDS(censor_mse_array, file = raw_censor_mse_out)
saveRDS(cov_array, file = raw_cov_out)
saveRDS(bandwidth_array, file = raw_bandwidth_out)

message("Wrote raw MSE array to ", raw_mse_out)
message("Wrote raw censored MSE array to ", raw_censor_mse_out)
message("Wrote raw coverage array to ", raw_cov_out)
message("Wrote raw bandwidth array to ", raw_bandwidth_out)

