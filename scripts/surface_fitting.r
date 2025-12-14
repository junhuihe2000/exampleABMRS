## Surface spline simulation
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

out_dir <- file.path(repo_root, "output", "surface_fitting")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

job_tag <- sprintf("job-%04d", job_id)
raw_mse_out <- file.path(out_dir, sprintf("surface_mse_raw_%s.rds", job_tag))
raw_censor_mse_out <- file.path(out_dir, sprintf("surface_mse_censor_raw_%s.rds", job_tag))
raw_cov_out <- file.path(out_dir, sprintf("surface_coverage_raw_%s.rds", job_tag))
raw_bandwidth_out <- file.path(out_dir, sprintf("surface_bandwidth_raw_%s.rds", job_tag))

# ---------- functions ----------
f1 <- function(x) {
  -40 * (exp(-(x[, 1] - 0.3)^2 / 0.1 - (x[, 2] - 0.3)^2 / 0.1) -
    exp(-(x[, 1] - 0.7)^2 / 0.1 - (x[, 2] - 0.7)^2 / 0.1))
}

f2 <- function(x) {
  xi_1 <- c(0.2, 0.5)
  xi_2 <- c(0.4, 0.4, 0.4, 0.4, 0.7)
  beta <- c(-1.21, 0.28, 1.08, -2.35, 0.43, 0.51, -0.57, -0.55, -0.56, -0.89, -0.48, -1.00, -0.78, 0.06, 0.96, -0.11, -0.51, -0.91, -0.84, 2.42, 0.13, -0.49, -0.44, 0.46, -0.69, -1.45, 0.57, -1.02, -0.02, -0.94, 1.10, -0.48, -0.71, -0.50, -1.63, -1.17, -2.18, -1.34, -0.29, -0.47)
  B <- tensor_spline(x, list(xi_1, xi_2), c(3, 3), c(FALSE, FALSE))
  c(20 * B %*% cbind(beta))
}

f3 <- function(x) {
  xi_1 <- c(0.2, 0.2, 0.2, 0.2, 0.6)
  xi_2 <- c(0.5, 0.5, 0.5, 0.5, 0.7)
  beta <- c(0.83, -0.28, -0.36, 0.09, 2.25, 0.83, 1.31, 2.50, 1.17, -0.43, -1.00, -1.11, -0.06, 1.17, 1.05, 0.06, -0.74, 0.93, 1.67, 0.56, -0.75, 1.26, 0.04, 0.19, 0.46, -0.43, 0.02, 0.70, 0.97, -0.62, -0.86, 0.07, -1.05, -2.75, -1.13, -0.86, 1.56, 1.02, 1.04, -1.12, -1.07, 0.97, 0.17, -0.90, 0.16, -0.50, -0.97, -0.11, 1.09, -1.21, -1.77, -0.49, 0.32, 1.46, 1.54, -0.34, -1.08, -1.49, -0.25, -0.12, -0.65, 0.31, 0.12, -0.84)
  B <- tensor_spline(x, list(xi_1, xi_2), c(3, 3), c(FALSE, FALSE))
  c(20 * B %*% cbind(beta))
}

fs <- list(f1, f2, f3)
noises <- c(0.8, 2, 2)

# ---------- simulation setup ----------
m_trains <- c(500, 1000)
m_tests <- c(500, 1000)
p <- 0.025

if (length(m_tests) != length(m_trains)) {
  stop("m_tests must have same length as m_trains")
}

level <- 0.95
upper_level <- (1 + level) / 2
lower_level <- (1 - level) / 2

cases <- c("Case 1", "Case 2", "Case 3")
methods <- c(
  "ABMRS (gamma=0.8)",
  "ABMRS (gamma=1)",
  "ABMRS (gamma=1.2)",
  "Thin Plate Spline",
  "Tensor Product Spline",
  "Gaussian Process",
  "BASS"
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

    for (i in seq_along(fs)) {
      row_idx <- row_offset + i
      message("  Running ", cases[i])
      f <- fs[[i]]
      noise <- noises[[i]]

      y_train <- f(x_train)
      y_obs <- y_train + rnorm(m_train, 0, noise)
      y_test <- f(x_test)

      censor_lo <- max(1L, round(p * m_test))
      censor_hi <- min(m_test, round((1 - p) * m_test))

      # ABMRS with gamma = 0.8
      mebars_1 <- mebars(x_train, y_obs, gamma = 0.8, times = c(2, 2))
      mebars_1$rjmcmc(burns = 20000, steps = 20000)
      pred_1 <- mebars_1$predict(x_test)
      y_mebars_1 <- rowMeans(pred_1)
      mse_array[row_idx, "ABMRS (gamma=0.8)", k] <- mean((y_mebars_1 - y_test)^2)
      censor_mse_array[row_idx, "ABMRS (gamma=0.8)", k] <- mean(sort((y_mebars_1 - y_test)^2)[censor_lo:censor_hi])
      lower_bound_1 <- apply(pred_1, 1, quantile, probs = lower_level)
      upper_bound_1 <- apply(pred_1, 1, quantile, probs = upper_level)
      cov_array[row_idx, "ABMRS (gamma=0.8)", k] <- mean((y_test >= lower_bound_1) & (y_test <= upper_bound_1))
      bandwidth_array[row_idx, "ABMRS (gamma=0.8)", k] <- mean(upper_bound_1 - lower_bound_1)

      # ABMRS with gamma = 1
      mebars_2 <- mebars(x_train, y_obs, gamma = 1, times = c(2, 2))
      mebars_2$rjmcmc(burns = 20000, steps = 20000)
      pred_2 <- mebars_2$predict(x_test)
      y_mebars_2 <- rowMeans(pred_2)
      mse_array[row_idx, "ABMRS (gamma=1)", k] <- mean((y_mebars_2 - y_test)^2)
      censor_mse_array[row_idx, "ABMRS (gamma=1)", k] <- mean(sort((y_mebars_2 - y_test)^2)[censor_lo:censor_hi])
      lower_bound_2 <- apply(pred_2, 1, quantile, probs = lower_level)
      upper_bound_2 <- apply(pred_2, 1, quantile, probs = upper_level)
      cov_array[row_idx, "ABMRS (gamma=1)", k] <- mean((y_test >= lower_bound_2) & (y_test <= upper_bound_2))
      bandwidth_array[row_idx, "ABMRS (gamma=1)", k] <- mean(upper_bound_2 - lower_bound_2)

      # ABMRS with gamma = 1.2
      mebars_3 <- mebars(x_train, y_obs, gamma = 1.2, times = c(2, 2))
      mebars_3$rjmcmc(burns = 20000, steps = 20000)
      pred_3 <- mebars_3$predict(x_test)
      y_mebars_3 <- rowMeans(pred_3)
      mse_array[row_idx, "ABMRS (gamma=1.2)", k] <- mean((y_mebars_3 - y_test)^2)
      censor_mse_array[row_idx, "ABMRS (gamma=1.2)", k] <- mean(sort((y_mebars_3 - y_test)^2)[censor_lo:censor_hi])
      lower_bound_3 <- apply(pred_3, 1, quantile, probs = lower_level)
      upper_bound_3 <- apply(pred_3, 1, quantile, probs = upper_level)
      cov_array[row_idx, "ABMRS (gamma=1.2)", k] <- mean((y_test >= lower_bound_3) & (y_test <= upper_bound_3))
      bandwidth_array[row_idx, "ABMRS (gamma=1.2)", k] <- mean(upper_bound_3 - lower_bound_3)

      # Thin plate spline
      train_df <- data.frame(x1 = x_train[, 1], x2 = x_train[, 2], y = y_obs)
      x_test_df <- data.frame(x1 = x_test[, 1], x2 = x_test[, 2])

      model_tps <- gam(y ~ s(x1, x2, k = 200), data = train_df)
      pred_tps <- predict(model_tps, newdata = x_test_df, se.fit = TRUE)
      y_tps <- pred_tps$fit
      lower_bound_tps <- pred_tps$fit - qnorm(upper_level) * pred_tps$se.fit
      upper_bound_tps <- pred_tps$fit + qnorm(upper_level) * pred_tps$se.fit
      mse_array[row_idx, "Thin Plate Spline", k] <- mean((y_tps - y_test)^2)
      censor_mse_array[row_idx, "Thin Plate Spline", k] <- mean(sort((y_tps - y_test)^2)[censor_lo:censor_hi])
      cov_array[row_idx, "Thin Plate Spline", k] <- mean((y_test >= lower_bound_tps) & (y_test <= upper_bound_tps))
      bandwidth_array[row_idx, "Thin Plate Spline", k] <- mean(upper_bound_tps - lower_bound_tps)

      # Tensor product spline
      model_tep <- gam(y ~ te(x1, x2, k = 15), data = train_df)
      pred_tep <- predict(model_tep, newdata = x_test_df, se.fit = TRUE)
      y_tep <- pred_tep$fit
      lower_bound_tep <- pred_tep$fit - qnorm(upper_level) * pred_tep$se.fit
      upper_bound_tep <- pred_tep$fit + qnorm(upper_level) * pred_tep$se.fit
      mse_array[row_idx, "Tensor Product Spline", k] <- mean((y_tep - y_test)^2)
      censor_mse_array[row_idx, "Tensor Product Spline", k] <- mean(sort((y_tep - y_test)^2)[censor_lo:censor_hi])
      cov_array[row_idx, "Tensor Product Spline", k] <- mean((y_test >= lower_bound_tep) & (y_test <= upper_bound_tep))
      bandwidth_array[row_idx, "Tensor Product Spline", k] <- mean(upper_bound_tep - lower_bound_tep)

      # Gaussian process regression
      gp_model <- gpkm(X = x_train, Z = y_obs, kernel = k_Matern52(D = 2))
      pred_gp <- gp_model$pred(x_test, se = TRUE, mean_dist = TRUE)
      mse_array[row_idx, "Gaussian Process", k] <- mean((pred_gp$mean - y_test)^2)
      censor_mse_array[row_idx, "Gaussian Process", k] <- mean(sort((pred_gp$mean - y_test)^2)[censor_lo:censor_hi])
      lower_bound_gp <- pred_gp$mean - qnorm(upper_level) * pred_gp$se
      upper_bound_gp <- pred_gp$mean + qnorm(upper_level) * pred_gp$se
      cov_array[row_idx, "Gaussian Process", k] <- mean((y_test >= lower_bound_gp) & (y_test <= upper_bound_gp))
      bandwidth_array[row_idx, "Gaussian Process", k] <- mean(upper_bound_gp - lower_bound_gp)

      # BASS
      bass_model <- bass(xx = x_train, y = y_obs, nmcmc = 10000, nburn = 9000, verbose = FALSE)
      pred_bass <- predict(bass_model, x_test)
      y_bass <- colMeans(pred_bass)
      mse_array[row_idx, "BASS", k] <- mean((y_bass - y_test)^2)
      censor_mse_array[row_idx, "BASS", k] <- mean(sort((y_bass - y_test)^2)[censor_lo:censor_hi])
      lower_bound_bass <- apply(pred_bass, 2, quantile, probs = lower_level)
      upper_bound_bass <- apply(pred_bass, 2, quantile, probs = upper_level)
      cov_array[row_idx, "BASS", k] <- mean((y_test >= lower_bound_bass) & (y_test <= upper_bound_bass))
      bandwidth_array[row_idx, "BASS", k] <- mean(upper_bound_bass - lower_bound_bass)
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

