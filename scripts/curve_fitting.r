# Univariate spline simulation
# Designed for SLURM array jobs:
# one array task runs one or more experiment replicates

suppressPackageStartupMessages({
  library(ABMRS)
  library(splines)
  library(GauPro)
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
repo_root  <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)
job_id <- get_job_id(args)
rep_count <- get_rep_count(args)

set.seed(1234 + job_id)

bars_path <- file.path(repo_root, "scripts", "bars.r")
if (file.exists(bars_path)) {
  source(bars_path)
} else {
  warning("bars.r not found at ", bars_path, "; BARS comparison will be skipped")
}

# output paths
out_dir <- file.path(repo_root, "output", "curve_fitting")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
job_tag <- sprintf("job-%04d", job_id)
raw_out <- file.path(out_dir, sprintf("curve_mse_raw_%s.rds", job_tag))
raw_censor_out <- file.path(out_dir, sprintf("curve_mse_censor_raw_%s.rds", job_tag))

# ---------- simulation setup ----------
p <- 0.025
m_trains <- c(200, 500)
m_tests  <- c(200, 500)
noises <- c(2, 4, 1)

f1 <- function(x) 10 * sin(2 * pi * x)
f2 <- function(x) {
  knot <- c(0.4, 0.4, 0.4, 0.4, 0.7)
  beta <- matrix(c(2, -5, 5, 2, -3, -1, 2), ncol = 1)
  B <- ns(x, knots = knot, intercept = TRUE, Boundary.knots = c(0, 1))
  c(5 * B %*% beta)
}
f3 <- function(x) 40 * (x %% 0.167)
fs <- list(f1, f2, f3)

cases <- c("Case 1", "Case 2", "Case 3")
methods <- c("ABMRS-1", "ABMRS-2", "ABMRS-3", "BARS", "Smooth Spline", "Gaussian Process")
row_labels <- paste(rep(cases, times = length(m_trains)), rep(m_trains, each = length(cases)), sep = "-")
mse_array <- array(NA_real_, dim = c(length(row_labels), length(methods), rep_count), dimnames = list(row_labels, methods, NULL))
censor_mse_array <- array(NA_real_, dim = c(length(row_labels), length(methods), rep_count), dimnames = list(row_labels, methods, NULL))

# ---------- run experiments ----------
message("Starting job ", job_id, " with ", rep_count, " replicate(s)")
for (k in seq_len(rep_count)) {
  message("Replicate ", k)

  for (s in seq_along(m_trains)) {
    m_train <- m_trains[s]
    m_test  <- m_tests[s]
    row_offset <- (s - 1L) * length(cases)

    x_train <- seq(0, 1, length.out = m_train)
    x_new <- runif(m_test, 0, 1)

    for (i in seq_along(fs)) {
      row_idx <- row_offset + i
      y_train <- fs[[i]](x_train)
      y_obs <- y_train + rnorm(m_train, 0, noises[[i]])
      y_new <- fs[[i]](x_new)

      # ABMRS gamma = 0.5
      ebars_1 <- ebars(x_train, y_obs, gamma = 0.5, times = 1, intercept = TRUE)
      ebars_1$rjmcmc(burns = 5000, steps = 5000)
      pred_1 <- ebars_1$predict(x_new)
      y_ebars_1 <- rowMeans(pred_1)
      mse_array[row_idx, "ABMRS-1", k] <- mean((y_new - y_ebars_1)^2)
      censor_mse_array[row_idx, "ABMRS-1", k] <- mean(sort((y_new - y_ebars_1)^2)[round(p * m_test):round((1 - p) * m_test)])

      # ABMRS gamma = 1
      ebars_2 <- ebars(x_train, y_obs, gamma = 1, times = 3, intercept = TRUE)
      ebars_2$rjmcmc(burns = 5000, steps = 5000)
      pred_2 <- ebars_2$predict(x_new)
      y_ebars_2 <- rowMeans(pred_2)
      mse_array[row_idx, "ABMRS-2", k] <- mean((y_new - y_ebars_2)^2)
      censor_mse_array[row_idx, "ABMRS-2", k] <- mean(sort((y_new - y_ebars_2)^2)[round(p * m_test):round((1 - p) * m_test)])

      # ABMRS gamma = 1.5
      ebars_3 <- ebars(x_train, y_obs, gamma = 1.5, times = 3, intercept = TRUE)
      ebars_3$rjmcmc(burns = 5000, steps = 5000)
      pred_3 <- ebars_3$predict(x_new)
      y_ebars_3 <- rowMeans(pred_3)
      mse_array[row_idx, "ABMRS-3", k] <- mean((y_new - y_ebars_3)^2)
      censor_mse_array[row_idx, "ABMRS-3", k] <- mean(sort((y_new - y_ebars_3)^2)[round(p * m_test):round((1 - p) * m_test)])

      # BARS
      if (exists("BARS") && exists("bars_fit")) {
        bars_conf <- BARS(x_train, y_obs, burns_in = 1000, effective_length = 1000)
        y_bars <- bars_fit(x_new, bars_conf)
        mse_array[row_idx, "BARS", k] <- mean((y_new - y_bars)^2)
        censor_mse_array[row_idx, "BARS", k] <- mean(sort((y_new - y_bars)^2)[round(p * m_test):round((1 - p) * m_test)])
      }

      # Smooth spline
      splinesmooth <- smooth.spline(x_train, y_obs)
      y_ss <- predict(splinesmooth, x_new)$y
      mse_array[row_idx, "Smooth Spline", k] <- mean((y_new - y_ss)^2)
      censor_mse_array[row_idx, "Smooth Spline", k] <- mean(sort((y_new - y_ss)^2)[round(p * m_test):round((1 - p) * m_test)])

      # Gaussian process regression
      gp_model <- gpkm(X = cbind(x_train), Z = y_obs, kernel = k_Gaussian(0.1))
      y_gp <- gp_model$pred(cbind(x_new))
      mse_array[row_idx, "Gaussian Process", k] <- mean((y_new - y_gp)^2)
      censor_mse_array[row_idx, "Gaussian Process", k] <- mean(sort((y_new - y_gp)^2)[round(p * m_test):round((1 - p) * m_test)])
    }
  }
}

# ---------- summarize & write ----------

saveRDS(mse_array, file = raw_out)
saveRDS(censor_mse_array, file = raw_censor_out)
message("Wrote raw MSE array to ", raw_out)
message("Wrote raw censored MSE array to ", raw_censor_out)