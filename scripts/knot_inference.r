## Knot inference simulation
## Designed for SLURM array jobs:
## one array task runs one or more experiment replicates

suppressPackageStartupMessages({
  library(ABMRS)
  library(splines)
  library(segmented)
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

# Source LSE function
lse_path <- file.path(repo_root, "scripts", "Linear_Spline_Estimation.R")
if (file.exists(lse_path)) {
  source(lse_path)
} else {
  warning("Linear_Spline_Estimation.R not found at ", lse_path)
}

out_dir <- file.path(repo_root, "output", "knot_inference")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

job_tag <- sprintf("job-%04d", job_id)
raw_error_out <- file.path(out_dir, sprintf("knot_error_raw_%s.rds", job_tag))
raw_quantiles_out <- file.path(out_dir, sprintf("knot_quantiles_raw_%s.rds", job_tag))

# ---------- functions ----------
f1 <- function(x) {
  knots <- c(0.5)
  beta <- c(1, -1, 1)
  B <- bs(x, knots = knots, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1))
  y <- B %*% beta
  return(y)
}

f2 <- function(x) {
  knots <- c(0.3, 0.7)
  beta <- c(2, -1, -2, -1)
  B <- bs(x, knots = knots, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1))
  y <- B %*% beta
  return(y)
}

f3 <- function(x) {
  knots <- c(0.2, 0.2, 0.5, 0.7)
  beta <- c(0, -1, 1, 0, 1, 0)
  B <- bs(x, knots = knots, degree = 1, intercept = TRUE, Boundary.knots = c(0, 1))
  y <- B %*% beta
  return(y)
}

# ---------- simulation setup ----------
noises <- c(0.4, 0.3, 0.4)
xis <- list(c(0.5), c(0.3, 0.7), c(0.2, 0.2, 0.5, 0.7))
fs <- list(f1, f2, f3)
ks <- c(1, 2, 4)
idxs <- list(c(1), c(2:3), c(4:7))

m_trains <- c(200, 500)
scenarios <- c("Case 1", "Case 2", "Case 3")
cases <- c("Case 1 Knot 1", "Case 2 Knot 1", "Case 2 Knot 2", "Case 3 Knot 1", "Case 3 Knot 2", "Case 3 Knot 3", "Case 3 Knot 4")
methods <- c("ABMRS", "LSE", "Segmented")
probs <- c(0.025, 0.125, 0.25, 0.75, 0.875, 0.975)
probs_names <- paste0(probs * 100, "%")

row_labels <- paste(rep(cases, times = length(m_trains)), rep(m_trains, each = length(cases)), sep = "-")

error_array <- array(
  NA_real_,
  dim = c(length(row_labels), length(methods), rep_count),
  dimnames = list(row_labels, methods, NULL)
)

knot_quantiles_array <- array(
  NA_real_,
  dim = c(length(row_labels), length(probs_names), rep_count),
  dimnames = list(row_labels, probs_names, NULL)
)

time_start <- Sys.time()
message("Starting job ", job_id, " with ", rep_count, " replicate(s)")

for (k in seq_len(rep_count)) {
  message("Replicate ", k)

  for (s in seq_along(m_trains)) {
    m_train <- m_trains[s]
    row_offset <- (s - 1L) * length(cases)
    message("Train size ", m_train)

    for (i in seq_along(fs)) {
      message("  Running ", scenarios[i])
      
      # generate training data
      k_knots <- ks[i]
      f <- fs[[i]]
      noise <- noises[[i]]
      xi <- xis[[i]]
      idx <- idxs[[i]]
      row_idx <- row_offset + idx
      
      x_train <- runif(m_train)
      y_train <- f(x_train)
      y_obs <- c(y_train + rnorm(m_train, 0, noise))

      # ABMRS
      ebars_linear <- ebars(x_train, y_obs, k = k_knots, times = 2, degree = 1, intercept = TRUE)
      ebars_linear$rjmcmc(burns = 20000, steps = 20000)
      knots_ebars <- rbind(sapply(ebars_linear$knots(), function(xi) {
        return(xi)
      }))
      xi_ebars <- rowMeans(knots_ebars)

      error_array[row_idx, "ABMRS", k] <- abs(xi - xi_ebars)
      knot_quantiles_array[row_idx, , k] <- t(apply(knots_ebars, 1, function(xi) {
        quantile(xi, probs = probs)
      }))

      # Linear Spline Estimation
      if (k_knots == 4) {
        init_lse <- quantile(x_train, probs = seq(0, 1, length.out = 5)[2:4])
        xi_lse <- sort(LSE(x_train, y_obs, init_lse))
        xi_lse <- c(xi_lse[1], xi_lse[1], xi_lse[2], xi_lse[3])
      } else {
        init_lse <- quantile(x_train, probs = seq(0, 1, length.out = (k_knots + 2))[2:(k_knots + 1)])
        xi_lse <- sort(LSE(x_train, y_obs, init_lse))
      }
      xi_lse <- sapply(xi_lse, function(x) {
        max(min(x, 1), 0)
      })
      error_array[row_idx, "LSE", k] <- abs(xi - xi_lse)

      # segmented method
      init_seg <- quantile(x_train, probs = seq(0, 1, length.out = (k_knots + 2))[2:(k_knots + 1)])
      fit.lm <- lm(y_obs ~ x_train)
      fit.seg <- tryCatch(
        {
          segmented(fit.lm, seg.Z = ~x_train, npsi = k_knots)
        },
        error = function(e) {
          warning("Segmented failed for ", scenarios[i], " (m=", m_train, "): ", e$message)
          return(NULL)
        }
      )
      
      if (!is.null(fit.seg) && !is.null(fit.seg$psi) && nrow(fit.seg$psi) == k_knots) {
        xi_seg <- sort(fit.seg$psi[, 2])
        error_array[row_idx, "Segmented", k] <- abs(xi - xi_seg)
      } else {
        # Segmented failed or returned wrong number of breakpoints
        error_array[row_idx, "Segmented", k] <- NA_real_
      }
    }
  }
}

time_end <- Sys.time()
message("Elapsed: ", format(time_end - time_start))

saveRDS(error_array, file = raw_error_out)
message("Wrote raw error array to ", raw_error_out)

saveRDS(knot_quantiles_array, file = raw_quantiles_out)
message("Wrote raw quantiles array to ", raw_quantiles_out)