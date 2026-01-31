# Share number of knots across dimensions simulation
# Designed for SLURM array jobs:
# one array task runs one or more experiment replicates

suppressPackageStartupMessages({
  library(ABMRS)
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

# output paths
out_dir <- file.path(repo_root, "output", "share_knots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
job_tag <- sprintf("job-%04d", job_id)
raw_out <- file.path(out_dir, sprintf("share_knots_raw_%s.rds", job_tag))

# ---------- simulation setup ----------

# Define the true function
# f(x) = sin(2πx1) + cos(8πx2)
f <- function(x) {
  sin(2 * pi * x[, 1]) + cos(8 * pi * x[, 2])
}

noise_sd <- 0.1

# Define the domain
xmin <- c(0, 0)
xmax <- c(1, 1)

# Simulation parameters
ms <- c(200, 500)
ks <- list(c(2, 8), c(2, 2), c(5, 5), c(8, 8))

# Fixed test set
set.seed(1234)
m_test <- 500
X_test <- matrix(runif(m_test * 2), ncol = 2)
y_test <- f(X_test)

# Results array
row_labels <- paste0("m=", ms)
col_labels <- paste0("ks=", sapply(ks, function(k) paste(k, collapse = ",")))
mse_array <- array(NA_real_, 
                   dim = c(length(ms), length(ks), rep_count), 
                   dimnames = list(row_labels, col_labels, NULL))

# ---------- run experiments ----------
message("Starting job ", job_id, " with ", rep_count, " replicate(s)")

for (r in seq_len(rep_count)) {
  message("Replicate ", r)
  
  for (i in seq_along(ms)) {
    m_train <- ms[i]
    
    # Generate training data
    X_train <- matrix(runif(m_train * 2), ncol = 2)
    y_train <- f(X_train)
    y_obs <- y_train + rnorm(m_train, sd = noise_sd)
    
    for (j in seq_along(ks)) {
      k <- ks[[j]]
      message(sprintf("  m=%d, k=(%d,%d)", m_train, k[1], k[2]))
      
      # Fit model
      mymebars <- mebars(
        x = X_train,
        y = y_obs,
        xmin = xmin,
        xmax = xmax,
        times = c(2, 2),
        ks = k
      )
      mymebars$rjmcmc(burns = 10000, steps = 10000)
      
      # Predict and compute MSE
      pred <- mymebars$predict(X_test)
      y_pred <- rowMeans(pred)
      mse_array[i, j, r] <- mean((y_pred - y_test)^2)
    }
  }
}

# ---------- summarize & write ----------

saveRDS(mse_array, file = raw_out)
message("Wrote raw MSE array to ", raw_out)
