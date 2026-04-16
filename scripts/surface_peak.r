## Surface peak detection simulation
## Designed for SLURM array jobs:
## one array task runs one or more experiment replicates

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
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
job_id <- get_job_id(args)
rep_count <- get_rep_count(args)

set.seed(1234 + job_id)

out_dir <- file.path(repo_root, "output", "surface_peak")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

job_tag <- sprintf("job-%04d", job_id)
peak_out <- file.path(out_dir, sprintf("surface_peak_%s.rds", job_tag))

# ---------- functions ----------
f1 <- function(x) {
  y = -40*(exp(-10*(x[,1]-0.3)^2-10*(x[,2]-0.3)^2) - exp(-10*(x[,1]-0.7)^2-10*(x[,2]-0.7)^2))
  return(y)
}

# ---------- simulation setup ----------
p_1 <- c(0.2864, 0.2864)
p_2 <- c(0.7136, 0.7136)

ms <- c(150, 260, 500)
noise <- 0.8
xmin <- c(0, 0)
xmax <- c(1, 1)

distances_array <- array(NA, dim = c(length(ms), 2, rep_count),
                         dimnames = list(
                           m = as.character(ms),
                           peak = c("negative", "positive"),
                           replicate = as.character(seq_len(rep_count))
                         ))

# ---------- simulation ----------
time_start <- Sys.time()
message("Starting job ", job_id, " with ", rep_count, " replicate(s)")

for (k in seq_len(rep_count)) {
  message("Replicate ", k)
  
  for (s in seq_along(ms)) {
    m <- ms[s]
    message("  Training size: ", m)
    
    # generate training samples
    x <- cbind(runif(m), runif(m))
    y <- f1(x)
    y_obs <- y + rnorm(m, 0, noise)
    
    # run surface ABMRS
    my_mebars <- mebars(x, y_obs, xmin = xmin, xmax = xmax, ks = c(5,5), times = c(2,2), degrees = c(1,1))
    my_mebars$rjmcmc(burns = 20000, steps = 20000)
    
    knots <- my_mebars$knots()
    l <- length(knots)
    
    peak_samples <- list()
    
    for (i in 1:l) {
      xi_1 <- knots[[i]][[1]]
      xi_2 <- knots[[i]][[2]]
      # fit linear spline regression
      z <- tensor_spline(x, xis = list(xi_1, xi_2), degrees = c(1,1), intercepts = c(TRUE,TRUE), xmin = xmin, xmax = xmax)
      reg_df <- data.frame(z, y = y_obs)
      lm_model <- lm(y ~ ., data = reg_df)
      
      xis <- as.matrix(expand.grid(xi_1, xi_2))
      is_peak <- array(NA, c(nrow(xis), 2))
      y_pred <- rep(NA, nrow(xis))
      for (j in 1:nrow(xis)) {
        xi <- xis[j, ]
        x_neigh <- as.matrix(expand.grid(seq(xi[1]-0.02, xi[1]+0.02, 0.01), 
                              seq(xi[2]-0.02, xi[2]+0.02, 0.01)))
        z_hat <- tensor_spline(rbind(xi, x_neigh), xis = list(xi_1, xi_2), degrees = c(1,1), intercepts = c(TRUE,TRUE), xmin = xmin, xmax = xmax)
        y_hat <- predict(lm_model, data.frame(z_hat))
        
        is_peak[j, 1] <- y_hat[1] == min(y_hat)
        is_peak[j, 2] <- y_hat[1] == max(y_hat)
        y_pred[j] <- y_hat[1]
      }
      
      peak_df <- data.frame(xis, is_peak, y_pred)
      colnames(peak_df) <- c("x_1", "x_2", "min", "max", "y_hat")
      
      peak_samples[[i]] <- peak_df
    }
    
    neg_peak_points <- list()
    pos_peak_points <- list()
    peak_count <- array(NA, c(l, 2))
    
    for (i in 1:l) {
      peak_df <- peak_samples[[i]]
      # negative peak
      neg_peak <- as.matrix(peak_df[peak_df[,3]==1, c(1:2,5)])
      neg_peak_points[[i]] <- neg_peak
      peak_count[i, 1] <- nrow(neg_peak)
      
      # positive peak
      pos_peak <- as.matrix(peak_df[peak_df[,4]==1, c(1:2,5)])
      pos_peak_points[[i]] <- pos_peak
      peak_count[i, 2] <- nrow(pos_peak)
    }
    
    negative_peaks <- c()
    positive_peaks <- c()
    for (i in 1:l) {
      if (peak_count[i,1]>0) {
        neg_peak <- neg_peak_points[[i]]
        negative_peaks <- rbind(negative_peaks, neg_peak[which.min(neg_peak[,3]),])
      }
      if (peak_count[i,2]>0) {
        pos_peak <- pos_peak_points[[i]]
        positive_peaks <- rbind(positive_peaks, pos_peak[which.max(pos_peak[,3]),])
      }
    }
    
    phat_1 <- colMeans(negative_peaks[,1:2])
    phat_2 <- colMeans(positive_peaks[,1:2])
    d_1 <- sqrt(sum((phat_1 - p_1)^2))
    d_2 <- sqrt(sum((phat_2 - p_2)^2))
    
    distances_array[s, 1, k] <- d_1
    distances_array[s, 2, k] <- d_2
    
    message("    Distance to negative peak: ", round(d_1, 4))
    message("    Distance to positive peak: ", round(d_2, 4))
  }
}

time_end <- Sys.time()
message("Elapsed: ", format(time_end - time_start))

saveRDS(distances_array, file = peak_out)
message("Wrote results to ", peak_out)