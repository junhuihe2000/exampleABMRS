# Set the CRAN mirror globally
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages = c("devtools", "oro.nifti", "neurobase", "mgcv")

# Identify any missing packages
missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]

# Install any that are missing
if(length(missing_packages)) {
  install.packages(missing_packages)
}

devtools::install_github("junhuihe2000/ABMRS")
