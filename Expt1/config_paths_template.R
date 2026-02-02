# Root project folder

# this object should be set to your overarching Expt 1 directory
# For example, my Expt 1 direcotry contains R/, data/, globals.R, and thesis_theme.R
path <- "path/to/Expt1/directory"

# this object should be the directory within the Expt 1 directory where all the data (csvs) are stored
data <- "path/to/Expt1/directory/datafolder"

# Folder for plots (absolute path)
# I prevent this folder from appearing on the repo by including it in the .gitignore
plots_folder <- file.path(path, "plots")

# this makes sure that the plots_folder will be created in case it doesn't already exist
if (!dir.exists(plots_folder)) {
  dir.create(plots_folder, showWarnings = FALSE, recursive = TRUE)
}

# sets the working directory
setwd(path)
