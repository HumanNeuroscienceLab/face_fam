# This convolves the trait and demographic features with the HRF


# Setup -------------------------------------------------------------------

# Path for my own packages
if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

# General workhorse function
library(plyr)

# Parallelization
suppressMessages(library(doMC))
registerDoMC(24)

# If we were to use ggplot2 and other plotting needs
library(RColorBrewer)
library(ggplot2)
source("/data1/famface01/command/encoding/FirstPass/plot_functions.R")

# This gets us the simple_lm function for faster GLMs
source("/data1/famface01/command/encoding/SubRois_Unfam/lm_functions.R")

## For heatmap.2
#suppressMessages(library(gplots))
## For our own custom convolution
#source("/data1/famface01/command/encoding/ShapesAnalysis/000_funs_load_convolve.R")

subjects <- sprintf("sub%02i", 1:6) # again skipping sub01
#indir    <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
rnames   <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", 'r.Amyg', 
              "l.vATL", "l.FFA", "l.OFA", "l.EBA", 'l.Amyg')



# Load -------------------------------------------------------------------

# ROI Data
load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects

# Features (don't do anything here with the low level features)
df.demos <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_df-demos.csv")[,-c(1:2)]
vnames <- read.table("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_demo-vnames.txt")
vnames <- as.character(vnames[,1])
all.mat <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-mat.csv")[,-1]

# Center all.mat columns
all.grps <- as.character(read.table("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-groups.txt")[,1])
all.cmat <- all.mat
all.cmat[,all.grps %in% c("traits","age","maketup")] <- scale(all.cmat[,all.grps %in% c("traits","age","maketup")], center=T, scale=F)


# Convolve ----------------------------------------------------------------

source("/data1/famface01/command/encoding/ShapesAnalysis/000_funs_load_convolve.R")
suppressMessages(library(bigmemory))

#ofile <- "/data1/famface01/analysis/misc/320_roi_task_activity/22_roi_predicts_vs_shapes_convs.rda"

system.time(conv.faces <- llply(dat.vols, function(dat) {
  convolve.features.worker(dat$basics$timing, dat$features$face, verbose=F, parallel=F)
}, .parallel=T))

system.time(conv.mats <- llply(dat.vols, function(dat) {
  # rearrange/assign
  ref.vnames <- as.character(dat$basics$timing$video)
  all.cmat2 <- matrix(0, length(ref.vnames), ncol(all.cmat))
  for (i in 1:nrow(all.cmat)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    for (ind in inds) all.cmat2[ind,] <- as.numeric(all.cmat[i,])
  }
  
  convs.comps <- convolve.features.byvid(dat$basics, all.cmat2, parallel=F, verbose=F)
  colnames(convs.comps) <- colnames(all.cmat)
  
  convs.comps
}, .parallel=T))

# get motion
load.mc <- function(subj) {
  funcdir  <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  df.paths <- read.table(file.path(funcdir, "df_paths.txt"), header=T)
  inds     <- df.paths$inindex[df.paths$name == "unfam_vids"]
  fpaths   <- sprintf("%s/mc/func_run%02i_dfile.1D", funcdir, inds)
  
  mcs <- ldply(fpaths, function(fpath) {
    x <- read.table(fpath)
    x <- as.matrix(x)
    x <- scale(x, scale=F, center=T)
    x
  })
  mcs <- as.matrix(mcs)
  
  colnames(mcs) <- c("roll", "pitch", "yaw", "dS", "dL", "dP")
  
  mcs
}
lst.mcs <- lapply(subjects, load.mc)
names(lst.mcs) <- subjects

# FDs
load.fds <- function(subj) {
  funcdir  <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  df.paths <- read.table(file.path(funcdir, "df_paths.txt"), header=T)
  inds     <- df.paths$inindex[df.paths$name == "unfam_vids"]
  
  fpaths   <- sprintf("%s/mc/func_run%02i_fd_rel.1D", funcdir, inds)
  fds_rel  <- ldply(fpaths, function(fpath) {
    x <- read.table(fpath)
    x <- as.matrix(x)
    x <- scale(x, scale=F, center=T)
    x
  })
  fds_rel <- as.matrix(fds_rel)
  
  fpaths   <- sprintf("%s/mc/func_run%02i_fd_abs.1D", funcdir, inds)
  fds_abs  <- ldply(fpaths, function(fpath) {
    x <- read.table(fpath)
    x <- as.matrix(x)
    x <- scale(x, scale=F, center=T)
    x
  })
  fds_abs <- as.matrix(fds_abs)
  
  cbind(rel=fds_rel, abs=fds_abs)
}
lst.fds <- lapply(subjects, load.fds)
names(lst.fds) <- subjects

# combine things
ssubs   <- rep(subjects, sapply(conv.faces, nrow))
subvars <- model.matrix(~ssubs+0)
sfaces  <- do.call(rbind, conv.faces)
smats   <- do.call(rbind, conv.mats)
smcs    <- do.call(rbind, lst.mcs)
sfds    <- do.call(rbind, lst.fds)
colnames(sfds) <- c("rel", "abs")

# these results remove the effect of each subject
subvars2 <- subvars[,-1]
sfaces2  <- do.call(rbind, lapply(conv.faces, scale, center=T, scale=F))
smats2   <- do.call(rbind, lapply(conv.mats, scale, center=T, scale=F))
smcs2    <- do.call(rbind, lapply(lst.mcs, scale, center=T, scale=F))
sfds2    <- do.call(rbind, lapply(lst.fds, scale, center=T, scale=F))
colnames(sfds2) <- c("rel", "abs")



# Save --------------------------------------------------------------------

ofile <- "/data1/famface01/analysis/misc/320_roi_task_activity/22_convolve_facefeats.rda"
save(ssubs, subvars, sfaces, smats, smcs, sfds, 
     subvars2, sfaces2, smats2, smcs2, sfds2, 
     all.grps, 
     file=ofile)


