
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



# Load --------------------------------------------------------------------

# Predictor/Residuals
load("/data1/famface01/analysis/misc/320_roi_task_activity/41_predictors.rda", verbose=T)
## combine pred stuff together
colnames(res.multinomial$predprobs) <- colnames(res.multinomial$residprobs)
preds1 <- cbind(res.gaussian$predvals, res.multinomial$predprobs, res.binomial$predprobs)
preds2 <- cbind(res.gaussian$predvals, res.multinomial$predclass, res.binomial$predclass)
preds2 <- model.matrix(~., data=preds2)
## combine resid stuff together
resids1 <- cbind(res.gaussian$residvals, res.multinomial$residprobs, res.binomial$residprobs)
resids2 <- cbind(res.gaussian$residvals, res.multinomial$residclass, res.binomial$residclass)

# Low-Level Features
low.info <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_sym-shape-texture.csv")
all.equal(vnames, as.character(low.info$X))
low.info <- low.info[,-1]

# Data
load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects



# Convolve ----------------------------------------------------------------

source("/data1/famface01/command/misc/face_representations/misc/convolve_funs.R")
suppressMessages(library(bigmemory))

system.time(conv.faces <- llply(dat.vols, function(dat) {
  convolve.features.worker(dat$basics$timing, dat$features$face, verbose=F, parallel=F)
}, .parallel=T))

#
# Predicted
#

system.time(conv.preds1 <- llply(dat.vols, function(dat, mat) {
  # Assign the values in order of the video presentation
  ref.vnames <- as.character(dat$basics$timing$video)
  full.mat <- matrix(0, length(ref.vnames), ncol(mat))
  for (i in 1:nrow(mat)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    full.mat[inds,] <- as.matrix(mat[i,])
  }
  colnames(full.mat) <- colnames(mat)
  full.mat <- scale(full.mat, center=T, scale=F)
  
  # Convolution
  convs.mat <- convolve.features.byvid(dat$basics, full.mat, parallel=F, verbose=F)
  colnames(convs.mat) <- colnames(mat)
  
  convs.mat
}, .parallel=T, mat=preds1))

system.time(conv.preds2 <- llply(dat.vols, function(dat, mat) {
  # Assign the values in order of the video presentation
  ref.vnames <- as.character(dat$basics$timing$video)
  full.mat <- matrix(0, length(ref.vnames), ncol(mat))
  for (i in 1:nrow(mat)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    full.mat[inds,] <- as.matrix(mat[i,])
  }
  colnames(full.mat) <- colnames(mat)
  full.mat <- scale(full.mat, center=T, scale=F)
  
  # Convolution
  convs.mat <- convolve.features.byvid(dat$basics, full.mat, parallel=F, verbose=F)
  colnames(convs.mat) <- colnames(mat)
  
  convs.mat
}, .parallel=T, mat=preds2))

#
# Residuals
#
system.time(conv.resids1 <- llply(dat.vols, function(dat, mat) {
  # Assign the values in order of the video presentation
  ref.vnames <- as.character(dat$basics$timing$video)
  full.mat <- matrix(0, length(ref.vnames), ncol(mat))
  for (i in 1:nrow(mat)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    full.mat[inds,] <- as.matrix(mat[i,])
  }
  colnames(full.mat) <- colnames(mat)
  full.mat <- scale(full.mat, center=T, scale=F)
  
  # Convolution
  convs.mat <- convolve.features.byvid(dat$basics, full.mat, parallel=F, verbose=F)
  colnames(convs.mat) <- colnames(mat)
  
  convs.mat
}, .parallel=T, mat=resids1))

system.time(conv.resids2 <- llply(dat.vols, function(dat, mat) {
  # Assign the values in order of the video presentation
  ref.vnames <- as.character(dat$basics$timing$video)
  full.mat <- matrix(0, length(ref.vnames), ncol(mat))
  for (i in 1:nrow(mat)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    full.mat[inds,] <- as.matrix(mat[i,])
  }
  colnames(full.mat) <- colnames(mat)
  full.mat <- scale(full.mat, center=T, scale=F)
  
  # Convolution
  convs.mat <- convolve.features.byvid(dat$basics, full.mat, parallel=F, verbose=F)
  colnames(convs.mat) <- colnames(mat)
  
  convs.mat
}, .parallel=T, mat=resids2))

# takes 22secs! but is ok cuz have others
system.time(convs.lows <- llply(dat.vols, function(dat) {
  # rearrange/assign
  ref.vnames <- as.character(dat$basics$timing$video)
  low.info2 <- matrix(0, length(ref.vnames), ncol(low.info))
  for (i in 1:nrow(low.info)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    low.info2[inds,] <- as.matrix(low.info[rep(i,2),])
  }
  convs.low <- convolve.features.byvid(dat$basics, low.info2, parallel=F, verbose=F)
  convs.low
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
ssubs    <- rep(subjects, sapply(conv.faces, nrow))
subvars  <- model.matrix(~ssubs+0)
sfaces   <- do.call(rbind, conv.faces)
spreds1  <- do.call(rbind, conv.preds1)
spreds2  <- do.call(rbind, conv.preds2)
sresids1 <- do.call(rbind, conv.resids1)
sresids2 <- do.call(rbind, conv.resids2)
slows    <- do.call(rbind, convs.lows)
smcs     <- do.call(rbind, lst.mcs)
sfds     <- do.call(rbind, lst.fds)
colnames(sfds) <- c("rel", "abs")

# these results remove the effect of each subject
scfaces   <- do.call(rbind, lapply(conv.faces, scale, center=T, scale=F))
scpreds1  <- do.call(rbind, lapply(conv.preds1, scale, center=T, scale=F))
scpreds2  <- do.call(rbind, lapply(conv.preds2, scale, center=T, scale=F))
scresids1 <- do.call(rbind, lapply(conv.resids1, scale, center=T, scale=F))
scresids2 <- do.call(rbind, lapply(conv.resids2, scale, center=T, scale=F))
sclows    <- do.call(rbind, lapply(convs.lows, scale, center=T, scale=F))
scmcs     <- do.call(rbind, lapply(lst.mcs, scale, center=T, scale=F))
scfds     <- do.call(rbind, lapply(lst.fds, scale, center=T, scale=F))
colnames(scfds) <- c("rel", "abs")

# Save
ofile <- "/data1/famface01/analysis/misc/320_roi_task_activity/24_predict_convolve.rda"
save(ssubs, subvars, sfaces, spreds1, spreds2, sresids1, sresids2, slows, smcs, sfds, 
     scfaces, scpreds1, scpreds2, scresids1, scresids2, sclows, scmcs, scfds, 
     file=ofile)
