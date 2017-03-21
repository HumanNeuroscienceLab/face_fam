
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



# Convolve ----------------------------------------------------------------

source("/data1/famface01/command/encoding/ShapesAnalysis/000_funs_load_convolve.R")
suppressMessages(library(bigmemory))

ofile <- "/data1/famface01/analysis/misc/320_roi_task_activity/22_roi_predicts_vs_shapes_convs.rda"

system.time(conv.faces <- llply(dat.vols, function(dat) {
  convolve.features.worker(dat$basics$timing, dat$features$face, verbose=F, parallel=F)
}, .parallel=T))

system.time(conv.fas <- llply(dat.vols, function(dat) {
  # rearrange/assign
  ref.vnames <- as.character(dat$basics$timing$video)
  fac.scores2 <- matrix(0, length(ref.vnames), ncol(fac.scores))
  for (i in 1:nrow(fac.scores)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    fac.scores2[inds,] <- fac.scores[rep(i,2),]
  }
  
  convs.comps <- convolve.features.byvid(dat$basics, fac.scores2, parallel=F, verbose=F)
  convs.comps
}, .parallel=T))

system.time(conv.preds <- llply(dat.vols, function(dat) {
  # rearrange/assign
  ref.vnames <- as.character(dat$basics$timing$video)
  fac.preds2 <- matrix(0, length(ref.vnames), ncol(fac.preds))
  for (i in 1:nrow(fac.scores)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    fac.preds2[inds,] <- fac.preds[rep(i,2),]
  }
  
  convs.comps <- convolve.features.byvid(dat$basics, fac.preds2, parallel=F, verbose=F)
  convs.comps
}, .parallel=T))

system.time(conv.resids <- llply(dat.vols, function(dat) {
  # rearrange/assign
  ref.vnames <- as.character(dat$basics$timing$video)
  fac.resids2 <- matrix(0, length(ref.vnames), ncol(fac.resids))
  for (i in 1:nrow(fac.scores)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    fac.resids2[inds,] <- fac.resids[rep(i,2),]
  }
  
  convs.comps <- convolve.features.byvid(dat$basics, fac.resids2, parallel=F, verbose=F)
  convs.comps
}, .parallel=T))

# takes 22secs! but is ok cuz have others
system.time(convs.lows <- llply(dat.vols, function(dat) {
  low.info <- pca.face.feats
  # rearrange/assign
  ref.vnames <- as.character(dat$basics$timing$video)
  low.info2 <- matrix(0, length(ref.vnames), ncol(low.info))
  for (i in 1:nrow(low.info)) {
    inds <- which(ref.vnames == vnames[i])
    if (length(inds) != 2) stop("error")
    low.info2[inds,] <- low.info[rep(i,2),]
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
ssubs   <- rep(subjects, sapply(conv.faces, nrow))
subvars <- model.matrix(~ssubs+0)
sfaces  <- do.call(rbind, conv.faces)
sfas    <- do.call(rbind, conv.fas)
spreds  <- do.call(rbind, conv.preds)
sresids <- do.call(rbind, conv.resids)
slows   <- do.call(rbind, convs.lows)
smcs    <- do.call(rbind, lst.mcs)
sfds    <- do.call(rbind, lst.fds)
colnames(sfds) <- c("rel", "abs")

# these results remove the effect of each subject
subvars2 <- subvars[,-1]
sfaces2  <- do.call(rbind, lapply(conv.faces, scale, center=T, scale=F))
sfas2    <- do.call(rbind, lapply(conv.fas, scale, center=T, scale=F))
spreds2  <- do.call(rbind, lapply(conv.preds, scale, center=T, scale=F))
sresids2 <- do.call(rbind, lapply(conv.resids, scale, center=T, scale=F))
slows2   <- do.call(rbind, lapply(convs.lows, scale, center=T, scale=F))
smcs2    <- do.call(rbind, lapply(lst.mcs, scale, center=T, scale=F))
sfds2    <- do.call(rbind, lapply(lst.fds, scale, center=T, scale=F))
colnames(sfds2) <- c("rel", "abs")

save(ssubs, subvars, sfaces, sfas, sresids, slows, smcs, sfds, spreds, 
     subvars2, sfaces2, sfas2, sresids2, slows2, smcs2, sfds2, spreds2, 
     file=ofile)


