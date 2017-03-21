
# This script will generate the beta-series for the average signal in a brain
# region


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

subjects <- sprintf("sub%02i", 1:6) # again skipping sub01
#indir    <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
rnames   <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", 'r.Amyg', 
              "l.vATL", "l.FFA", "l.OFA", "l.EBA", 'l.Amyg')



# Load --------------------------------------------------------------------

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

# Std Data
load("/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda", verbose=T)



# Functions ---------------------------------------------------------------

source("/data1/famface01/command/encoding/ShapesAnalysis/000_funs_load_convolve.R")

#' mags => magnitudes, can be a value, vector, or matrix. if matrix each column is some vector
#' frate indicates how much to upsample
gen.events <- function(onsets, durs, tot.time, mags=1, frate=1/24)
{
  # For comparison:
  # neuRosim::stimfunction(tot.time, onsets, durs, frate)
  
  if (length(mags) == 1) mags <- rep(mags, length(onsets))
  if (length(durs) == 1) durs <- rep(durs, length(onsets))
  if (length(onsets) != length(durs)) stop("onsets and durs must be same length")
  mags <- as.matrix(mags)
  nc <- ncol(mags)
  
  # Upsample the features (faster way to do this?)
  ## setup upsampled data
  up.feats <- matrix(0, round(tot.time/frate), nc)
  ## find new indices to put things
  onsets  <- round(onsets/frate + 1)
  durs    <- round(durs/frate)
  offsets <- sapply(1:length(onsets), function(i) max(onsets[i] + durs[i] - 1, onsets[i]))
  offsets[offsets > nrow(up.feats)] <- nrow(up.feats)
  ns      <- offsets - onsets + 1
  up.inds <- unlist(lapply(1:length(onsets), function(i) onsets[i]:offsets[i]))
  ## corresponding indices in regularly sampled data
  reg.inds<- rep(1:length(onsets), ns)
  ## upsample
  up.feats[up.inds,] <- mags[reg.inds,]
  
  return(up.feats)
}

library(neuRosim)
apply.convolve <- function(up.feats, runs, hrf_fun=canonicalHRF, frate=1/24, parallel=F, ...) {
  # Get the HRF to convolve
  uruns    <- sort(unique(runs))
  nruns    <- max(uruns)
  ntpts    <- sum(runs==uruns[1])
  tpts     <- seq(frate, ntpts, frate) # 318 * 16runs = 2288 * 24frames
  chrf     <- hrf_fun(tpts, ...)
  up.runs <- rep(1:nruns, each=length(tpts))
  
  # Convolve
  conv.up.feats <- llply(1:nruns, function(irun) {
    ys <- up.feats[up.runs==irun,,drop=F]
    ys <- ys[nrow(ys):1,,drop=F] # same as rev to each column
    s.convs <- convolve.mfftw(chrf, ys)
    for (i in 1:ncol(s.convs)) {
      s.convs[,i] <- s.convs[,i]/max(s.convs[,i]) # normalize
      s.convs[is.na(s.convs[,i]),i] <- 0
    }
    s.convs
  }, .parallel=parallel) # TODO: test if parallel is actually faster!
  conv.up.feats <- do.call(rbind, conv.up.feats)
  conv.up.feats <- as.matrix(conv.up.feats)
  
  return(conv.up.feats)
}

downsample <- function(up.dat, tr=1, frate=1/24) {
  ix <- seq(1, nrow(up.dat), tr/frate)
  down.dat <- up.dat[ix,,drop=F]
  return(down.dat)
}

convolve.hrf <- function(onsets, durs, tot.time, runs, 
                         mags=1, frate=1/24, tr=1, parallel=F)
{
  s1 <- gen.events(onsets, durs, tot.time, mags, frate)
  c1 <- apply.convolve(s1, runs, frate=frate, parallel=parallel, verbose=F)
  d1 <- downsample(c1, tr, frate)
  return(d1)
}

convolve.spm.hrf <- function(onsets, durs, tot.time, runs, nparams=3, 
                             mags=1, frate=1/24, tr=1, parallel=F)
{
  s1 <- gen.events(onsets, durs, tot.time, mags, frate)
  d1s   <- sapply(c(spmt, dspmt, ddspmt)[1:nparams], function(hfun) {
    c1 <- apply.convolve(s1, runs, hrf_fun=hfun, frate=frate, parallel=parallel)
    d1 <- downsample(c1, tr, frate)
    d1
  })
  return(d1s)
}

convolve.afni.hrf <- function(onsets, durs, tot.time, runs, afni.ret, 
                              mags=1, frate=1/24, tr=1, parallel=F)
{
  afni_fun <- function(tpts, ci=1) {
    afni.tpts <- vector("numeric", length(tpts))
    afni.tpts[1:nrow(afni.ret)] <- afni.ret[,ci]
    afni.tpts
  }
  s1 <- gen.events(onsets, durs, tot.time, mags, frate)
  d1s   <- sapply(1:ncol(afni.ret), function(ci) {
    c1 <- apply.convolve(s1, runs, hrf_fun=afni_fun, frate=frate, parallel=parallel, ci=ci)
    d1 <- downsample(c1, tr, frate)
    d1
  })
  return(d1s)
}

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



# Beta-Series -------------------------------------------------------------

bvnames <- sort(levels(dat.vols$sub01$basics$timing$video))

# Do the regular one
sub.bs1 <- laply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  #rdats <- sapply(dat$fmri$dat, rowMeans)
  rdats <- std.rdats2[ssubs==subj,]
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  ## Covars (same throughout)
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.hrf(onsets3, durs3, tot.time, runs, frate=1/12)
  # Day of scan
  days <- factor(runs)
  days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  # Repeat days (on the second and fourth day we repeat the stimuli)
  repeat.day <- factor(vid.timing$run)
  repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  onsets4    <- vid.timing$onset[repeat.day==2]
  durs4      <- vid.timing$duration[repeat.day==2]
  repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs, frate=1/12)
  # Motion
  mc <- load.mc(subj)
  fds <- load.fds(subj)
  
  # Loop through each video and process
  system.time(vid.bs <- laply(vnames, function(vname) {
    # Particular trial
    onsets1   <- vid.timing$onset[vid.timing$video == vname]
    durs1     <- vid.timing$duration[vid.timing$video == vname]
    cur.trial <- convolve.hrf(onsets1, durs1, tot.time, runs, frate=1/12)
    
    # All other trials
    onsets2 <- vid.timing$onset[vid.timing$video != vname]
    durs2   <- vid.timing$duration[vid.timing$video != vname]
    other.trials <- convolve.hrf(onsets2, durs2, tot.time, runs, frate=1/12)
    
    # Run regression
    dmat <- model.matrix(~cur.trial + other.trials + quests + days + repeat.faces + mc + fds)
    ret.fit <- simple_lm(rdats, dmat)
    
    # Get the betas and tvals for the current trial estimate
    bs <- ret.fit$b[2,]
    ts <- ret.fit$tvals[2,]
    
    cbind(betas=bs, tvals=ts)
  }, .parallel=T))
  dimnames(vid.bs)[[1]] <- vnames
  names(dimnames(vid.bs)) <- c("vid", "roi", "measure")
  
  return(vid.bs)
}, .progress="text")

# group ave
grp.bs1 <- apply(sub.bs1, 2:4, mean)


# Save -------------------------------------------------------------

save(grp.bs1, sub.bs1, bvnames, file="/data1/famface01/analysis/misc/320_roi_task_activity/dgamma_betas.rda")


save(sub.bs1, file="/data1/famface01/analysis/misc/320_roi_task_activity/dgamma_betas.rda")

