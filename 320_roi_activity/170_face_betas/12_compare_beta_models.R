
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

# Features
indir <- "/data1/famface01/analysis/misc/320_roi_task_activity"
load(file.path(indir, "20_predict_face_feats.rda"), verbose=T)



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

# So now finally getting our hands dirty

subj <- "sub01"
dat <- dat.vols[[subj]]

rdats <- sapply(dat$fmri$dat, rowMeans)

runs <- dat$fmri$runs
tot.time <- length(runs)

vid.timing <- dat$basics$timing
vnames <- levels(vid.timing$video)

## Covars (same throughout)
# Questions
onsets3 <- vid.timing$onset[vid.timing$question!="none"]
durs3   <- 4
quests  <- convolve.hrf(onsets3, durs3, tot.time, runs)
# Day of scan
days <- factor(runs)
days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
# Repeat days (on the second and fourth day we repeat the stimuli)
repeat.day <- factor(vid.timing$run)
repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
onsets4    <- vid.timing$onset[repeat.day==2]
durs4      <- vid.timing$duration[repeat.day==2]
repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
# Motion
mc <- load.mc(subj)
fds <- load.fds(subj)

vname <- vnames[1]

# Particular trial
onsets1   <- vid.timing$onset[vid.timing$video == vname]
durs1     <- vid.timing$duration[vid.timing$video == vname]
cur.trial <- convolve.spm.hrf(onsets1, durs1, tot.time, runs)

# All other trials
onsets2 <- vid.timing$onset[vid.timing$video != vname]
durs2   <- vid.timing$duration[vid.timing$video != vname]
other.trials <- convolve.spm.hrf(onsets2, durs2, tot.time, runs)

# Run regression
dmat <- model.matrix(~cur.trial[,1:2] + other.trials[,1:2] + quests + days + repeat.faces + mc + fds)
ret.fit <- simple_lm(rdats, dmat)

# Do the regular one
sub.bs1 <- laply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  rdats <- sapply(dat$fmri$dat, rowMeans)
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  ## Covars (same throughout)
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.hrf(onsets3, durs3, tot.time, runs)
  # Day of scan
  days <- factor(runs)
  days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  # Repeat days (on the second and fourth day we repeat the stimuli)
  repeat.day <- factor(vid.timing$run)
  repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  onsets4    <- vid.timing$onset[repeat.day==2]
  durs4      <- vid.timing$duration[repeat.day==2]
  repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
  # Motion
  mc <- load.mc(subj)
  fds <- load.fds(subj)
  
  # Loop through each video and process
  system.time(vid.bs <- laply(vnames, function(vname) {
    # Particular trial
    onsets1   <- vid.timing$onset[vid.timing$video == vname]
    durs1     <- vid.timing$duration[vid.timing$video == vname]
    cur.trial <- convolve.hrf(onsets1, durs1, tot.time, runs)
    
    # All other trials
    onsets2 <- vid.timing$onset[vid.timing$video != vname]
    durs2   <- vid.timing$duration[vid.timing$video != vname]
    other.trials <- convolve.hrf(onsets2, durs2, tot.time, runs)
    
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


sub.bs2 <- laply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  rdats <- sapply(dat$fmri$dat, rowMeans)
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  ## Covars (same throughout)
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.hrf(onsets3, durs3, tot.time, runs)
  # Day of scan
  days <- factor(runs)
  days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  # Repeat days (on the second and fourth day we repeat the stimuli)
  repeat.day <- factor(vid.timing$run)
  repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  onsets4    <- vid.timing$onset[repeat.day==2]
  durs4      <- vid.timing$duration[repeat.day==2]
  repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
  # Motion
  mc <- load.mc(subj)
  fds <- load.fds(subj)
  
  # Loop through each video and process
  system.time(vid.bs <- laply(vnames, function(vname) {
    # Particular trial
    onsets1   <- vid.timing$onset[vid.timing$video == vname]
    durs1     <- vid.timing$duration[vid.timing$video == vname]
    cur.trial <- convolve.spm.hrf(onsets1, durs1, tot.time, runs, nparams=2, frate=1/12)[,1:2]
    
    # All other trials
    onsets2 <- vid.timing$onset[vid.timing$video != vname]
    durs2   <- vid.timing$duration[vid.timing$video != vname]
    other.trials <- convolve.spm.hrf(onsets2, durs2, tot.time, runs, nparams=2, frate=1/12)[,1:2]
    
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


## this returns the derivative
sub.bs3 <- laply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  rdats <- sapply(dat$fmri$dat, rowMeans)
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  ## Covars (same throughout)
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.hrf(onsets3, durs3, tot.time, runs)
  # Day of scan
  days <- factor(runs)
  days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  # Repeat days (on the second and fourth day we repeat the stimuli)
  repeat.day <- factor(vid.timing$run)
  repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  onsets4    <- vid.timing$onset[repeat.day==2]
  durs4      <- vid.timing$duration[repeat.day==2]
  repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
  # Motion
  mc <- load.mc(subj)
  fds <- load.fds(subj)
  
  # Loop through each video and process
  system.time(vid.bs <- laply(vnames, function(vname) {
    # Particular trial
    onsets1   <- vid.timing$onset[vid.timing$video == vname]
    durs1     <- vid.timing$duration[vid.timing$video == vname]
    cur.trial <- convolve.spm.hrf(onsets1, durs1, tot.time, runs, nparams=2, frate=1/12)[,1:2]
    
    # All other trials
    onsets2 <- vid.timing$onset[vid.timing$video != vname]
    durs2   <- vid.timing$duration[vid.timing$video != vname]
    other.trials <- convolve.spm.hrf(onsets2, durs2, tot.time, runs, nparams=2, frate=1/12)[,1:2]
    
    # Run regression
    dmat <- model.matrix(~cur.trial + other.trials + quests + days + repeat.faces + mc + fds)
    ret.fit <- simple_lm(rdats, dmat)
    
    # Get the betas and tvals for the current trial estimate
    bs <- ret.fit$b[3,]
    ts <- ret.fit$tvals[3,]
    
    cbind(betas=bs, tvals=ts)
  }, .parallel=T))
  dimnames(vid.bs)[[1]] <- vnames
  names(dimnames(vid.bs)) <- c("vid", "roi", "measure")
  
  return(vid.bs)
}, .progress="text")


# Compare correlation btw subjects for the two
betas.comp <- sapply(1:11, function(i) {
  mean((cor(t(sub.bs1[,,i,1])) - cor(t(sub.bs2[,,i,1])))[lower.tri(matrix(0,6,6))])
})
table(sign(betas.comp))
tvals.comp <- sapply(1:11, function(i) {
  mean((cor(t(sub.bs1[,,i,1])) - cor(t(sub.bs2[,,i,1])))[lower.tri(matrix(0,6,6))])
})
table(sign(tvals.comp))
## compare the tvals...it's always better with the 2param model
table(sign(colMeans(sapply(1:11, function(i) rowMeans(sub.bs1[,,i,2] - sub.bs2[,,i,2])))))



# Regressions -------------------------------------------------------------

# Get group average
grp.bs1 <- apply(sub.bs1, 2:4, mean)
grp.bs2 <- apply(sub.bs2, 2:4, mean)
grp.bs3 <- apply(sub.bs3, 2:4, mean)

# Prep the factors to regress
vnames <- dimnames(sub.bs2)$vid
inds <- sapply(vnames, function(vname) which(demo.vnames==vname))
fac.scores2 <- fac.scores[inds,]
fac.preds2  <- fac.preds[inds,]
fac.resids2 <- fac.resids[inds,]
pca.feats2 <- pca.face.feats[inds,]


# Do a regression
fit1t <- simple_lm(grp.bs1[,,2], fac.scores2)
fit2t <- simple_lm(grp.bs2[,,2], fac.scores2)
fit1b <- simple_lm(grp.bs1[,,1], fac.scores2)
fit2b <- simple_lm(grp.bs2[,,1], fac.scores2)


fit1t <- lm(grp.bs1[,,2] ~ fac.scores2)
fit2t <- lm(grp.bs2[,,2] ~ fac.scores2)
fit3t <- lm(grp.bs3[,,2] ~ fac.scores2)
fit1b <- lm(grp.bs1[,,1] ~ fac.scores2)
fit2b <- lm(grp.bs2[,,1] ~ fac.scores2)
fit3b <- lm(grp.bs3[,,1] ~ fac.scores2)
sfit1t <- summary(fit1t)
sfit2t <- summary(fit2t)
sfit3t <- summary(fit3t)
sfit1b <- summary(fit1b)
sfit2b <- summary(fit2b)
sfit3b <- summary(fit3b)

## oh full model fit is better for the regular HRF!!!
table(sign(sapply(sfit1t, function(x) x$adj.r.squared) - sapply(sfit2t, function(x) x$adj.r.squared)))
table(sign(sapply(sfit1b, function(x) x$adj.r.squared) - sapply(sfit2b, function(x) x$adj.r.squared)))
## compare beta to tval....tval wins!
table(sign(sapply(sfit1b, function(x) x$adj.r.squared) - sapply(sfit1t, function(x) x$adj.r.squared)))

round(sapply(sfit1t, function(x) x$adj.r.squared), 3)


# okay so now with the tvals try to do an aov to compare models
bs.tvals <- grp.bs1[,,2]
fit <- aov(bs.tvals ~ fac.preds2 + fac.resids2 + pca.feats2)
sfit <- summary(fit)
sfit[[3]]
fit <- lm(bs.tvals ~ fac.preds2 + fac.resids2 + pca.feats2)
summary(fit)[[3]]

for (i in 1:6) {
  cat(i,"\n")
  bs.tvals <- sub.bs1[i,,,1]
  print(summary(aov(bs.tvals ~ fac.preds2 + fac.resids2 + pca.feats2))[[3]])
}




# Classify ----------------------------------------------------------------

library(glmnet)

get_bestfit <- function(cvfit, type.measure, exclude.zero=FALSE) {
  bestouts <- cvfit$measures[[type.measure]]
  extreme  <- ifelse(type.measure == "rmse", min, max)
  
  if (exclude.zero) {
    val <- extreme(bestouts[cvfit$nzero>0])
  } else {
    val <- extreme(bestouts)
  }
  
  ind <- which(bestouts == val)
  
  bestfit  <- list(
    measure = type.measure, 
    val     = val, 
    ind     = ind, 
    lam     = cvfit$lambda[ind], 
    preval  = cvfit$fit.preval[,ind], 
    nzero   = cvfit$nzero[ind], 
    coef    = coef(cvfit, s=cvfit$lambda[ind])
  )
  bestfit
}

run_cvglmnet <- function(X, y, keep=T, parallel=T, type.measure="rsq", exclude.zero=FALSE, ...) 
{
  if (!(type.measure %in% c("rsq", "r", "rmse"))) stop("unknown type.measure: ", type.measure)
  
  rmse <- function(x1, x2) sqrt(mean((x1-x2)^2))
  
  #cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel)
  cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel, ...)
  
  rs     <- sapply(1:length(cvfit$lambda), function(i) cor(cvfit$fit.preval[,i], y))
  rsqs   <- rs^2
  rmses  <- sapply(1:length(cvfit$lambda), function(i) rmse(cvfit$fit.preval[,i], y))
  cvfit$measures <- list(
    r = rs, 
    rsq = rsqs, 
    rmse = rmses
  )
  
  cvfit$bestfit <- get_bestfit(cvfit, type.measure, exclude.zero)
  
  return(cvfit)
}

i <- 1
bs.tvals <- grp.bs1[,,2]
for (i in 1:6) {
  bs.tvals <- sub.bs2[i,,,2]
  cvfit1 <- run_cvglmnet(fac.preds2, bs.tvals[,3])
  print(cvfit1$bestfit$val)
}

for (i in 1:6) {
  bs.tvals <- sub.bs1[i,,,2]
  cvfit1 <- run_cvglmnet(fac.preds2, bs.tvals[,3])
  print(cvfit1$bestfit$val)
}

# oh wait this btw subj works!
tmp <- unlist(lapply(1:6, function(i) sub.bs1[i,,3,2]))
subs <- rep(1:6, each=864)
cvfit1 <- run_cvglmnet(fac.preds2[rep(1:nrow(fac.preds2),6),], tmp, foldid=subs)
cvfit1$bestfit$val

cvfit2 <- run_cvglmnet(fac.resids2[rep(1:nrow(fac.resids2),6),], tmp, foldid=subs)
cvfit2$bestfit$val

cvfit3 <- run_cvglmnet(pca.feats2[rep(1:nrow(pca.feats2),6),], tmp, foldid=subs)
cvfit3$bestfit$val

X <- cbind(fac.preds2, fac.resids2, pca.feats2)[rep(1:nrow(pca.feats2),6),]
cvfit <- run_cvglmnet(X, tmp, foldid=subs)
cvfit$bestfit$val

X <- cbind(fac.resids2, pca.feats2)[rep(1:nrow(pca.feats2),6),]
cvfit <- run_cvglmnet(X, tmp, foldid=subs)
cvfit$bestfit$val

X <- cbind(fac.preds2, pca.feats2)[rep(1:nrow(pca.feats2),6),]
cvfit <- run_cvglmnet(X, tmp, foldid=subs)
cvfit$bestfit$val

X <- cbind(fac.preds2, fac.resids2)[rep(1:nrow(pca.feats2),6),]
cvfit <- run_cvglmnet(X, tmp, foldid=subs, alpha=0)
cvfit$bestfit$val
