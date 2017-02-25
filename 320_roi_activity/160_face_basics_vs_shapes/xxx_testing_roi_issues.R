# See the results with our ROIs

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(niftir)


# ROIs --------------------------------------------------------------------

# ROI Information
roi.names <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", "r.Amyg", 
               "l.vATL", "l.FFA", "l.OFA", "l.EBA", "l.Amyg")
nrois <- length(roi.names)
roi.names1 <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", 
                "l.vATL", "l.FFA", "l.OFA", "l.EBA")
nrois1 <- length(roi.names1)
roi.names2 <- c("r.Amyg", "l.Amyg")
nrois2 <- length(roi.names2)

mrois1 <- read.mask("/data1/famface01/analysis/roi/Functional_v3/group_rois.nii.gz", NULL)
mask1 <- mrois1 != 0

mrois2 <- read.mask("/data1/famface01/analysis/roi/Functional_v3/amygdala/group_rois_select.nii.gz", NULL)
mask2 <- mrois2 != 0

mrois <- mrois1 * 0
for (i in 1:nrois) {
  if (any(roi.names[i] == roi.names1)) {
    ri <- which(roi.names[i] == roi.names1)
    mrois[mrois1 == ri] <- i
  } else if (any(roi.names[i] == roi.names2)) {
    ri <- which(roi.names[i] == roi.names2)
    mrois[mrois2 == ri] <- i
  } else {
    stop("error")
  }
}
rnames <- roi.names
mask <- mrois!=0
rois <- mrois[mask]
roi.rnames <- factor(rois, levels=1:length(rnames), labels=rnames)



# Load Data ---------------------------------------------------------------

# Cluster corrected
indir <- "/data1/famface01/analysis/task_activity/group/face_basics_unfam/traitsfa_givenshape.reml/easythresh"
infiles <- list.files(indir, pattern="combined_thresh_zstat.*trait")
innames <- sub("_", ".", sub("[.]nii[.]gz", "", sub("combined_thresh_zstat_", "", infiles)))
voxdats <- sapply(infiles, function(ifile) {
  read.mask(file.path(indir, ifile), NULL)[mask]
})
colnames(voxdats) <- innames

# Raw
indir <- "/data1/famface01/analysis/task_activity/group/face_basics_unfam/traitsfa_givenshape.reml"
infiles <- list.files(indir, pattern="^zstat.*trait")
voxdats2 <- sapply(infiles, function(ifile) {
  read.mask(file.path(indir, ifile), NULL)[mask]
})
colnames(voxdats2) <- innames


# Check Voxelwise ---------------------------------------------------------

# cluster corrected
ret <- apply(voxdats!=0, 2, function(x) tapply(x, roi.rnames, mean))
colnames(ret) <- c(sprintf("resid%02i", 1:6), sprintf("shape%02i", 1:6))
barplot(t(ret), beside = T, legend=T)
round(ret, 3)*100

# not cluster corrected
ret2 <- apply(abs(voxdats2)>1.96, 2, function(x) tapply(x, roi.rnames, mean))
colnames(ret2) <- c(sprintf("resid%02i", 1:6), sprintf("shape%02i", 1:6))
#barplot(t(ret2), beside = T, legend=T)
round(ret2, 3)*100


# So the above results are nice because they show that only for trait 5 and 6
# which are about typicality and attractiveness; and to some extent the memorable
# dimension but that's only in the left hemisphere. Nothing is good with #2.
# 
# TODO: get photos of each...
# 
# The residuals are extensive for 1,2,and4.


# Compare with ROIs ---------------------------------------------------------

# What we want to see are ROI effects similar to with AFNI.

subjects <- sprintf("sub%02i", 1:6)

# Load the voxelwise ROI data
load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects

# Average in each ROI
lst.rdats <- lapply(dat.vols, function(dat) {
  rdats <- sapply(dat$fmri$dat, rowMeans)
  rdats
})
rdats <- do.call(rbind, lst.rdats)

# Load the XMAT information
xmat_labs <- function(fn) {
  str <- system(sprintf("grep ColumnLabels %s | sed s/'#  ColumnLabels = '//", fn), intern=T)
  str <- gsub("\"", "", str)
  cols <- strsplit(str, ' ; ')[[1]]
  cols
}
xfiles <- sprintf("/data1/famface01/analysis/task_activity/%s/face_basics_unfam/traitsfa_givenshape.reml/xmat.1D", subjects)
lst.xmats <- lapply(xfiles, function(xfile) {
  xmat <- read.table(xfile)
  colnames(xmat) <- xmat_labs(xfile)
  colnames(xmat) <- sub("#0", "", colnames(xmat))
  xmat
})
xmats  <- do.call(rbind, lst.xmats)
xmats2 <- xmats[,-c(1:16)]
ssubs <- rep(subjects, sapply(lst.xmats, nrow))

# Run the regressions
fit <- lm(rdats ~ ssubs + ., data=xmats[,-1])
sfit <- summary(fit)
sfit[[4]]

## try without the runs
fit <- lm(rdats ~ ssubs + ., data=xmats2)
sfit <- summary(fit)
sfit[[4]]

runs   <- as.matrix(xmats[,1:16])
shapes <- as.matrix(xmats2[,3:8])
resids <- as.matrix(xmats2[,9:14])
motion <- as.matrix(xmats2[,15:20])
fit <- aov(rdats ~ ssubs + runs + faces + quests + shapes + resids + motion, data=xmats2)
sfit <- summary(fit)
sfit[[3]]



# Convolve on Own ---------------------------------------------------------

# So we now must do the convolutions on our own. Als

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

mc.sub01 <- load.mc("sub01")

mc2.sub01 <- xmats2[ssubs=="sub01",15:20]



# Let's check the face signal

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

subj <- "sub01"
dat <- dat.vols$sub01
rdats <- sapply(dat$fmri$dat, rowMeans)

vid.timing <- dat$basics$timing
vnames <- levels(vid.timing$video)

runs <- dat$fmri$runs
tot.time <- length(runs)

# All faces
onsets0 <- vid.timing$onset
durs0   <- vid.timing$duration
faces.all <- convolve.hrf(onsets0, durs0, tot.time, runs)
faces.all2 <- convolve.spm.hrf(onsets0, durs0, tot.time, runs, nparams=1)


afni.faces <- xmats$faces[ssubs=="sub01"]

plot.ts(cbind(faces.all, afni.faces, faces.all2*2.5)[1:100,], plot.type="single", col=1:3)

faces.allx <- faces.all*2.5
summary(lm(rdats[,3] ~ faces.all + afni.faces))




# Compare motion ----------------------------------------------------------

# We want to see using our updated motion parameters vs the truncated ones used
# in the voxelwise analyses

### Let's add the correct motion information and see what we get
# Average in each ROI
lst.rdats <- lapply(dat.vols, function(dat) {
  rdats <- sapply(dat$fmri$dat, rowMeans)
  rdats
})
rdats <- do.call(rbind, lst.rdats)

lst.mcs <- lapply(subjects, load.mc)
lst.mcs2 <- lapply(lst.mcs, scale, center=T, scale=F)
names(lst.mcs) <- subjects
mcs <- do.call(rbind, lst.mcs)
mcs2 <- do.call(rbind, lst.mcs2)

# Plots
plot.ts(xmats2[1:5000,15:20])
plot.ts(mcs[1:5000,])

# Original
runs   <- as.matrix(xmats[,1:16])
shapes <- as.matrix(xmats2[,3:8])
resids <- as.matrix(xmats2[,9:14])
motion <- as.matrix(xmats2[,15:20])
fit <- aov(rdats ~ ssubs + runs + faces + quests + shapes + resids + motion, data=xmats2)
sfit <- summary(fit)
sfit[[3]]

# Update the motion
fit <- aov(rdats ~ ssubs + runs + faces + quests + shapes + resids + mcs2, data=xmats2)
sfit <- summary(fit)
sfit[[3]]



# Compare convolutions -----------------------------------------------------

# I want to see now if using a different convolution approach would change anything.

# Get the trait information
traitsfa <- read.csv("/data1/famface01/command/misc/face_representations/300_task_activity/150_face_basics_unfam/z_traitsfa_givenshape.csv")
tvids <- as.character(traitsfa$X)
traitsfa <- traitsfa[,-1]
shape.traits <- traitsfa[,grep("shape", colnames(traitsfa))]
resid.traits <- traitsfa[,grep("resid", colnames(traitsfa))]

vnames <- as.character(dat.vols$sub01$basics$timing$video)
head(vnames)
match.inds <- sapply(vnames, function(vname) which(tvids==vname))
all.equal(vnames, tvids[match.inds])
traitsfa2 <- traitsfa[match.inds,] # this is double the length due to repeats


# Do the regular one
reg.xmats <- llply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  rdats <- sapply(dat$fmri$dat, rowMeans)
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  # Faces
  onsets0 <- vid.timing$onset
  durs0   <- vid.timing$duration
  faces   <- convolve.hrf(onsets0, durs0, tot.time, runs)
  
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.hrf(onsets3, durs3, tot.time, runs)
  
  ## Day of scan
  #days <- factor(runs)
  #days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  #days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  
  ## Repeat days (on the second and fourth day we repeat the stimuli)
  #repeat.day <- factor(vid.timing$run)
  #repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  #onsets4    <- vid.timing$onset[repeat.day==2]
  #durs4      <- vid.timing$duration[repeat.day==2]
  #repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
  
  # Motion
  mc <- load.mc(subj)
  #fds <- load.fds(subj)
  
  # Traits
  vnames <- as.character(dat$basics$timing$video)
  match.inds <- sapply(vnames, function(vname) which(tvids==vname))
  shape.traits2 <- shape.traits[match.inds,]
  resid.traits2 <- resid.traits[match.inds,]
  shape.traits3 <- convolve.hrf(onsets0, durs0, tot.time, runs, shape.traits2)
  resid.traits3 <- convolve.hrf(onsets0, durs0, tot.time, runs, resid.traits2)
  
  ret <- cbind(faces=faces, quests=quests, mc=mc, shape.traits=shape.traits3, 
                resid.traits=resid.traits3)
  colnames(ret) <- c("faces", "quests", colnames(mc), sprintf("shape%02i", 1:ncol(shape.traits3)), sprintf("resids%02i", 1:ncol(resid.traits3)))
  ret
}, .progress="text")
reg.xmats2 <- do.call(rbind, reg.xmats)


runs   <- as.matrix(xmats[,1:16])
shapes <- as.matrix(xmats2[,3:8])
resids <- as.matrix(xmats2[,9:14])
fit <- aov(rdats ~ ssubs + runs + faces + quests + shapes + resids + mcs, data=xmats2)
sfit <- summary(fit)
sfit[[3]]

shapes <- reg.xmats2[,grep("shape", colnames(reg.xmats2))]
resids <- reg.xmats2[,grep("resid", colnames(reg.xmats2))]
fit <- aov(rdats ~ ssubs + runs + faces + quests + shapes + resids + mcs, data=as.data.frame(reg.xmats2))
sfit <- summary(fit)
sfit[[3]]



spm.xmats <- llply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  rdats <- sapply(dat$fmri$dat, rowMeans)
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  # Faces
  onsets0 <- vid.timing$onset
  durs0   <- vid.timing$duration
  faces   <- convolve.spm.hrf(onsets0, durs0, tot.time, runs, nparams=2)
  
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.spm.hrf(onsets3, durs3, tot.time, runs, nparams=2)
  
  ## Day of scan
  #days <- factor(runs)
  #days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  #days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  
  ## Repeat days (on the second and fourth day we repeat the stimuli)
  #repeat.day <- factor(vid.timing$run)
  #repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  #onsets4    <- vid.timing$onset[repeat.day==2]
  #durs4      <- vid.timing$duration[repeat.day==2]
  #repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
  
  # Motion
  mc <- load.mc(subj)
  #fds <- load.fds(subj)
  
  # Traits
  vnames <- as.character(dat$basics$timing$video)
  match.inds <- sapply(vnames, function(vname) which(tvids==vname))
  shape.traits2 <- shape.traits[match.inds,]
  resid.traits2 <- resid.traits[match.inds,]
  tmp <- lapply(1:ncol(shape.traits2), function(i) convolve.spm.hrf(onsets0, durs0, tot.time, runs, shape.traits2[,i], nparams = 2))
  shape.traits3 <- do.call(cbind, tmp)
  tmp <- lapply(1:ncol(resid.traits2), function(i) convolve.spm.hrf(onsets0, durs0, tot.time, runs, resid.traits2[,i], nparams = 2))
  resid.traits3 <- do.call(cbind, tmp)
  
  ret <- cbind(faces, quests, mc, shape.traits3, resid.traits3)
  colnames(ret) <- c(rep("faces",2), rep("quests",2), colnames(mc), rep(sprintf("shape%02i", 1:ncol(shape.traits2)),2), rep(sprintf("resids%02i", 1:ncol(resid.traits2)),2))
  ret
}, .progress="text")
spm.xmats2 <- do.call(rbind, spm.xmats)

shapes <- reg.xmats2[,grep("shape", colnames(reg.xmats2))]
resids <- reg.xmats2[,grep("resid", colnames(reg.xmats2))]
fit <- aov(rdats ~ ssubs + runs + faces + quests + shapes + resids + mcs, data=as.data.frame(reg.xmats2))
sfit <- summary(fit)
sfit[[3]]

faces <- spm.xmats2[,grep("faces", colnames(spm.xmats2))]
faces1 <- faces[,1]
faces2 <- faces[,2]
quests <- spm.xmats2[,grep("quests", colnames(spm.xmats2))]
quests1 <- quests[,1]
quests2 <- quests[,2]
shapes <- spm.xmats2[,grep("shape", colnames(spm.xmats2))]
shapes1 <- shapes[,seq(1,12,by=2)]
colnames(shapes1) <- sprintf("shape%02i", 1:ncol(shapes1))
shapes2 <- shapes[,seq(2,12,by=2)]
resids <- spm.xmats2[,grep("resid", colnames(spm.xmats2))]
resids1 <- resids[,seq(1,12,by=2)]
colnames(resids1) <- sprintf("shape%02i", 1:ncol(resids1))
resids2 <- resids[,seq(2,12,by=2)]
fit2 <- aov(rdats ~ ssubs + runs + faces1 + faces2 + quests1 + quests2 + shapes1 + shapes2 + resids1 + resids2 + mcs)
sfit <- summary(fit2)
sfit[[3]]

fit1 <- aov(rdats ~ ssubs + runs + faces1 + quests1 + shapes1 + resids1 + mcs)
sfit <- summary(fit1)
sfit[[3]]

anova(fit1, fit2)

# regression
fit1 <- lm(rdats[,3] ~ ssubs + runs + faces1 + quests1 + shapes1 + resids1 + mcs)
summary(fit1)



# Compare Added Factors ---------------------------------------------------

# We see the effect of a model with a larger factor model
# Features
indir <- "/data1/famface01/analysis/misc/320_roi_task_activity"
load(file.path(indir, "20_predict_face_feats.rda"), verbose=T)

# Convolve the fac.preds and fac.resids
spm.xmats.more <- llply(subjects, function(subj) {
  dat <- dat.vols[[subj]]
  
  rdats <- sapply(dat$fmri$dat, rowMeans)
  
  runs <- dat$fmri$runs
  tot.time <- length(runs)
  
  vid.timing <- dat$basics$timing
  vnames <- sort(levels(vid.timing$video))
  
  # Faces
  onsets0 <- vid.timing$onset
  durs0   <- vid.timing$duration
  faces   <- convolve.spm.hrf(onsets0, durs0, tot.time, runs, nparams=2)
  
  # Questions
  onsets3 <- vid.timing$onset[vid.timing$question!="none"]
  durs3   <- 4
  quests  <- convolve.spm.hrf(onsets3, durs3, tot.time, runs, nparams=2)
  
  ## Day of scan
  #days <- factor(runs)
  #days <- mapvalues(days, from=1:16, to=rep(1:4,each=4))
  #days <- model.matrix(~days)[,-1] # matrix of 1s 0s ... no intercept
  
  ## Repeat days (on the second and fourth day we repeat the stimuli)
  #repeat.day <- factor(vid.timing$run)
  #repeat.day <- mapvalues(repeat.day, from=1:16, to=rep(c(1,2,1,2),each=4))
  #onsets4    <- vid.timing$onset[repeat.day==2]
  #durs4      <- vid.timing$duration[repeat.day==2]
  #repeat.faces <- convolve.hrf(onsets4, durs4, tot.time, runs)
  
  # Motion
  mc <- load.mc(subj)
  #fds <- load.fds(subj)
  
  # Traits
  vnames <- as.character(dat$basics$timing$video)
  match.inds <- sapply(vnames, function(vname) which(demo.vnames==vname))
  shape.traits2 <- fac.preds[match.inds,]
  resid.traits2 <- fac.resids[match.inds,]
  tmp <- lapply(1:ncol(shape.traits2), function(i) convolve.spm.hrf(onsets0, durs0, tot.time, runs, shape.traits2[,i], nparams = 2))
  shape.traits3 <- do.call(cbind, tmp)
  tmp <- lapply(1:ncol(resid.traits2), function(i) convolve.spm.hrf(onsets0, durs0, tot.time, runs, resid.traits2[,i], nparams = 2))
  resid.traits3 <- do.call(cbind, tmp)
  
  ret <- cbind(faces, quests, mc, shape.traits3, resid.traits3)
  colnames(ret) <- c(rep("faces",2), rep("quests",2), colnames(mc), rep(sprintf("shape%02i", 1:ncol(shape.traits2)),2), rep(sprintf("resids%02i", 1:ncol(resid.traits2)),2))
  ret
}, .progress="text")
spm.xmats.more2 <- do.call(rbind, spm.xmats.more)


# Do full model with derivative
faces <- spm.xmats.more2[,grep("faces", colnames(spm.xmats.more2))]
faces1 <- faces[,1]
faces2 <- faces[,2]
quests <- spm.xmats.more2[,grep("quests", colnames(spm.xmats.more2))]
quests1 <- quests[,1]
quests2 <- quests[,2]
shapes <- spm.xmats.more2[,grep("shape", colnames(spm.xmats.more2))]
shapes1 <- shapes[,seq(1,36,by=2)]
colnames(shapes1) <- sprintf("shape%02i", 1:ncol(shapes1))
shapes2 <- shapes[,seq(2,36,by=2)]
resids <- spm.xmats.more2[,grep("resid", colnames(spm.xmats.more2))]
resids1 <- resids[,seq(1,36,by=2)]
colnames(resids1) <- sprintf("resids%02i", 1:ncol(resids1))
resids2 <- resids[,seq(2,36,by=2)]
fit2 <- aov(rdats ~ ssubs + runs + faces1 + faces2 + quests1 + quests2 + shapes1 + shapes2 + resids1 + resids2 + mcs)
sfit <- summary(fit2)
sfit[[3]]

# No derivative
fit1 <- aov(rdats ~ ssubs + runs + faces1 + quests1 + shapes1 + resids1 + mcs)
sfit <- summary(fit1)
sfit[[3]]

# Let's run the regular regression to see the fits
fit1 <- lm(rdats[,3] ~ ssubs + runs + faces1 + quests1 + shapes1 + resids1 + mcs)
summary(fit1)

# Try anova on subset
summary(aov(rdats[,3] ~ ssubs + runs + faces1 + quests1 + shapes1[,c(1:5,9)] + resids1[,c(1:5,9)] + mcs))

# Try regression on the amygdala
summary(lm(rdats[,6] ~ ssubs + runs + faces1 + quests1 + shapes1[,c(1:5,9)] + resids1[,c(1:5,9)] + mcs))
summary(lm(rdats[,11] ~ ssubs + runs + faces1 + quests1 + shapes1[,c(1:5,9)] + resids1[,c(1:5,9)] + mcs))
rnames[6]

# SIDE SHOW
# Can we try the cancor approach?
library(PMA)
perm.out <- CCA.permute(rdats, cbind(faces1,quests1,shapes1[,c(1:5,9)]), typex="standard", typez="standard", nperms=25)
out <- CCA(rdats, cbind(faces1,quests1,shapes1[,c(1:5,9)]), typex="standard",typez="standard", K=6, penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)
rownames(out$u) <- rnames
rownames(out$v) <- c("face", "quest", sprintf("shape%0i", 1:6))
round(out$u, 3)
round(out$v, 3)




# View All ROI Results ----------------------------------------------------

fits <- aov(rdats ~ ssubs + runs + faces1 + quests1 + shapes1[,c(1:5,9)] + resids1[,c(1:5,9)] + mcs)
sfits <- summary(fits)

# Get the fstats as tstats
tstats <- sapply(sfits, function(sfit) {
  fstats <- sfit$`F value`
  fstats <- fstats[-length(fstats)]
  tstats <- sqrt(fstats)
  names(tstats) <- c("sub", "run", "face", "quest", "shapes", "resids", "motion")
  tstats
})

# Format the roi names and add as colnames
rnames2 <- sub("[.]", " ", rnames)
rnames2 <- sub("^r", "R", rnames2)
rnames2 <- sub("^l", "L", rnames2)
colnames(tstats) <- rnames2

# Remove sub, run, and quest
tstats <- tstats[-c(1:2,4),]

# Add a gap between the left and the right
tstats2 <- cbind(tstats[,1:6], 0, tstats[,7:11])

# Plot with the faces
col <- brewer.pal(4, "Spectral")
barplot(tstats2, beside=T, legend=T, col=col, ylab=bquote(Tstat==sqrt(Fstat)))
abline(h=1.65, lty=3)

# Plot without faces
barplot(tstats2[-1,], beside=T, legend=T, col=col[-1], ylab=bquote(Tstat==sqrt(Fstat)))
abline(h=1.65, lty=3)



# Hierarchical Classification  ---------------------------------------------

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


# Setup the basics for the classification
snums <- as.numeric(factor(ssubs))
runs  <- dat.vols$sub01$fmri$runs
sshapes1 <- shapes1[,c(1:5,9)]
sresids1 <- resids1[,c(1:5,9)]
sshapes2 <- shapes2[,c(1:5,9)]
sresids2 <- resids2[,c(1:5,9)]


# Let's first get the reduced vs full model results for each average ROI
# for each subject
modelr2.rets <- laply(subjects, function(subj) {
  X1 <- model.matrix(~faces1 + quests1 + sshapes1 + sresids1 + mcs)
  X2 <- model.matrix(~faces1 + quests1 + sresids1 + mcs)
  X3 <- model.matrix(~quests1 + sshapes1 + sresids1 + mcs)
  
  rets <- sapply(1:ncol(rdats), function(i) {
    cvfit1 <- run_cvglmnet(X1[ssubs==subj,], rdats[ssubs==subj,i], foldid=runs, parallel=F)
    cvfit2 <- run_cvglmnet(X2[ssubs==subj,], rdats[ssubs==subj,i], foldid=runs, parallel=F)
    cvfit3 <- run_cvglmnet(X3[ssubs==subj,], rdats[ssubs==subj,i], foldid=runs, parallel=F)
    c(full=cvfit1$bestfit$val, reduced=cvfit2$bestfit$val, face.reduced=cvfit3$bestfit$val)
  })
  
  rets
}, .parallel=T)
dimnames(modelr2.rets)[[1]] <- subjects
dimnames(modelr2.rets)[[3]] <- rnames2
# a negative value means that the reduced model is better
# a positive value indicates the amount of additional information from having that variable
modelr2.full       <- modelr2.rets[,1,]
modelr2.diff.shape <- modelr2.rets[,1,] - modelr2.rets[,2,]
modelr2.diff.face  <- modelr2.rets[,1,] - modelr2.rets[,3,] 

# Plot
barplot(colMeans(modelr2.full), ylab="Lasso r-squared", main="Predicting ROI data with face+factors+etc model (Full Model)")
barplot(colMeans(modelr2.diff.shape), ylab="Full - Reduced Model (r-squared)", main="Remove Predicted Factor Scores", sub="(Higher is Better)")
abline(h=0)
barplot(colMeans(modelr2.diff.face), ylab="Full - Reduced Model (r-squared)", main="Remove Face Regressor", sub="(Higher is Better)")
abline(h=0)
## also plot the # of subjects that show a positive effect
barplot(colMeans(sign(modelr2.diff.shape))*100, ylab="% of Subjects with Full - Reduced Model (r-squared)", main="Remove Predicted Factor Scores", sub="(Higher is Better)")
abline(h=0)




# Voxelwise ---------------------------------------------------------------

# Here we test out running the analyses voxelwise and compare those results to
# the average ROI as well as the results from AFNI


#subj <- "sub01"

ave.sigs <- laply(subjects, function(subj) {
  cat(subj, "\n")
  sdat <- dat.vols[[subj]]
  laply(sdat$fmri$dat, function(sroi) {
    sroi <- as.matrix(sroi)
    
    X1 <- model.matrix(~faces1 + quests1 + sshapes1 + sresids1 + mcs)
    inds <- ssubs==subj
    
    fits <- aov(sroi ~ faces1[inds] + quests1[inds] + sshapes1[inds,] + sresids1[inds,] + mcs[inds,])
    sfits <- summary(fits)
    
    ststats <- sapply(sfits, function(sfit) {
      fstats <- sfit$`F value`
      fstats <- fstats[-length(fstats)]
      tstats <- sqrt(fstats)
      names(tstats) <- c("face", "quests", "shapes", "resids", "motion")
      tstats
    })
    
    rowMeans(ststats>1.65)
  }, .parallel=T)
})
dimnames(ave.sigs)[[1]] <- subjects
dimnames(ave.sigs)[[2]] <- rnames2
round(ave.sigs, 3)*100

round(apply(ave.sigs, 2:3, mean), 3)*100



ave.sigs2 <- laply(subjects, function(subj) {
  cat(subj, "\n")
  sdat <- dat.vols[[subj]]
  laply(sdat$fmri$dat, function(sroi) {
    sroi <- as.matrix(sroi)
    
    X1 <- model.matrix(~faces1 + quests1 + sshapes1 + sresids1 + mcs)
    inds <- ssubs==subj
    
    fits <- lm(sroi ~ faces1[inds] + quests1[inds] + sshapes1[inds,] + sresids1[inds,] + mcs[inds,])
    sfits <- summary(fits)
    tstats <- sapply(sfits, function(sfit) sfit$coefficients[,3])
    
    rowMeans(tstats>1.96)
  }, .parallel=T)
})
dimnames(ave.sigs2)[[1]] <- subjects
dimnames(ave.sigs2)[[2]] <- rnames2
tmp <- ave.sigs2[,,4:15]
tmp <- apply(tmp, 2:3, mean)
colnames(tmp) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(tmp, 3)*100







# Load the XMAT information
xmat_labs <- function(fn) {
  str <- system(sprintf("grep ColumnLabels %s | sed s/'#  ColumnLabels = '//", fn), intern=T)
  str <- gsub("\"", "", str)
  cols <- strsplit(str, ' ; ')[[1]]
  cols
}
xfiles <- sprintf("/data1/famface01/analysis/task_activity/%s/face_basics_unfam/traitsfa_givenshape.reml/xmat.1D", subjects)
lst.xmats <- lapply(xfiles, function(xfile) {
  xmat <- read.table(xfile)
  colnames(xmat) <- xmat_labs(xfile)
  colnames(xmat) <- sub("#0", "", colnames(xmat))
  xmat
})
xmats  <- do.call(rbind, lst.xmats)
xmats2 <- xmats[,-c(1:16)]
ssubs <- rep(subjects, sapply(lst.xmats, nrow))


ave.sigs3 <- laply(subjects, function(subj) {
  cat(subj, "\n")
  sdat <- dat.vols[[subj]]
  laply(sdat$fmri$dat, function(sroi) {
    sroi <- as.matrix(sroi)
    
    fits <- lm(sroi ~ ., data=xmats[ssubs==subj,-1])
    sfits <- summary(fits)
    tstats <- sapply(sfits, function(sfit) sfit$coefficients[,3])
    
    rowMeans(abs(tstats)>1.96)
  }, .parallel=T)
})
dimnames(ave.sigs3)[[1]] <- subjects
dimnames(ave.sigs3)[[2]] <- rnames2
tmp <- ave.sigs3[,,19:30]
tmp <- apply(tmp, 2:3, mean)
colnames(tmp) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(tmp, 3)*100


# Let's plot the afni results again
ret <- apply(voxdats!=0, 2, function(x) tapply(x, roi.rnames, mean))
ret <- ret[,c(7:12,1:6)]
colnames(ret) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(ret, 3)*100


# Plot the regression of the average ROIs
fits <- lm(rdats ~ ssubs + ., data=xmats[,-1])
sfits <- summary(fits)
tvals <- sapply(sfits, function(sfit) sfit$coefficients[-c(1:21),3][3:14])
tvals <- t(tvals)
rownames(tvals) <- rnames2
colnames(tvals) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(tvals, 3)
round(tvals*(abs(tvals)>1.96), 2)
round(tvals*(abs(tvals)>1.96), 2)[,c(7:12,1:6)]

# Note try above with the proper motion
fits <- lm(rdats ~ ssubs + . + mcs, data=xmats[,-c(1,31:36)])
sfits <- summary(fits)
tvals <- sapply(sfits, function(sfit) sfit$coefficients[-c(1:21),3][3:14])
tvals <- t(tvals)
rownames(tvals) <- rnames2
colnames(tvals) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(tvals, 3)
round(tvals*(abs(tvals)>1.96), 2)
round(tvals*(abs(tvals)>1.96), 2)[,c(7:12,1:6)]


# Standard ROIs -----------------------------------------------------------

# We shall try to load the ROIs in standard space for each participant
 

# ROI Information
roi.names <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", "r.Amyg", 
               "l.vATL", "l.FFA", "l.OFA", "l.EBA", "l.Amyg")
nrois <- length(roi.names)
roi.names1 <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", 
                "l.vATL", "l.FFA", "l.OFA", "l.EBA")
nrois1 <- length(roi.names1)
roi.names2 <- c("r.Amyg", "l.Amyg")
nrois2 <- length(roi.names2)

mrois1 <- read.mask("/data1/famface01/analysis/roi/Functional_v3/group_rois.nii.gz", NULL)
mask1 <- mrois1 != 0

mrois2 <- read.mask("/data1/famface01/analysis/roi/Functional_v3/amygdala/group_rois_select.nii.gz", NULL)
mask2 <- mrois2 != 0

mrois <- mrois1 * 0
for (i in 1:nrois) {
  if (any(roi.names[i] == roi.names1)) {
    ri <- which(roi.names[i] == roi.names1)
    mrois[mrois1 == ri] <- i
  } else if (any(roi.names[i] == roi.names2)) {
    ri <- which(roi.names[i] == roi.names2)
    mrois[mrois2 == ri] <- i
  } else {
    stop("error")
  }
}
rnames <- roi.names
mask <- mrois!=0
rois <- mrois[mask]
roi.rnames <- factor(rois, levels=1:length(rnames), labels=rnames)


# Now that we have the rois, we need to load the ROI data
load.fmri <- function(subj, rois, mask) {
  preproc.base <- "/data1/famface01/analysis/preprocessed"
  
  # Masks
  rois      <- rois[mask]
  
  # Functional
  func.fnames <- Sys.glob(sprintf("%s/%s/func/unfam_vids/std_filtered_func_run*.nii.gz", preproc.base, subj))
  func.fnames <- sort(func.fnames)
  
  brain.dat <- ldply(func.fnames, function(funcfile) {
    func <- read.big.nifti(funcfile)
    func <- func[,mask]
    func <- scale(func, center=T, scale=F)
    func
  }, .progress="text")
  
  rois.dat <- llply(1:nrois, function(i) {
    brain.dat[,rois == i]
  }, .progress="text")
  rm(brain.dat)
  names(rois.dat) <- roi.names
  
  # also save which time-point is which run
  runs <- llply(1:length(func.fnames), function(i) {
    hdr <- read.nifti.header(func.fnames[1])
    rep(i, hdr$dim[4])
  }, .progress="text")
  runs <- unlist(runs)
  
  list(mask=mask, rois=rois, dat=rois.dat, runs=runs)
}
std.lst.rdats <- lapply(subjects, function(subj) {
  cat(subj, "\n")
  load.fmri(subj, mrois, mask)
})
names(std.lst.rdats) <- subjects
## collapse across subjects
std.rdats.all <- llply(rnames, function(rname) {
  ldat <- lapply(std.lst.rdats, function(s.rdats) {
    s.rdats$dat[[rname]]
  })
  do.call(rbind, ldat)
}, .progress="text")
names(std.rdats.all) <- rnames
## get average voxelwise data
std.rdats <- sapply(std.rdats.all, rowMeans)
## save for future
save(rnames, roi.rnames, mrois, mask, rois, std.rdats.all, std.rdats, 
     file="/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda")
## load
load("/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda", verbose=T)
std.rdats <- sapply(std.rdats.all, rowMeans)

# Ok so now let's try the regression approach again
fits2 <- lm(std.rdats ~ ssubs + . + mcs2, data=xmats[,-c(1,31:36)])
sfits2 <- summary(fits2)
tvals2 <- sapply(sfits2, function(sfit) sfit$coefficients[-c(1:21),3][3:14])
tvals2 <- t(tvals2)
rownames(tvals2) <- rnames2
colnames(tvals2) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(tvals2, 3)
round(tvals2*(abs(tvals2)>1.96), 2)
round(tvals2*(abs(tvals2)>1.96), 2)[,c(7:12,1:6)]

# as a comparison, our native space ROIs
fits <- lm(rdats ~ ssubs + . + mcs, data=xmats[,-c(1,31:36)])
sfits <- summary(fits)
tvals <- sapply(sfits, function(sfit) sfit$coefficients[-c(1:21),3][3:14])
tvals <- t(tvals)
rownames(tvals) <- rnames2
colnames(tvals) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
round(tvals, 3)
round(tvals*(abs(tvals)>1.96), 2)
round(tvals*(abs(tvals)>1.96), 2)[,c(7:12,1:6)]

# compare the two
round(tvals - tvals2, 2)
tmp <- sign(tvals - tvals2)
table(tmp) # overall
apply(tmp, 1, table)
apply(tmp, 2, table)



# Cluster etc -------------------------------------------------------------

library(dynamicTreeCut)
rdat <- as.matrix(std.rdats.all$r.pFFA)

mean.rdat <- rowMeans(rdat)

cmat <- cor(rdat)
dmat <- sqrt(2*(1-cmat)^2)
d <- as.dist(dmat)
hc <- hclust(d, method="ward.D2")
cl <- cutreeDynamic(hc, distM=dmat)
table(cl)
clmeans <- sapply(sort(unique(cl)), function(i) rowMeans(rdat[,cl==i]))

dmat <- Rfast::Dist(t(rdat))
d <- as.dist(dmat)
hc <- hclust(d, method="ward.D2")
cl2 <- cutreeDynamic(hc, distM=dmat)
table(cl2)
clmeans2 <- sapply(sort(unique(cl2)), function(i) rowMeans(rdat[,cl2==i]))

plot.ts(clmeans[1:200,1:4])
plot.ts(clmeans2[1:200,1:4])

get.tvals <- function(y, thr=T) {
  fits <- lm(y ~ ssubs + . + mcs, data=xmats[,-c(1,31:36)])
  sfits <- summary(fits)
  tvals <- sapply(sfits, function(sfit) sfit$coefficients[-c(1:21),3][3:14])
  tvals <- t(tvals)
  rownames(tvals) <- colnames(y)
  colnames(tvals) <- c(sprintf("shape%02i", 1:6), sprintf("resid%02i", 1:6))
  if (thr) {
    round(tvals*(abs(tvals)>1.96), 2)[,c(7:12,1:6)]
  } else {
    return(tvals[,c(7:12,1:6)])
  }
}

get.tvals(cbind(one=mean.rdat,two=mean.rdat))
get.tvals(clmeans)
get.tvals(clmeans2) # should choose this
mean(abs(get.tvals(clmeans)))
mean(abs(get.tvals(clmeans2)))


# now try the PCA (no this thing is confusing. the clusters though could be ok)
library(doMC)
registerDoMC(30)
comps <- getcomps.iter(rdat) # whoops taking too long...says 22
eigs <- eigen(cov(rdat), only.values=T)
nFactors::nScree(x=eigs$values, cor=F) # says 42 or 34
pca <- prcomp(rdat, retx=T)

get.tvals(pca$x[,1:22])

ica <- fastICA::fastICA(rdat, 22, method="C")
get.tvals(ica$S)

mean(abs(get.tvals(ica$S)))
mean(abs(get.tvals(pca$x[,1:34])))
