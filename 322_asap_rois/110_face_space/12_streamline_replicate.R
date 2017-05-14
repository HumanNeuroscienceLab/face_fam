
# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(doMC)
registerDoMC(30)
library(plyr)
library(bigmemory)
library(nlme)

library(ggplot2)
library(ggthemr)
library(RColorBrewer)

ggthemr('pale', type='outer', text_size=14, layout='plain')

subjects <- sprintf("sub%02i", 1:6)



# Load --------------------------------------------------------------------

setwd("/data1/famface01/command/misc/face_representations/322_asap_rois/")

# regressors etc
fdf <- read.csv("120_gender/z_df_info.csv")
fdf <- fdf[,-1]
head(fdf)

# face features
bdir   <- "/data1/famface01/analysis/misc/openface/concatenate_layers"
dfile  <- sprintf("clayer%02i.desc", 26)
bm     <- attach.big.matrix(file.path(bdir, dfile))
feats  <- as.matrix(bm)

# reduce to only vids
vdf    <- fdf[fdf$fr==3,]
vdf$ind<- 1:nrow(vdf)
vfeats <- apply(feats, 2, function(x) tapply(x, fdf$vid, mean))

# read in average
base <- "/data1/famface01/analysis/misc/openface"
labels <- read.csv(sprintf('%s/masked_labels_redo.csv', base), header=F)
features <- read.csv(sprintf('%s/masked_reps_redo.csv', base), header=F)
## find the average face
aind <- grep("average_face", labels$V2)
avg.all <- as.numeric(features[aind,])

# read in gender average
base <- "/data1/famface01/analysis/misc/openface/prototypes_gender"
labels <- read.csv(sprintf('%s/labels.csv', base), header=F)
features <- read.csv(sprintf('%s/reps.csv', base), header=F)
avg.male <- as.numeric(features[1,])
avg.female <- as.numeric(features[2,])

# remove the effect of the average face

# cor the prototypes
round(cor(cbind(avg.all, avg.male, avg.female)), 3)

# read in pca stuff
basedir      <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
## shapes
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)
oinds <- sapply(vdf$vid, function(x) which(shape.vnames == x))
all.equal(as.character(vdf$vid), shape.vnames[oinds])
sym.shapes <- sym.shapes[oinds,]
sym.shapes <- sym.shapes[,1:50]
## texture
sym.textures <- read.csv(file.path(basedir, "texture_sym_pca_scores.csv"), row.names = 1)
all.equal(rownames(sym.textures), shape.vnames)
sym.textures <- sym.textures[oinds,1:200]
## tot
pca.feats <- cbind(sym.shapes, sym.textures)


# Load the time-series ----------------------------------------------------

predir   <- "/data1/famface01/analysis/preprocessed"

# get the roi names
peak.tab <- read.csv("/data1/famface01/command/misc/face_representations/240_roi/z_asap_allpeaks.csv")
rnames <- paste(peak.tab$hemi, peak.tab$name, sep=".")
rords  <- rep(1:10, 2)

# load the rois
ts.mats <- ldply(subjects, function(subj) {
  tsdir   <- sprintf("%s/%s/func/unfam_vids/rois_asap_ventral_peaks", predir, subj)
  tsfile  <- file.path(tsdir, "asap_ventral_peaks_all.1D")
  ts <- read.table(tsfile)
  colnames(ts) <- rnames
  ts
}, .progress="text")
ts.mats <- as.matrix(ts.mats)


# and the timing
load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}
# Load timing information
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
head(dat.vols$sub01$basics$timing)



# Regression Functions ----------------------------------------------------

xmat_labs <- function(fn) {
  str <- system(sprintf("grep ColumnLabels %s | sed s/'#  ColumnLabels = '//", fn), intern=T)
  str <- gsub("\"", "", str)
  cols <- strsplit(str, ' ; ')[[1]]
  cols
}
read.xmat <- function(fn, rm.nums=T) {
  xmat <- read.table(fn)
  cnames <- xmat_labs(fn)
  if (rm.nums) {
    cnames <- sub("#0$", "", cnames)
  }
  colnames(xmat) <- cnames
  xmat
}

afni.timing.amp <- function(onsets, amps, runs, nruns, center=F, scale=F) {
  amps <- scale(amps, center=center, scale=scale)
  lines <- sapply(1:nruns, function(i) {
    inds <- runs == i
    if (any(inds)) {
      amp  <- amps[inds]
      onset<- onsets[inds]
      line <- paste(sprintf("%.5f*%.8f", onset, amp), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  
  as.character(lines)
}

afni.convolve.regressors <- function(regressors, vdf, subj, timing, nruns=max(timing$run), ...) {
  regressors <- as.matrix(regressors)
  rlabs <- colnames(regressors)
  if (is.null(rlabs)) rlabs <- sprintf("Stim%02i", 1:ncol(regressors))
  
  # Reorder regressors
  cur.vids <- as.character(vdf$vid)
  ref.vids <- as.character(timing$video)
  oinds   <- sapply(ref.vids, function(vname) which(vname==cur.vids))
  if (!all.equal(cur.vids[oinds], ref.vids)) stop("ordering issue")
  ro.regressors <- regressors[oinds,,drop=F]
  
  # Get the timing information
  lst.amps <- lapply(1:ncol(ro.regressors), function(i) {
    onsets <- timing$onset
    runs   <- timing$run
    amps   <- afni.timing.amp(onsets, ro.regressors[,i], runs, nruns, center=T)
    amps
  })
  
  # Save timing information to temp files
  amp.fnames <- sapply(1:length(lst.amps), function(i) {
    fname <- tempfile(pattern=sprintf("afni_timing_am%02i_", i), fileext=".1D")
    writeLines(lst.amps[[i]], fname)
    fname
  })
  
  # Run 3dDeconvolve
  infiles <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/unfam_vids/rois_asap_ventral_peaks/asap_ventral_peaks_run*.nii.gz", subj)
  ofile <- tempfile(pattern="xmat_", fileext=".1D")
  cmd0 <- "3dDeconvolve -global_times -input '%s' -force_TR 1 -polort -1 -x1D_stop -num_stimts %i %s -x1D %s"
  stim_cmds <- sapply(1:length(amp.fnames), function(i) {
    fname <- amp.fnames[i]
    lab   <- rlabs[i]
    sprintf("-stim_times_AM1 %i '%s' 'SPMG1(2)' -stim_label %i %s", i, fname, i, lab)
  })
  stim_cmds <- paste(stim_cmds, collapse=" ")
  cmd <- sprintf(cmd0, infiles, ncol(ro.regressors), stim_cmds, ofile)
  cat(cmd, "\n")
  retcode <- system(cmd, ...)
  if (retcode != 0) stop("ERROR: did not run properly\n")
  
  # Read in xmat
  xmat <- read.xmat(ofile)
  
  return(xmat)
}

gen.xmat <- function(design_mat) {
  # Get 
  taskdir <- "/data1/famface01/analysis/task_activity"
  df.xmats <- ldply(subjects, function(subj) {
    cat(subj, "\n")
    ## get runs/faces/quests/motion regressors
    xfile   <- sprintf("%s/%s/face_deviations_unfam/nnet2_only_avgdist2.reml/xmat.1D", taskdir, subj)
    xmat    <- read.xmat(xfile, rm.nums=T)
    xmat    <- xmat[,colnames(xmat)!="avg_dist"]
    # change the runs to one column indicating the run
    rinds <- grep("^Run.*Pol", colnames(xmat))
    rnums <- rowSums(sweep(xmat[,rinds], 2, 1:length(rinds), FUN="*"))
    runs  <- sprintf("run%02i", rnums)
    # add the subject and runs
    xmat1 <- cbind(subject=subj, run=runs, xmat[,-rinds])
    
    ## get regressors of interest
    timing <- dat.vols[[subj]]$basics$timing
    xmat2  <- afni.convolve.regressors(design_mat, vdf, subj, timing, max(timing$run), 
                                       ignore.stdout=T, ignore.stderr=T)
    
    xmat    <- cbind(xmat1, xmat2)
    
    cat("\n")
    
    return(xmat)
  })
  df.xmats
}

get.fits <- function(f, df.xmats) {
  sfits <- llply(1:ncol(ts.mats), function(ri) {
    cat(ri, "\n")
    
    fit = lme(f, random = ~ 1|subject/run, 
              data=cbind(y=ts.mats[,ri], df.xmats), 
              control = lmeControl(opt = "optim"))
    sfit <- summary(fit)
    sfit
  }, .parallel=T)
  names(sfits) <- rnames
  
  sfits
}

get.sfits <- function(sfits) {
  # put them in table form
  df.sfit <- ldply(rnames, function(rname) {
    # use all the p-values for the fdr (soon)
    tTable <- sfits[[rname]]$tTable
    
    # subset for reporting
    tTable <- tTable[-c(1,3,4:9),] # rm covariates
    ts <- tTable[,4]
    zs <- qt(tTable[,5], Inf, lower.tail=F)
    ps <- tTable[,5]
    hemi <- sub("[.].*", "", rname)
    name <- sub("[lr]h.", "", rname)
    ord  <- rords[rnames==rname]
    sdf <- data.frame(hemi=hemi, roi=name, ord=ord, 
                      measure=rownames(tTable), 
                      tval=ts, pval=ps, zval=zs)
    sdf
  })
  df.sfit
}

# classify
run_caret <- function(X, y, mthd, tlen=20) {
  colnames(X) <- sprintf("feat%02i", 1:ncol(X))
  
  nrepeats <- 10
  nfolds   <- 8
  fitControl <- trainControl(
    method = "repeatedcv", 
    number = nfolds, 
    repeats = nrepeats, 
    returnResamp = "final", 
    savePredictions = "final", 
    classProbs = T, 
    allowParallel = T
  )
  
  fit <- train(X, y, method = mthd, 
               trControl = fitControl, 
               preProcess = c("center","scale"), 
               tuneLength = tlen)
  ri <- as.numeric(rownames(fit$bestTune))
  print(fit$results[ri,])
  
  return(fit)  
}



# Prototype --------------------------------------------------------------------

# calculate distance to average prototypes
dist.to.avg    <- Rfast::Dist(rbind(avg.all, vfeats))[-1,1]
dist.to.avg2    <- Rfast::Dist(rbind(face.norm, vfeats))[-1,1]

# calculate framewise pose
pmat <- vdf[,c("pitch","yaw","roll")]
pmat <- rbind(c(0,0,0), pmat)
fd.pose <- sapply(1:(nrow(pmat)-1), function(ri) {
  #Rfast::dista(pmat[ri,], pmat[ri+1,])
  sqrt(sum((pmat[ri,] - pmat[ri+1,])^2))
})

# get design
design_mat <- cbind(avg=dist.to.avg, pose=vdf$fd, fd.pose=fd.pose, 
                    luminance=vdf$luminance)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
sfits1 <- get.fits(y ~ faces + quests + motion + avg, df.xmats)
df.sfits1 <- get.sfits(sfits1)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

sdf <- subset(df.sfits1, measure != "faces" & hemi == "lh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



# kNN ---------------------------------------------------------------------

dmat <- Rfast::Dist(vfeats)

# simple knn measure
k <- 10
kdists10 <- sapply(1:nrow(dmat), function(i) {
  mean(sort(as.numeric(dmat[i,-i]))[1:k])
})
ks <- c(seq(10,50,10), 75, 100, 150, 200)
lst.kdists <- lapply(ks, function(k) {
  sapply(1:nrow(dmat), function(i) {
    mean(sort(as.numeric(dmat[i,-i]))[1:k])
  })
})
names(lst.kdists) <- sprintf("k%03i", ks)

# get design
design_mat <- sapply(lst.kdists, function(x) x)
colnames(design_mat) <- sprintf("kdists%03i", ks)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
sfits1 <- get.fits(y ~ faces + quests + motion + kdists010, df.xmats)
df.sfits1 <- get.sfits(sfits1)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

sdf <- subset(df.sfits1, measure != "faces" & hemi == "lh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



# Compare ------------------------------------------------------------------

# get design
design_mat <- cbind(avg.dist=dist.to.avg, knn.dist=lst.kdists$k010)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
sfits1 <- get.fits(y ~ faces + quests + motion + avg.dist + knn.dist, df.xmats)
df.sfits1 <- get.sfits(sfits1)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

sdf <- subset(df.sfits1, measure != "faces" & hemi == "lh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))


## save the vals
subset(df.sfits1, measure!="faces" & roi%in%c("OFA","mFFA-1") & hemi=="rh")





# Prototype (PCA) -----------------------------------------------------------

# calculate distance to average
dist.to.avg    <- Rfast::Dist(rbind(rep(0,ncol(pca.feats)), pca.feats))[-1,1]

# get design
design_mat <- cbind(avg=dist.to.avg, luminance=vdf$luminance)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
sfits1 <- get.fits(y ~ faces + quests + motion + avg, df.xmats)
df.sfits1 <- get.sfits(sfits1)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

sdf <- subset(df.sfits1, measure != "faces" & hemi == "lh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



# kNN (PCA) -----------------------------------------------------------------

dmat <- Rfast::Dist(pca.feats)

# simple knn measure
k <- 10
kdists10 <- sapply(1:nrow(dmat), function(i) {
  mean(sort(as.numeric(dmat[i,-i]))[1:k])
})
ks <- c(seq(10,50,10), 75, 100, 150, 200)
lst.kdists <- lapply(ks, function(k) {
  sapply(1:nrow(dmat), function(i) {
    mean(sort(as.numeric(dmat[i,-i]))[1:k])
  })
})
names(lst.kdists) <- sprintf("k%03i", ks)

# get design
design_mat <- sapply(lst.kdists, function(x) x)
colnames(design_mat) <- sprintf("kdists%03i", ks)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
sfits1 <- get.fits(y ~ faces + quests + motion + kdists050, df.xmats)
df.sfits1 <- get.sfits(sfits1)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

sdf <- subset(df.sfits1, measure != "faces" & hemi == "lh")
ggplot(sdf, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(sdf$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

