
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

# cor the prototypes
round(cor(cbind(avg.all, avg.male, avg.female)), 3)


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




# Prototype 1 -------------------------------------------------------------

# This uses the above avgs as the prototypes

# calculate distance to averages
dist.to.avg    <- Rfast::Dist(rbind(avg.all, vfeats))[-1,1]
dist.to.male   <- Rfast::Dist(rbind(avg.male, vfeats))[-1,1]
dist.to.female <- Rfast::Dist(rbind(avg.female, vfeats))[-1,1]
dist.btw.gender<- Rfast::Dist(rbind(avg.male, avg.female))[1,2]

# can there be some category measure? closeness to one category than another
dist.to.male/dist.to.female
(dist.to.male - dist.to.female)/dist.btw.gender
## want to ask if stuff close to each of the prototypes differs from the group further away
## so we can get the closest 25%

# run a regression to see how well these things fit the data
fit <- glm(vdf$gender ~ dist.to.avg + dist.to.male + dist.to.female, family=binomial(link='logit'))
summary(fit)



# Convolve regressors -----------------------------------------------------

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

# Get 
taskdir <- "/data1/famface01/analysis/task_activity"
regs <- cbind(avg=dist.to.avg, male=dist.to.male, female=dist.to.female, diff=abs((dist.to.male - dist.to.female)/dist.btw.gender))
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
  xmat2  <- afni.convolve.regressors(regs, vdf, subj, timing, max(timing$run))
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})

rregs <- regs
rregs[,2:3] <- lm(regs[,2:3] ~ regs[,1])$residuals
df.xmats2 <- ldply(subjects, function(subj) {
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
  xmat2  <- afni.convolve.regressors(rregs, vdf, subj, timing, max(timing$run))
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})


rregs2 <- cbind(gender=as.numeric(vdf$gender)-1.5, avg=rregs[,1], lm(rregs[,2:3] ~ vdf$gender)$residuals)
df.xmats3 <- ldply(subjects, function(subj) {
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
  xmat2  <- afni.convolve.regressors(rregs2, vdf, subj, timing, max(timing$run))
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})


# Run regression ----------------------------------------------------------

# Note: that gender effects present in just the average 
summary(lm(dist.to.avg ~ gender, data=vdf))

library(nlme)

lme.to.sdf <- function(sfit) {
  ts <- sfit$tTable[,4]
  zs <- qt(sfit$tTable[,5], Inf, lower.tail=F)
  ps <- sfits[[rname]]$tTable[,5]
  hemi <- sub("[.].*", "", rname)
  name <- sub("[lr]h.", "", rname)
  ord  <- rords[rnames==rname]
  sdf <- data.frame(hemi=hemi, roi=name, ord=ord, 
                    measure=rownames(sfits[[rname]]$tTable), 
                    tval=ts, pval=ps, zval=zs)
  sdf[-1,] # remove intercept
}

# get the lme results
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  motion <- as.matrix(df.xmats[,5:10])
  prots  <- as.matrix(df.xmats[,11:13])
  fit = lme(y ~ faces + quests + motion + prots,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames


ri <- 5
motion <- as.matrix(df.xmats[,5:10])
prots  <- lm(as.matrix(df.xmats[,12:13]) ~ df.xmats[,11])$residuals
fit1 = lme(y ~ faces + quests + motion + avg + prots,
           random = ~ 1|subject/run, 
           data=cbind(y=ts.mats[,ri], df.xmats), 
           control = lmeControl(opt = "optim"))
sfit1 <- summary(fit1)

avg <- df.xmats[,11]
prots  <- as.matrix(df.xmats[,12:13])
fit2 = lme(y ~ faces + quests + motion + avg + prots,
          random = ~ 1|subject/run, 
          data=cbind(y=ts.mats[,ri], df.xmats), 
          control = lmeControl(opt = "optim"))
sfit2 <- summary(fit2)

sfit2

tmp <- anova(sfit1)
tmp$`p-value`

anova(sfit2)
anova(sfit3)
anova(fit1, fit2)

motion <- as.matrix(df.xmats[,5:10])
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  #prots  <- as.matrix(df.xmats2[,12:13])
  fit = lme(y ~ faces + quests + motion + avg + diff,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

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

# add on the fdr values
library(fdrtool)
ret <- fdrtool(c(0.01/2,0.05/2,df.sfit$pval), statistic="pvalue", plot=F)
fdr.thr1 <- ret$qval[1]
fdr.thr2 <- ret$qval[2]
ret <- fdrtool(df.sfit$pval, statistic="pvalue", plot=F)
df.sfit$fdr.pval <- ret$qval

# get t-thr
dof <- sfits[[1]]$tTable[2,3]
fdr.tthr1 <- qt(fdr.thr1, dof, lower.tail=F)
fdr.tthr2 <- qt(fdr.thr2, dof, lower.tail=F)

head(df.sfit)


library(ggplot2)
library(ggthemr)
library(RColorBrewer)

ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
df.sfit2 <- subset(df.sfit, measure != "faces" & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



motion <- as.matrix(df.xmats3[,5:10])
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + motion + gender + avg + male,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats3), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames
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
df.sfit2 <- subset(df.sfit, measure %in% c("avg", "male", "gender") & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



# Random Prototypes --------------------------------------------------------

dmat <- Rfast::Dist(vfeats)

set.seed(42)
rand1 <- sample(vdf$ind[vdf$gender=="Male"], 1)
rand2 <- sample(vdf$ind[vdf$gender=="Female"], 1)

dist.to.rand1 <- dmat[rand1,-rand1]
dist.to.rand2 <- dmat[rand2,-rand2]

# Get 
taskdir <- "/data1/famface01/analysis/task_activity"
regs <- cbind(avg=dist.to.avg, male=dist.to.male, female=dist.to.female, rand1=dist.to.rand1, rand2=dist.to.rand2)
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
  xmat2  <- afni.convolve.regressors(regs, vdf, subj, timing, max(timing$run))
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})

motion <- as.matrix(df.xmats[,5:10])
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  #prots  <- as.matrix(df.xmats2[,12:13])
  fit = lme(y ~ faces + quests + motion + avg + rand2,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

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

ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
df.sfit2 <- subset(df.sfit, measure != "faces" & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))


# What if we don't do a random subset but instead we take the k-means
library(ClusterR)

# this says 2 (hmm so just showing that the neural network clearly clusters into male/female)
opt = Optimal_Clusters_KMeans(vfeats, max_clusters=20, plot_clusters = T, verbose = F, 
                              initializer = "optimal_init", criterion = 'distortion_fK')
opt = Optimal_Clusters_KMeans(vfeats, max_clusters=20, plot_clusters = T, verbose = F, 
                              initializer = "kmeans++", criterion = 'distortion_fK')

opt = Optimal_Clusters_KMeans(vfeats, max_clusters=20, plot_clusters = T, verbose = F, 
                              initializer = "optimal_init", criterion = 'silhouette')
opt = Optimal_Clusters_KMeans(vfeats, max_clusters=20, plot_clusters = T, verbose = F, 
                              initializer = "kmeans++", criterion = 'silhouette')

opt = Optimal_Clusters_KMeans(vfeats, max_clusters=20, plot_clusters = T, verbose = F, 
                              initializer = "optimal_init", criterion = 'variance_explained')

opt = Optimal_Clusters_KMeans(vfeats, max_clusters=30, plot_clusters = T, verbose = F, 
                              initializer = "optimal_init", criterion = 'Adjusted_Rsquared')


#criterion = "variance_explained", initializer = "optimal_init"
# 
#?Optimal_Clusters_KMeans

# Nearest Neighbor --------------------------------------------------------

dmat <- Rfast::Dist(vfeats)

k <- 30
gender.knn.probs <- t(sapply(1:nrow(dmat), function(i) {
  ninds <- order(dmat[i,])[-i][1:k]
  prop.table(table(vdf$gender[ninds]))
}))
gender.knn.dists <- sapply(1:nrow(dmat), function(i) {
  ninds <- order(dmat[i,])[-i][1:k]
  tab <- table(vdf$gender[ninds])
  sel <- names(which.max(tab))
  ds  <- Rfast::Dist(rbind(vfeats[i,], vfeats[ninds,][vdf$gender[ninds]==sel,]))[1,-1]
  mean(ds)
})
cor(cbind(gender.knn.probs, gender.knn.dists, dist.to.avg, dist.to.male, dist.to.female))



# Get 
taskdir <- "/data1/famface01/analysis/task_activity"


ps <- apply(gender.knn.probs, 1, min)
regs <- cbind(avg=dist.to.avg, avg.male=dist.to.male, avg.female=dist.to.female, 
              knn.prob=ps, knn.dist=gender.knn.dists)
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
  xmat2  <- afni.convolve.regressors(regs, vdf, subj, timing, max(timing$run), 
                                     ignore.stdout=T, ignore.stderr=T)
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})

motion <- as.matrix(df.xmats[,5:10])
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + motion + avg + avg.male + knn.prob + knn.dist,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

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

df.sfit2 <- subset(df.sfit, measure != "faces" & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))

## manual knn

library(caret)
library(glmnet)
registerDoMC(24)

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

extract_probs <- function(fit) {
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], vdf$gender)
  mean.preds <- ddply(preds, .(rowIndex), colwise(mean, .(Female,Male)))
  return(mean.preds)
}

X <- vfeats; y <- vdf$gender

library(kernlab)
fit <- run_caret(X, y, "svmLinear2", tlen=10)
fit2 <- run_caret(X, y, "svmLinear3", tlen=10)
mean.preds <- extract_probs(fit)
fit$finalModel$alpha


knnFit <- run_caret(X, y, "knn")
mean.preds <- extract_probs(knnFit)

# Get 
taskdir <- "/data1/famface01/analysis/task_activity"

ps <- apply(gender.knn.probs, 1, min)
ps2 <- apply(mean.preds[,-1], 1, min)
regs <- cbind(avg=dist.to.avg, avg.male=dist.to.male, avg.female=dist.to.female, 
              knn.prob=ps, knn.dist=gender.knn.dists, knn.pred=ps2)
round(cor(regs), 3)

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
  xmat2  <- afni.convolve.regressors(regs, vdf, subj, timing, max(timing$run), 
                                     ignore.stdout=T, ignore.stderr=T)
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})

motion <- as.matrix(df.xmats[,5:10])
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + motion + avg + avg.male + knn.prob + knn.dist,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

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

fdr.tthr1 <- 2.7
df.sfit2 <- subset(df.sfit, measure != "faces" & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))


# note: with higher k the class prob is higher!
# todo: maybe do some comparison here. a higher k, might have to do more with some category?





# Category ----------------------------------------------------------------

library(e1071)

# Get distance from hyperplane as measure of category
## from: http://mvpa.blogspot.com/2014/02/code-snippet-extracting-weights-from.html

X <- vfeats; y <- vdf$gender
fit <- svm(X, y, type="C-classification", kernel="sigmoid", cost=1, scale=FALSE, probability=T)
ret <- predict(fit, test.data, decision.values=T, probability = T)
head(attr(ret, "decision.values")) # same as manual calc below
#tmp1 <- attr(ret, "decision.values")
#tmp2 <- attr(ret, "decision.values")

w <- t(fit$coefs) %*% fit$SV;  
b <- -1 * fit$rho # (sometimes called w0)

test.data <- X
dist.hyp <- as.numeric((w %*% t(test.data)) + b) / sqrt(w %*% t(w))
labs <- sign((w %*% t(test.data)) + b)
table(vdf$gender, labs) # mostly gets it!

dist.hyp <- attr(ret, "decision.values")

cor(cbind(dist.hyp, gender.knn.dists, dist.to.male)) # ok related

regs <- cbind(avg=dist.to.avg, avg.male=dist.to.male, 
              dist.hyp=dist.hyp, dist.hyp2=dist.hyp^2, dist.inv.hyp=1/dist.hyp)
round(cor(regs), 3)

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
  xmat2  <- afni.convolve.regressors(regs, vdf, subj, timing, max(timing$run), 
                                     ignore.stdout=T, ignore.stderr=T)
  
  xmat    <- cbind(xmat1, xmat2)
  
  cat("\n")
  
  return(xmat)
})

motion <- as.matrix(df.xmats[,5:10])
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + motion + avg + avg.male + dist.hyp,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

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

fdr.tthr1 <- 2.7
df.sfit2 <- subset(df.sfit, measure != "faces" & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(aes(fill=roi), color="black", size=3, shape = 21) + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))





# Viz ---------------------------------------------------------------------

library(MASS)

# MDS
dmat <- Rfast::Dist(vfeats)
mds <- isoMDS(dmat, k=2)

plot(mds$points, col=as.numeric(vdf$gender)+1)

d <- as.dist(dmat)
tsne1 <- tsne::tsne(d)
tsne2 <- Rtsne::Rtsne(vfeats)

plot(tsne1, col=as.numeric(vdf$gender)+1)
plot(tsne2$Y, col=as.numeric(vdf$gender)+1)

pca <- prcomp(vfeats, retx=T)
plot(pca$x[,1], pca$x[,2], col=as.numeric(vdf$gender)+1)

library(e1071)
gender <- vdf$gender
dat <- tsne1; colnames(dat) <- c("x", "y"); dat <- as.data.frame(dat)
fit = svm(gender~., data=dat, cost=0.25, type='C-classification', kernel='linear')
fit
plot(fit, dat, x~y)
?plot.svm

