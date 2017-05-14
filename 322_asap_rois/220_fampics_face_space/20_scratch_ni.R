
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

library(xlsx)

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



# Load nnet + behav -------------------------------------------------------

indir <- "/data1/famface01/command/misc/face_representations/322_asap_rois/210_fampics_setup"

## nnets
face.norm <- faces.avg <- faces.exp <- faces.fullset <- faces.exp.targ <- faces.exp2 <- faces.avg.targ <- NULL
load(file.path(indir, "z_face_nnet_feats.rda"), verbose=T)

## behav
behav.df <- read.csv(file.path(indir, "z_likeness+rt+gender+set.csv"))[,-1]

vdf2 <- faces.exp2$df



# Load the time-series / timing ----------------------------------------------

predir   <- "/data1/famface01/analysis/preprocessed"

# get the roi names
peak.tab <- read.csv("/data1/famface01/command/misc/face_representations/240_roi/z_asap_allpeaks.csv")
rnames <- paste(peak.tab$hemi, peak.tab$name, sep=".")
rords  <- rep(1:10, 2)

# load the rois
ts.mats <- ldply(subjects, function(subj) {
  tsdir   <- sprintf("%s/%s/func/fam_pics/rois_asap_ventral_peaks", predir, subj)
  tsfile  <- file.path(tsdir, "asap_ventral_peaks_all.1D")
  ts <- read.table(tsfile)
  colnames(ts) <- rnames
  ts
}, .progress="text")
ts.mats <- as.matrix(ts.mats)

# and the timing
load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/SubjRois_FamPics/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  ret <- NULL
  load(infile)
  ret
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

afni.convolve.regressors <- function(regressors, subj, timing, vdf=NULL, nruns=max(timing$run), ...) {
  regressors <- as.matrix(regressors)
  rlabs <- colnames(regressors)
  if (is.null(rlabs)) rlabs <- sprintf("Stim%02i", 1:ncol(regressors))
  
  # Reorder regressors
  if (is.null(vdf)) {
    ro.regressors <- regressors
  } else {
    cur.vids <- as.character(vdf$pic)
    ref.vids <- as.character(timing$pic)
    oinds   <- sapply(ref.vids, function(vname) which(vname==cur.vids))
    if (!all.equal(cur.vids[oinds], ref.vids)) stop("ordering issue")
    ro.regressors <- regressors[oinds,,drop=F]
  }
  
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
  infiles <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/fam_pics/rois_asap_ventral_peaks/asap_ventral_peaks_run*.nii.gz", subj)
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

gen.xmat <- function(design_mat, include.identity=T, vdf0=NULL) {
  # Get 
  taskdir <- "/data1/famface01/analysis/task_activity"
  df.xmats <- ldply(subjects, function(subj) {
    cat(subj, "\n")
    ## get runs/faces/quests/motion regressors
    xfile   <- sprintf("%s/%s/fampics/fampics_main+nnet_effect.reml/xmat.1D", taskdir, subj)
    xmat    <- read.xmat(xfile, rm.nums=T)
    xmat    <- xmat[,!(colnames(xmat)%in%c("avg_diff","rt_dur"))]
    # change the runs to one column indicating the run
    rinds <- grep("^Run.*Pol", colnames(xmat))
    rnums <- rowSums(sweep(xmat[,rinds], 2, 1:length(rinds), FUN="*"))
    runs  <- sprintf("run%02i", rnums)
    # add the subject and runs
    xmat1 <- cbind(subject=subj, run=runs, xmat[,-rinds])
    
    ## get timing
    timing <- dat.vols[[subj]]$basics$timing
    timing$pic <- sub(".jpg", "", basename(as.character(timing$pic)))
    
    ## add the identity
    male   <- c("Justin_Timberlake", "Will_Smith", "Brad_Pitt", "Johnny_Depp", 
                "Barack_Obama", "Leonardo_DiCaprio", "Michael_Jackson", 
                "Tom_Cruise", "Will_Smith", "Bill_Clinton", "Dwayne_Johnson", 
                "George_Clooney", "Johnny_Depp", "Morgan_Freeman")
    timing$gender <- factor((timing$name %in% male)*1, levels=c(0,1), 
                            labels=c("female", "male"))
    
    timing$name <- factor(timing$name)
    timing$run  <- factor(timing$run)
    
    targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
    sinds <- (timing$name %in% targets)
    timing$name2 <- as.character(timing$name)
    timing$name2[!sinds] <- "Distractors"
    timing$name2 <- factor(timing$name2)
    
    tmp <- model.matrix(~ run + name2 + gender - 1, data=timing)
    id.mat <- tmp[,-grep("run", colnames(tmp))]
    
    timing$run  <- as.numeric(timing$run)
    xmat2 <- afni.convolve.regressors(id.mat, subj, timing, nruns=max(timing$run), ignore.stdout=T, ignore.stderr=T)
  
    ## get regressors of interest
    if (is.null(vdf0)) {
      vdf0 <- vdf2
      
      targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
      timing <- subset(timing, name %in% targets)
    }
    xmat3  <- afni.convolve.regressors(design_mat, subj, timing, vdf0, max(timing$run), 
                                       ignore.stdout=T, ignore.stderr=T)
    
    xmat    <- cbind(xmat1, xmat2, xmat3)
    
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
run_caret <- function(X, y, mthd, tlen=20, ...) {
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
               tuneLength = tlen, ...)
  ri <- as.numeric(rownames(fit$bestTune))
  print(fit$results[ri,])
  
  return(fit)  
}



# Prototype ------------------------------------------------------------------

# load in the other set
labs <- read.csv("/data1/famface01/data/stimuli/pics/landmark_facedb/nnet/labels.csv", header=F)
reps <- read.csv("/data1/famface01/data/stimuli/pics/landmark_facedb/nnet/reps.csv", header=F)
## select images only (not vid frames)
fnames <- basename(as.character(labs$V2))
inds <- grep("^img_", fnames)
reps <- reps[inds,]

faces.exp2$df$person <- factor(faces.exp2$df$person)
trgts <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
             "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")

# calculate distance to exemplars
k <- 10
dmat <- Rfast::Dist(rbind(faces.exp2$feats, vfeats))[1:nrow(faces.exp2$feats),-c(1:nrow(faces.exp2$feats))]
kdist0 <- sapply(1:nrow(dmat), function(i) {
  mean(sort(dmat[i,])[1:k])
})
#dmat <- Rfast::Dist(rbind(faces.exp2$feats, reps))[1:nrow(faces.exp2$feats),-c(1:nrow(faces.exp2$feats))]
#kdist0 <- sapply(1:nrow(dmat), function(i) {
#  mean(sort(dmat[i,])[1:k])
#})

# calculate distance to average prototypes
#dist.to.avg0    <- Rfast::Dist(rbind(rep(0,128), faces.exp2$feats))[-1,1]
dist.to.avg0    <- Rfast::Dist(rbind(colMeans(vfeats), faces.exp2$feats))[-1,1]
#dist.to.avg0    <- Rfast::Dist(rbind(colMeans(reps), faces.exp2$feats))[-1,1]
dist.to.avg1    <- Rfast::Dist(rbind(avg.all, faces.exp2$feats))[-1,1]
dist.to.avg2    <- Rfast::Dist(rbind(face.norm, faces.exp2$feats))[-1,1]

# can do distance to each cluster center
avg.clusts <- apply(faces.exp2$feats, 2, function(x) tapply(x, faces.exp2$df$person, mean))
cor(t(avg.clusts), t(faces.avg.targ$feats))
## other independent version of clust average
fullfeats <- faces.fullset$feats[faces.fullset$df$person %in% trgts,]
fullppls  <- factor(faces.fullset$df$person[faces.fullset$df$person %in% trgts])
avg.clusts <- apply(fullfeats, 2, function(x) tapply(x, fullppls, mean))
## now calc
dist.to.prots   <- Rfast::Dist(rbind(avg.clusts, faces.exp2$feats))[-(1:nrow(avg.clusts)),1:nrow(avg.clusts)]
colnames(dist.to.prots) <- rownames(avg.clusts)
dist.to.prots2  <- sapply(rownames(avg.clusts), function(pname) {
  inds <- faces.exp2$df$person == pname
  x <- dist.to.prots[,pname == colnames(dist.to.prots)]
  x[!inds] <- 0
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
dist.to.prots3 <- rowSums(dist.to.prots2)


# get design
design_mat <- cbind(kdist0=kdist0, avg0=dist.to.avg0, avg1=dist.to.avg1, avg2=dist.to.avg2, davg3=dist.to.prots3, dist.to.prots2)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')
fdr.tthr1 <- 2.7
#sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")

sdf <- subset(df.sfits1, measure%in%c("avg0","davg3"))
sdf$measure <- factor(sdf$measure)
levels(sdf$measure) <- c("Average Face", "Person Face")
ggplot(sdf, aes(x=ord, y=tval)) + 
  geom_hline(yintercept=0, color="grey25") + 
  geom_line(aes(color=hemi)) + 
  scale_color_hue() + 
  geom_point(color="black", size=3, shape = 21) + 
  geom_hline(yintercept=c(-1,1)*fdr.tthr1, linetype='dotted', size=0.5) + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  facet_grid(measure~.) + 
  expand_limits(y=range(sdf$tval)*1.05) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



# Exemplar to Prot --------------------------------------------------------

# This will run the exemplar in a way similar to the prototype

## exemplar distances (dist to everything and then average)
n <- nrow(faces.exp2$feats)
kdist.to.avg0    <- Rfast::Dist(rbind(faces.exp2$feats, vfeats))
kdist.to.avg0    <- kdist.to.avg0[1:n,-c(1:n)]
kdist.to.avg0    <- rowSums(kdist.to.avg0)

## exemplar person distance
fullfeats <- faces.fullset$feats[faces.fullset$df$person %in% trgts,]
fullppls  <- factor(faces.fullset$df$person[faces.fullset$df$person %in% trgts])
avg.clusts <- apply(fullfeats, 2, function(x) tapply(x, fullppls, mean))

dist.to.exemps <- Rfast::Dist(rbind(faces.exp2$feats, fullfeats))
dist.to.exemps <- dist.to.exemps[1:n,-c(1:n)]
dist.to.exemps <- sapply(1:n, function(i) {
  pname <- as.character(faces.exp2$df$person[i])
  mean(dist.to.exemps[i,fullppls == pname])
})
dist.to.exemps2  <- sapply(levels(faces.exp2$df$person), function(pname) {
  inds <- faces.exp2$df$person == pname
  x <- dist.to.exemps
  x[!inds] <- 0
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
dist.to.exemps3 <- rowSums(dist.to.exemps2)

## do an exemplar distance here
edist.to.avg0 <- exp(-dist.to.avg0)

# get design
design_mat <- cbind(kdist0=kdist0, edavg=edist.to.avg0, davg=dist.to.avg0, kavg=kdist.to.avg0, davg3=dist.to.prots3, kavg3=dist.to.exemps3, dist.to.prots2)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## try with prot dists first
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + davg + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)
## try with exemp dists second
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + kavg + kavg3, df.xmats)
df.sfits2 <- get.sfits(sfits2)
anova(sfits2$`rh.mFFA-1`)
anova(sfits2$`lh.mFFA-1`)
## try with exp (same!)
sfits3 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + edavg + davg3, df.xmats)
df.sfits3 <- get.sfits(sfits3)
anova(sfits3$`rh.mFFA-1`)
anova(sfits3$`lh.mFFA-1`)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')
fdr.tthr1 <- 2.7
#sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")

sdf <- subset(df.sfits1, measure%in%c("avg0","davg3"))
sdf$measure <- factor(sdf$measure)
levels(sdf$measure) <- c("Average Face", "Person Face")
ggplot(sdf, aes(x=ord, y=tval)) + 
  geom_hline(yintercept=0, color="grey25") + 
  geom_line(aes(color=hemi)) + 
  scale_color_hue() + 
  geom_point(color="black", size=3, shape = 21) + 
  geom_hline(yintercept=c(-1,1)*fdr.tthr1, linetype='dotted', size=0.5) + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  facet_grid(measure~.) + 
  expand_limits(y=range(sdf$tval)*1.05) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))



# Relative Prototype ------------------------------------------------------

# Here we see if the relative summed similarity work
# get the avg.clusts. compare to yourself and to your set.

targetsA <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith")
targetsB <- c("Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")

# Collect ref set
fullfeats <- faces.fullset$feats[faces.fullset$df$person %in% trgts,]
fullppls  <- factor(faces.fullset$df$person[faces.fullset$df$person %in% trgts])
avg.clusts <- apply(fullfeats, 2, function(x) tapply(x, fullppls, mean))
## now calc
dist.to.prots   <- Rfast::Dist(rbind(avg.clusts, faces.exp2$feats))[-(1:nrow(avg.clusts)),1:nrow(avg.clusts)]
colnames(dist.to.prots) <- rownames(avg.clusts)
rdist.to.prots2  <- sapply(rownames(avg.clusts), function(pname) {
  # Compare to prototype
  inds <- faces.exp2$df$person == pname
  x1 <- dist.to.prots[,pname == colnames(dist.to.prots)]
  # Compare to other four prots
  if (any(targetsA == pname)) {
    cinds <- colnames(dist.to.prots) %in% targetsA
  } else {
    cinds <- colnames(dist.to.prots) %in% targetsB
  }
  x2 <- dist.to.prots[,cinds]
  # Compare self to others
  x <- x1/rowMeans(x2)
  # Scale
  x[!inds] <- 0
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
rdist.to.prots3 <- rowSums(rdist.to.prots2)

# get design
design_mat <- cbind(kdist0=kdist0, davg=dist.to.avg0, davg3=dist.to.prots3, ravg3=rdist.to.prots3, ravg3=rdist.to.prots3, dist.to.prots2)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## try with prot dists first
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + davg + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)
## try with relative dists second
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + davg + ravg3, df.xmats)
df.sfits2 <- get.sfits(sfits2)
anova(sfits2$`rh.mFFA-1`)
anova(sfits2$`lh.mFFA-1`)
## try with exp (same!)
sfits3 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + edavg + davg3, df.xmats)
df.sfits3 <- get.sfits(sfits3)
anova(sfits3$`rh.mFFA-1`)
anova(sfits3$`lh.mFFA-1`)

# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')
fdr.tthr1 <- 2.7
#sdf <- subset(df.sfits1, measure != "faces" & hemi == "rh")

sdf <- subset(df.sfits1, measure%in%c("avg0","davg3"))
sdf$measure <- factor(sdf$measure)
levels(sdf$measure) <- c("Average Face", "Person Face")
ggplot(sdf, aes(x=ord, y=tval)) + 
  geom_hline(yintercept=0, color="grey25") + 
  geom_line(aes(color=hemi)) + 
  scale_color_hue() + 
  geom_point(color="black", size=3, shape = 21) + 
  geom_hline(yintercept=c(-1,1)*fdr.tthr1, linetype='dotted', size=0.5) + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=sdf$ord, labels=sdf$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  facet_grid(measure~.) + 
  expand_limits(y=range(sdf$tval)*1.05) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))


# Classification ----------------------------------------------------------

# This is to weight the prototypes by the important features

# calculate the class probabilities
X.train <- faces.fullset$feats; y.train <- factor(faces.fullset$df$person)
X.train <- X.train[y.train %in% trgts,]; y.train <- factor(y.train[y.train %in% trgts])
X.test <- faces.exp2$feats; y.test <- factor(faces.exp2$df$person)

library(caret)
library(glmnet)
library(doMC)
registerDoMC(30)

tuneGrid <- ldply(c(0,0.5,1), function(alpha) {
  lambdas <- glmnet(as.matrix(X.train), y.train, family="multinomial", nlambda=50, alpha=0.5)$lambda
  data.frame(alpha=alpha, lambda=lambdas)
}, .parallel=T)
glmFit <- run_caret(X.train, y.train, "glmnet", tuneGrid=tuneGrid, family="multinomial")
ret <- predict(glmFit$finalModel, newx=as.matrix(X.test), s=glmFit$bestTune$lambda, type="response")
prob.to.prots <- ret[,,1]
prob.to.prots2  <- sapply(colnames(prob.to.prots), function(pname) {
  inds <- faces.exp2$df$person == pname
  x <- prob.to.prots[,pname == colnames(prob.to.prots)]
  x[!inds] <- 0
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
prob.to.prots3 <- rowSums(prob.to.prots2)
cor(prob.to.prots3, dist.to.prots3) # negatively related


###
# CLASSIFY: weighted distance
###

# okay but now get the coefficients
vimps <- varImp(glmFit)$importance
cfs <- predict(glmFit$finalModel, s=glmFit$bestTune$lambda, type="coefficients")
wts <- sapply(cfs, function(x) x[-1])
## can we use these to calculate the distances?
head(wts)

wdist.to.prots <- laply(1:nrow(avg.clusts), function(ai) {
  avg.vec <- avg.clusts[ai,]
  wt.vec  <- wts[,ai]
  vimp.vec<- vimps[,ai]
  laply(1:nrow(faces.exp2$feats), function(fi) {
    vec     <- faces.exp2$feats[fi,]
    sqrt(sum(  (wt.vec*(vec - avg.vec))^2  ))
  })
}, .parallel=T)
wdist.to.prots <- t(wdist.to.prots)
colnames(wdist.to.prots) <- rownames(avg.clusts)
wdist.to.prots2  <- sapply(rownames(avg.clusts), function(pname) {
  inds <- faces.exp2$df$person == pname
  x <- wdist.to.prots[,pname == colnames(dist.to.prots)]
  x[!inds] <- 0
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
wdist.to.prots3 <- rowSums(wdist.to.prots2)

# varimp is very similar
#vdist.to.prots <- laply(1:nrow(avg.clusts), function(ai) {
#  avg.vec <- avg.clusts[ai,]
#  wt.vec  <- vimps[,ai]
#  vimp.vec<- vimps[,ai]
#  laply(1:nrow(faces.exp2$feats), function(fi) {
#    vec     <- faces.exp2$feats[fi,]
#    sqrt(sum(  (wt.vec*(vec - avg.vec))^2  ))
#  })
#}, .parallel=T)
#vdist.to.prots <- t(vdist.to.prots)
#colnames(vdist.to.prots) <- rownames(avg.clusts)
#vdist.to.prots2  <- sapply(rownames(avg.clusts), function(pname) {
#  inds <- faces.exp2$df$person == pname
#  x <- vdist.to.prots[,pname == colnames(dist.to.prots)]
#  x[!inds] <- 0
#  x[inds] <- scale(x[inds], center=T, scale=F)
#  x
#})
#vdist.to.prots3 <- rowSums(vdist.to.prots2)
#cor(wdist.to.prots3, vdist.to.prots3) # highly similar

# get design
design_mat <- cbind(kdist0=kdist0, avg0=dist.to.avg0, avg1=dist.to.avg1, avg2=dist.to.avg2, davg3=dist.to.prots3, wavg3=wdist.to.prots3, probs=prob.to.prots3, dist.to.prots2)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## old stuff
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)
## weighted
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + wavg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)
## probs
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + probs, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)

## plot weighted compared to not weighted
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + wavg3, df.xmats)
df.sfits2 <- get.sfits(sfits2)
df.subset <- rbind(
  #subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits1, measure == "davg3" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits2, measure == "wavg3" & roi %in% c("OFA", "mFFA-1"))
)
df.subset$measure <- factor(df.subset$measure)
df.subset$roi <- factor(df.subset$roi)
levels(df.subset$measure) <- c("Prototype", "Weighted Prototype")
levels(df.subset$hemi) <- c("Right Hemi", "Left Hemi")
levels(df.subset$roi) <- c("OFA", "FFA")

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.subset, aes(x=roi, y=tval, fill=roi, alpha=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.6,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("T-Statistic") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(.~hemi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))#, 
        #axis.text.x = element_text(angle = 45, hjust = 1))

## plot the probs compared to the weights
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + probs, df.xmats)
df.sfits1 <- get.sfits(sfits1)
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + wavg3, df.xmats)
df.sfits2 <- get.sfits(sfits2)
df.subset <- rbind(
  #subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits1, measure == "probs" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits2, measure == "wavg3" & roi %in% c("OFA", "mFFA-1"))
)
df.subset$measure <- factor(df.subset$measure)
df.subset$roi <- factor(df.subset$roi)
levels(df.subset$measure) <- rev(c("Weighted Prototype", "Probability"))
levels(df.subset$hemi) <- c("Right Hemi", "Left Hemi")
levels(df.subset$roi) <- c("OFA", "FFA")

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.subset, aes(x=roi, y=abs(tval), fill=roi, alpha=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.6,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("Abs(T-Stat)") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(.~hemi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))#, 
#axis.text.x = element_text(angle = 45, hjust = 1))

## same as above but in same model
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + probs + wavg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + wavg3 + probs, df.xmats)
df.sfits2 <- get.sfits(sfits2)
df.subset <- rbind(
  #subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits1, measure == "probs" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits2, measure == "wavg3" & roi %in% c("OFA", "mFFA-1"))
)
df.subset$measure <- factor(df.subset$measure)
df.subset$roi <- factor(df.subset$roi)
levels(df.subset$measure) <- rev(c("Weighted Prototype", "Probability"))
levels(df.subset$hemi) <- c("Right Hemi", "Left Hemi")
levels(df.subset$roi) <- c("OFA", "FFA")

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.subset, aes(x=roi, y=tval, fill=roi, alpha=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.6,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("T-Stat") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(.~hemi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))#, 


###
# CLASSIFY: weighted prototype
###

train.df <- faces.fullset$df[faces.fullset$df$person %in% trgts,]

tmp <- predict(glmFit, type="prob")
tmp <- ddply(glmFit$pred, .(Resample), function(x) {
  x[order(x$rowIndex),]
})
tmp2 <- daply(tmp, .(rowIndex), function(x) {
  colMeans(x[,-c(1:5,14)])
})
wavg.clusts <- sapply(colnames(tmp2), function(pname) {
  wts <- tmp2[train.df$person == pname, pname]
  wts <- wts/sum(wts)
  xmat <- X.train[train.df$person == pname,]
  ave.vec <- apply(xmat, 2, function(x) weighted.mean(x, wts))
  ave.vec
})
wavg.clusts <- t(wavg.clusts)

pdist.to.prots   <- Rfast::Dist(rbind(wavg.clusts, faces.exp2$feats))[-(1:nrow(wavg.clusts)),1:nrow(wavg.clusts)]
colnames(pdist.to.prots) <- rownames(wavg.clusts)
pdist.to.prots2  <- sapply(rownames(wavg.clusts), function(pname) {
  inds <- faces.exp2$df$person == pname
  x <- dist.to.prots[,pname == colnames(pdist.to.prots)]
  x[!inds] <- 0
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
pdist.to.prots3 <- rowSums(pdist.to.prots2)


# get design
design_mat <- cbind(kdist0=kdist0, avg0=dist.to.avg0, avg1=dist.to.avg1, avg2=dist.to.avg2, davg3=dist.to.prots3, wavg3=wdist.to.prots3, pavg3=pdist.to.prots3, dist.to.prots2)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## old stuff
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)
## weighted prototype
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + pavg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)



# kNN --------------------------------------------------------------------

library(caret)
registerDoMC(30)

extract_probs <- function(fit) {
  ps <- factor(faces.exp2$df$person)
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], ps)
  mean.preds <- ddply(preds, .(rowIndex), colwise(mean, ps))
  return(mean.preds)
}

extract_preds <- function(fit) {
  ps <- factor(faces.exp2$df$person)
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], ps)
  mode.preds <- ddply(preds, .(rowIndex), function(x) {
    c(pred=names(which.max(table(x$pred))))
  })
  mode.preds
}


# calculate the class probabilities
X.train <- faces.fullset$feats; y.train <- factor(faces.fullset$df$person)
X.train <- X.train[y.train %in% trgts,]; y.train <- factor(y.train[y.train %in% trgts])
X.test <- faces.exp2$feats; y.test <- factor(faces.exp2$df$person)

knnFit <- run_caret(X.train, y.train, "knn", tlen=20)
probs <- predict(knnFit$finalModel, X.test)
preds <- predict(knnFit$finalModel, X.test, type="class")

# not always right
table(comp=preds, ref=factor(faces.exp2$df$person))

# collapse probs
probs2 <- sapply(1:ncol(probs), function(i) {
  cname <- colnames(probs)[i]
  inds <- as.character(y.test) == cname
  x <- probs[,i] 
  x[!inds] <- 0
  x[inds] <- scale(x[inds], scale=F, center=T)
  x
})
probs3 <- rowSums(probs2)

# get design
design_mat <- cbind(kdist0=kdist0, avg0=dist.to.avg0, avg1=dist.to.avg1, 
                    davg3=dist.to.prots3, wavg3=wdist.to.prots3, 
                    knn.prob=probs3, probs)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

## individual models
## covars
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## fit
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + knn.prob, df.xmats)
df.sfits1 <- get.sfits(sfits1)
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + wavg3, df.xmats)
df.sfits2 <- get.sfits(sfits2)
df.subset <- rbind(
  #subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits1, measure == "knn.prob" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits2, measure == "wavg3" & roi %in% c("OFA", "mFFA-1"))
)
df.subset$measure <- factor(df.subset$measure)
df.subset$roi <- factor(df.subset$roi)
levels(df.subset$measure) <- rev(c("Weighted Prototype", "Exemplar Probability"))
levels(df.subset$hemi) <- c("Right Hemi", "Left Hemi")
levels(df.subset$roi) <- c("OFA", "FFA")

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.subset, aes(x=roi, y=abs(tval), fill=roi, alpha=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.6,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("Abs(T-Stat)") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(.~hemi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))#, 


## OTHER STUFF
## do for avg and avg
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + davg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + davg3 + avg0, df.xmats)
df.sfits2 <- get.sfits(sfits2)
df.subset <- rbind(
  #subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits2, measure == "davg3" & roi %in% c("OFA", "mFFA-1"))
)
df.subset$measure <- factor(df.subset$measure)
df.subset$roi <- factor(df.subset$roi)
levels(df.subset$measure) <- c("Average Face", "Person Face")
levels(df.subset$hemi) <- c("Right Hemi", "Left Hemi")
levels(df.subset$roi) <- c("OFA", "FFA")

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.subset, aes(x=roi, y=abs(tval), fill=roi, alpha=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.6,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("Abs(T-Stat)") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(.~hemi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))#, 



## combined model
## covars
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## fit
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + knn.prob + wavg3, df.xmats)
df.sfits1 <- get.sfits(sfits1)
sfits2 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + wavg3 + knn.prob, df.xmats)
df.sfits2 <- get.sfits(sfits2)
df.subset <- rbind(
  #subset(df.sfits1, measure == "avg0" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits1, measure == "knn.prob" & roi %in% c("OFA", "mFFA-1")), 
  subset(df.sfits2, measure == "wavg3" & roi %in% c("OFA", "mFFA-1"))
)
df.subset$measure <- factor(df.subset$measure)
df.subset$roi <- factor(df.subset$roi)
levels(df.subset$measure) <- rev(c("Weighted Prototype", "Exemplar Probability"))
levels(df.subset$hemi) <- c("Right Hemi", "Left Hemi")
levels(df.subset$roi) <- c("OFA", "FFA")

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.subset, aes(x=roi, y=abs(tval), fill=roi, alpha=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.6,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("Abs(T-Stat)") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(.~hemi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))#, 





# SVM --------------------------------------------------------------------

library(caret)
library(e1071)
registerDoMC(30)

# Data
X.train <- faces.fullset$feats; y.train <- factor(faces.fullset$df$person)
X.train <- X.train[y.train %in% trgts,]; y.train <- factor(y.train[y.train %in% trgts])
X.test <- faces.exp2$feats; y.test <- factor(faces.exp2$df$person)

# classify
svmFit <- run_caret(X.train, y.train, "svmLinear2", tlen=20)

# preds/probs (not used just to see)
svm.probs <- extract_probs(svmFit)
svm.preds <- extract_preds(svmFit)
table(svm.preds$pred, y.train)

# test
svmFit2 <- run_caret(X.train, (y.train=="Angelina_Jolie")*1, "svmLinear2", tlen=10)
pred2 <- predict(svmFit2$finalModel, X.test, 
                 decision.values=T, probability=T)

# test preds/probs (will save to run regression)
pred <- predict(svmFit$finalModel, X.test, 
                decision.values=T, probability=T)
head(attr(pred, "decision.values"))
head(attr(pred, "probabilities"))
svm.dvals <- attr(pred, "decision.values")
svm.probs <- attr(pred, "probabilities")

# combine dvals
# there are 7 dvals per identity
unames <- levels(faces.exp2$df$person)
svm.dvals.ave <- sapply(1:length(unames), function(i) {
  tmp <- svm.dvals[,grep(unames[i], colnames(svm.dvals))]
  rowMeans(tmp)
})
colnames(svm.dvals.ave) <- unames
svm.dvals.min <- sapply(1:length(unames), function(i) {
  tmp <- svm.dvals[,grep(unames[i], colnames(svm.dvals))]
  apply(tmp, 1, function(x) x[which.min(x)])
})
colnames(svm.dvals.min) <- unames
svm.dvals.in <- t(sapply(1:nrow(svm.dvals), function(i) {
  cname <- unames[faces.exp2$df$person[i]==unames]
  vec <- vector("numeric", 8)
  
  x <- svm.dvals[i,grep(cname, colnames(svm.dvals))]
  
  cur.names <- sub("/", "", sub(cname, "", names(x)))
  inds <- sapply(cur.names, function(cn) which(unames==cn))
  vec[inds] <- x
  names(vec) <- unames
  
  vec
}))
colnames(svm.dvals.in) <- unames

# collapse probs
svm.probs2 <- sapply(colnames(svm.probs), function(name) {
  x       <- svm.probs[,colnames(svm.probs)==name]
  inds    <- (faces.exp2$df$person==name)
  x       <- x * inds
  x[inds] <- scale(x[inds], center=T, scale=T)
  x
})
svm.probs3 <- rowSums(svm.probs2)

# collapse dvals
svm.dvals.ave2 <- sapply(colnames(svm.dvals.ave), function(name) {
  x       <- svm.dvals.ave[,colnames(svm.dvals.ave)==name]
  inds    <- (faces.exp2$df$person==name)
  x       <- x * inds
  x[inds] <- scale(x[inds], center=T, scale=T)
  x
})
svm.dvals.ave2 <- rowSums(svm.dvals.ave2)
svm.dvals.min2 <- sapply(colnames(svm.dvals.min), function(name) {
  x       <- svm.dvals.ave[,colnames(svm.dvals.min)==name]
  inds    <- (faces.exp2$df$person==name)
  x       <- x * inds
  x[inds] <- scale(x[inds], center=T, scale=T)
  x
})
svm.dvals.min2 <- rowSums(svm.dvals.min2)

# get design
design_mat <- cbind(avg0=dist.to.avg0, avg1=dist.to.avg1, 
                    knn.prob=probs3, 
                    svm.prob=svm.probs3, svm.dval1=svm.dvals.ave2, svm.dval2=svm.dvals.min2, 
                    svm.probs)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
## covars
motion <- as.matrix(df.xmats[,6:11])
distracts <- as.numeric(df.xmats$name2Distractors)
ids <- as.matrix(df.xmats[,grep("^name2", colnames(df.xmats))])
ids <- ids[,colnames(ids) != "name2Distractors"]
dist.ids <- as.matrix(df.xmats[,trgts])
## probs
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + svm.prob, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)
## svm dist
sfits1 <- get.fits(y ~ faces + incorrect + motion + distracts + ids + gendermale + avg0 + svm.dval2, df.xmats)
df.sfits1 <- get.sfits(sfits1)
anova(sfits1$`rh.mFFA-1`)
anova(sfits1$`lh.mFFA-1`)


# plot
ggthemr('pale', type='outer', text_size=14, layout='plain')

# removed the gender from probs
fdr.tthr1 <- 2.7
sdf <- subset(df.sfits5, measure != "faces" & hemi == "rh")
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

sdf <- subset(df.sfits5, measure != "faces" & hemi == "lh")
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

