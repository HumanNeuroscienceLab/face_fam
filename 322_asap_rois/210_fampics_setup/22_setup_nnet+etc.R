
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
avg.all2 <- as.numeric(features[aind,])
avg.all  <- colMeans(vfeats)


# Load nnet + behav -------------------------------------------------------

indir <- "/data1/famface01/command/misc/face_representations/322_asap_rois/210_fampics_setup"

## nnets
face.norm <- faces.avg <- faces.exp <- faces.fullset <- faces.exp.targ <- faces.exp2 <- faces.avg.targ <- NULL
load(file.path(indir, "z_face_nnet_feats.rda"), verbose=T)

## behav
behav.df <- read.csv(file.path(indir, "z_likeness+rt+gender+set.csv"))[,-1]

## Save the face exp
vdf2    <- faces.exp2$df

## target faces
trgts <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")

# calculate the prototypes
fullfeats <- faces.fullset$feats[faces.fullset$df$person %in% trgts,]
fullppls  <- factor(faces.fullset$df$person[faces.fullset$df$person %in% trgts])
avg.clusts <- apply(fullfeats, 2, function(x) tapply(x, fullppls, mean))

# calculate distance to exemplars
k <- 10
dmat <- Rfast::Dist(rbind(faces.exp2$feats, vfeats))[1:nrow(faces.exp2$feats),-c(1:nrow(faces.exp2$feats))]
kdist <- sapply(1:nrow(dmat), function(i) {
  mean(sort(dmat[i,])[1:k])
})


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


## SAVE THE WORKSPACE
flist <- c("afni.convolve.regressors", "afni.timing.amp", "avg.all", "avg.all2", "behav.df", ls()[grep("face", ls())], "fdf", "feats", "gen.xmat", "get.fits", "get.sfits", "indir", "base", "labels", "load.cur", "peak.tab", "predir", "read.xmat", "rnames", "rords", "run_caret", "subjects", "ts.mats", "vdf", "vdf2", "vfeats", "xmat_labs", "trgts", "fullfeats", "fullppls", "avg.clusts", "kdist")
save(list=flist, file="/data1/famface01/command/misc/face_representations/322_asap_rois/210_fampics_setup/z_fampics_vars+funs.Rdata")
#save.image(file="/data1/famface01/command/misc/face_representations/322_asap_rois/210_fampics_setup/z_fampics_vars+funs.Rdata")
