
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

extract_probs <- function(fit) {
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], vdf$gender)
  mean.preds <- ddply(preds, .(rowIndex), colwise(mean, .(Female,Male)))
  return(mean.preds)
}

extract_preds <- function(fit) {
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], vdf$gender)
  mode.preds <- ddply(preds, .(rowIndex), function(x) {
    c(pred=names(which.max(table(x$pred))))
  })
  mode.preds
}



# Prototype --------------------------------------------------------------------

# calculate distance to average prototypes
dist.to.avg    <- Rfast::Dist(rbind(avg.all, vfeats))[-1,1]
dist.to.male   <- Rfast::Dist(rbind(avg.male, vfeats))[-1,1]
dist.to.female <- Rfast::Dist(rbind(avg.female, vfeats))[-1,1]
dist.to.prots  <- Rfast::Dist(rbind(avg.all, avg.male, avg.female))

# remove the gender effect from distances
dist.to.male2  <- lm(dist.to.male ~ gender, data=vdf)$residuals
dist.to.female2 <- lm(dist.to.female ~ gender, data=vdf)$residuals

# remove the gender effect from distances
dist.to.male3  <- lm(dist.to.male ~ gender, data=vdf)$residuals * (vdf$gender == "Male")
dist.to.female3 <- lm(dist.to.female ~ gender, data=vdf)$residuals * (vdf$gender == "Female")

# get design
design_mat <- cbind(avg=dist.to.avg, avg.male=dist.to.male2, avg.female=dist.to.female2,
                    gender.m=(vdf$gender=="Male")*1, gender.f=(vdf$gender=="Female")*1, 
                    avg.male2=dist.to.male3, avg.female2=dist.to.female3)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
prots  <- as.matrix(df.xmats[,c("avg.male", "avg.female")])
sfits1 <- get.fits(y ~ faces + quests + motion + prots + gender.m, df.xmats)
df.sfits1 <- get.sfits(sfits1)

# get the fits for removing gender
motion    <- as.matrix(df.xmats[,5:10])
prots     <- as.matrix(df.xmats[,c("avg.male2", "avg.female2")])
sfits2    <- get.fits(y ~ faces + quests + motion + prots + gender.m, df.xmats)
df.sfits2 <- get.sfits(sfits2)

# plot (regular)
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

# plot (remove gender effects)
ggthemr('pale', type='outer', text_size=14, layout='plain')

fdr.tthr1 <- 2.7
sdf <- subset(df.sfits2, measure != "faces" & hemi == "rh")
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



# kNN --------------------------------------------------------------------

library(caret)
registerDoMC(30)

# calculate the class probabilities
X <- vfeats; y <- vdf$gender
knnFit <- run_caret(X, y, "knn", tlen=20)
mean.preds <- extract_probs(knnFit)
mode.preds <- extract_preds(knnFit)

# not always right
table(comp=mode.preds$pred, ref=vdf$gender)

# split the probs into two regressors, one for male and one for female
probs.male   <- mean.preds$Male * (vdf$gender == "Male")
probs.female <- mean.preds$Female * (vdf$gender == "Female")

# get design
tmp1 <- probs.male; tmp1[probs.male!=0] <- scale(probs.male[probs.male!=0], scale=F)
tmp2 <- probs.female; tmp1[probs.female!=0] <- scale(probs.female[probs.female!=0], scale=F)
design_mat <- cbind(pmale=probs.male, pfemale=probs.female, 
                    male=(vdf$gender=="Male")*1, female=(vdf$gender=="Female")*1, 
                    pmale2=tmp1, pfemale2=tmp2)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
probs  <- as.matrix(df.xmats[,c("pmale", "pfemale")])
sfits1 <- get.fits(y ~ faces + quests + motion + probs, df.xmats)
df.sfits1 <- get.sfits(sfits1)
## and combine probs with main effect of gender
sfits3 <- get.fits(y ~ faces + quests + motion + probs + male, df.xmats)
df.sfits3 <- get.sfits(sfits3)
sfits4 <- get.fits(y ~ faces + quests + motion + probs + female, df.xmats)
df.sfits4 <- get.sfits(sfits4)
## and last one
probs2  <- as.matrix(df.xmats[,c("pmale2", "pfemale2")])
sfits5 <- get.fits(y ~ faces + quests + motion + male + probs2, df.xmats)
df.sfits5 <- get.sfits(sfits5)


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



# Category Boundary: SVM --------------------------------------------------

library(caret)
library(e1071)
registerDoMC(30)

# get optimal vals
X <- vfeats; y <- vdf$gender
svmFit <- run_caret(X, y, "svmLinear2", tlen=20)
mean.preds <- extract_probs(svmFit)
mode.preds <- extract_preds(svmFit)
dvs <- svmFit$finalModel$decision.values

# not always right
table(comp=mode.preds$pred, ref=vdf$gender)

# split the probs into two regressors, one for male and one for female
probs.male   <- lm(mean.preds$Male ~ vdf$gender)$residuals * (vdf$gender == "Male")
probs.female <- lm(mean.preds$Female ~ vdf$gender)$residuals * (vdf$gender == "Female")

# split the decision values too
dvs.male   <- lm(dvs ~ vdf$gender)$residuals * (vdf$gender == "Male")
dvs.female <- lm(dvs ~ vdf$gender)$residuals * (vdf$gender == "Female")

# get design
design_mat <- cbind(pmale=probs.male, pfemale=probs.female, 
                    dmale=dvs.male, dfemale=dvs.female, 
                    male=(vdf$gender=="Male")*1, female=(vdf$gender=="Female")*1)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
probs  <- as.matrix(df.xmats[,c("pmale", "pfemale")])
sfits1 <- get.fits(y ~ faces + quests + motion + probs + male, df.xmats)
df.sfits1 <- get.sfits(sfits1)
## decision vals
dvals  <- as.matrix(df.xmats[,c("dmale", "dfemale")])
sfits2 <- get.fits(y ~ faces + quests + motion + dvals + male, df.xmats)
df.sfits2 <- get.sfits(sfits2)

# probs
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

# dvals
fdr.tthr1 <- 2.7
sdf <- subset(df.sfits2, measure != "faces" & hemi == "rh")
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



# Compare Classifiers -----------------------------------------------------

# accuracy is fairly similar
knnFit$results[rownames(knnFit$bestTune),]
svmFit$results[rownames(svmFit$bestTune),]

# get probs/preds
knn.probs <- extract_probs(knnFit)
knn.preds <- extract_preds(knnFit)
svm.probs <- extract_probs(svmFit)
svm.preds <- extract_preds(svmFit)

# cross...table (mostly the same)
table(knn.preds$pred, svm.preds$pred)
cor(knn.probs$Female, svm.probs$Female) # and highly related (strange)
cor(lm(knn.probs$Female ~ vdf$gender)$residuals, 
    lm(svm.probs$Female ~ vdf$gender)$residuals) # less so if remove main effect of gender


# Combine -----------------------------------------------------------------

# Can we determine which of the measures is important?

# kNN
mean.preds <- extract_probs(knnFit)
mode.preds <- extract_preds(knnFit)
probs.male   <- lm(mean.preds$Male ~ vdf$gender)$residuals * (vdf$gender == "Male")
probs.female <- lm(mean.preds$Female ~ vdf$gender)$residuals * (vdf$gender == "Female")

# svm
mean.preds <- extract_probs(svmFit)
mode.preds <- extract_preds(svmFit)
dvs <- svmFit$finalModel$decision.values
probs.male2   <- lm(mean.preds$Male ~ vdf$gender)$residuals * (vdf$gender == "Male")
probs.female2 <- lm(mean.preds$Female ~ vdf$gender)$residuals * (vdf$gender == "Female")

# combine
design_mat <- cbind(avg=dist.to.avg, avg.male=dist.to.male2, avg.female=dist.to.female2,
                    gender.m=(vdf$gender=="Male")*1, gender.f=(vdf$gender=="Female")*1, 
                    knn.probs.m=probs.male, knn.probs.f=probs.female, 
                    svm.probs.m=probs.male2, svm.probs.f=probs.female2, 
                    svm.dmale=dvs.male, svm.dfemale=dvs.female)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion  <- as.matrix(df.xmats[,5:10])
pdist   <- as.matrix(df.xmats[,c("avg.male", "avg.female")])
kprobs  <- as.matrix(df.xmats[,c("knn.probs.m", "knn.probs.f")])
sprobs  <- as.matrix(df.xmats[,c("svm.probs.m", "svm.probs.f")])
sdists  <- as.matrix(df.xmats[,c("svm.dmale", "svm.dfemale")])

###
# Do each set individually
###

sfits1 <- get.fits(y ~ faces + quests + motion + gender.m + pdist, df.xmats)
sfits2 <- get.fits(y ~ faces + quests + motion + gender.m + kprobs, df.xmats)
sfits3 <- get.fits(y ~ faces + quests + motion + gender.m + sprobs, df.xmats)
sfits4 <- get.fits(y ~ faces + quests + motion + gender.m + sdists, df.xmats)

# So across all see that it's knn, prototype, svm
anova(sfits1[[5]]) # F = 11
anova(sfits2[[5]]) # F = 23
anova(sfits3[[5]]) # F = 7.8
anova(sfits4[[5]]) # F = 1

anova(sfits1[[2]]) # F = 13
anova(sfits2[[2]]) # F = 18
anova(sfits3[[2]]) # F = 7
anova(sfits4[[2]]) # F = 2

###
# Plot the f-stats
###

df.fs1 <- data.frame(
  hemi = "Right Hemi", 
  roi = rep(c("OFA", "FFA"), each=4), 
  measure = rep(c("Protoype: Dist", "Exemplar: Prob", "Category: Prob", "Category: Dist"), 2), 
  fval = c(
  c(anova(sfits1[[2]])[6,3], anova(sfits2[[2]])[6,3], anova(sfits3[[2]])[6,3], anova(sfits4[[2]])[6,3]), 
  c(anova(sfits1[[5]])[6,3], anova(sfits2[[5]])[6,3], anova(sfits3[[5]])[6,3], anova(sfits4[[5]])[6,3])
  )
)
df.fs2 <- data.frame(
  hemi = "Left Hemi", 
  roi = rep(c("OFA", "FFA"), each=4), 
  measure = rep(c("Protoype: Dist", "Exemplar: Prob", "Category: Prob", "Category: Dist"), 2), 
  fval = c(
    c(anova(sfits1[[12]])[6,3], anova(sfits2[[12]])[6,3], anova(sfits3[[12]])[6,3], anova(sfits4[[12]])[6,3]), 
    c(anova(sfits1[[15]])[6,3], anova(sfits2[[15]])[6,3], anova(sfits3[[15]])[6,3], anova(sfits4[[15]])[6,3])
  )
)
df.fs <- rbind(df.fs1, df.fs2)
df.fs$measure <- factor(df.fs$measure, levels=c("Protoype: Dist", "Exemplar: Prob", "Category: Prob", "Category: Dist"))

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2^2
ggplot(df.fs, aes(x=measure, y=fval, fill=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("F-Stats") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(hemi~roi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1))



###
# Plot the t-stats
###

# summarize
df.sfits1 <- get.sfits(sfits1)
df.sfits2 <- get.sfits(sfits2)
df.sfits3 <- get.sfits(sfits3)
df.sfits4 <- get.sfits(sfits4)
df.sfits  <- rbind(
  subset(df.sfits1, !(measure %in% c("faces","gender.m"))),
  subset(df.sfits2, !(measure %in% c("faces","gender.m"))), 
  subset(df.sfits3, !(measure %in% c("faces","gender.m"))), 
  subset(df.sfits4, !(measure %in% c("faces","gender.m")))
)
df.sfits <- df.sfits[order(df.sfits$hemi, df.sfits$roi, df.sfits$ord, df.sfits$measure),]
df.sfits$measure <- factor(df.sfits$measure)
df.sfits$mname   <- factor(df.sfits$measure, levels=levels(df.sfits$measure), 
                           labels=c("Prototype: Dist", "Prototype: Dist", 
                                    "Exemplar: Prob", "Exemplar: Prob", 
                                    "Category: Prob", "Category: Prob", 
                                    "Category: Dist", "Category: Dist"))
df.sfits$contrast <- factor(df.sfits$measure, levels=levels(df.sfits$measure), 
                           labels=c("Female", "Male", "Female", "Male", 
                                    "Female", "Male", "Female", "Male"))
df.sfits$mname <- factor(df.sfits$mname)
df.sfits$contrast <- factor(df.sfits$contrast)

# select only OFA and FFA
df.sfits <- subset(df.sfits, roi %in% c("OFA", "mFFA-2"))
df.sfits$roi <- factor(df.sfits$roi)
levels(df.sfits$roi)[2] <- "FFA"
levels(df.sfits$hemi) <- c("Right Hemi", "Left Hemi")

# plot
ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.sfits, aes(x=mname, y=tval, fill=mname, alpha=contrast)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.5,1)) + 
  geom_hline(yintercept=c(1,-1)*fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("T-Values") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(hemi~roi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1))


###
# Combine
###

# ANOVA

library(lme4)
library(lmerTest)

get.fits4 <- function(f, df.xmats) {
  sfits <- llply(1:ncol(ts.mats), function(ri) {
    cat(ri, "\n")
    
    fit = lmer(f, data=cbind(y=ts.mats[,ri], df.xmats))
    sfit <- summary(fit)
    sfit
  }, .parallel=T)
  names(sfits) <- rnames
  
  sfits
}

sfits1 <- get.fits4(y ~ faces + quests + motion + gender.m + kprobs + pdist + sprobs + sdists + (1|subject/run), df.xmats)
sfits2 <- get.fits(y ~ faces + quests + motion + gender.m + kprobs + sprobs + pdist + sdists, df.xmats)
sfits3 <- get.fits(y ~ faces + quests + motion + gender.m + kprobs + sdists + sprobs + pdist, df.xmats)
sfits4 <- get.fits(y ~ faces + quests + motion + gender.m + kprobs + sdists + pdist + sprobs, df.xmats)
sfits5 <- get.fits(y ~ faces + quests + motion + gender.m + pdist + kprobs + sprobs + sdists, df.xmats)

ri <- 5
#kprobs + pdist + sprobs + sdists + 
fit1 = lmer(y ~ faces + quests + motion + gender.m + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit2 = lmer(y ~ faces + quests + motion + gender.m + kprobs + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit3 = lmer(y ~ faces + quests + motion + gender.m + kprobs + pdist + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit4 = lmer(y ~ faces + quests + motion + gender.m + kprobs + sprobs + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit5 = lmer(y ~ faces + quests + motion + gender.m + kprobs + sdists + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
anova(fit1, fit2, fit3, fit4, fit5)
anova(fit1, fit2, fit5)

ri <- 
#kprobs + pdist + sprobs + sdists + 
fit1 = lmer(y ~ faces + quests + motion + gender.m + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit2 = lmer(y ~ faces + quests + motion + gender.m + kprobs + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit3 = lmer(y ~ faces + quests + motion + gender.m + kprobs + pdist + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit4 = lmer(y ~ faces + quests + motion + gender.m + kprobs + sprobs + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
fit5 = lmer(y ~ faces + quests + motion + gender.m + kprobs + sprobs + sdists + (1|subject/run), data=cbind(y=ts.mats[,ri], df.xmats))
anova(fit1, fit2)
anova(fit1, fit2, fit3)
anova(fit1, fit2, fit4)
anova(fit1, fit2, fit4, fit5)



anova(sfits1[[2]]) # most parsimonious
anova(sfits2[[2]])
anova(sfits3[[2]])
anova(sfits4[[2]])

anova(sfits1[[5]]) # close
anova(sfits2[[5]])
anova(sfits3[[5]])
anova(sfits4[[5]]) # most parsimonious

# combine all of them together

df.fs1 <- data.frame(
  hemi = "Right Hemi", 
  roi = rep(c("OFA", "FFA"), each=5), 
  measure = rep(c("Gender", "Exemplar: Prob", "Protoype: Dist", "Category: Prob", "Category: Dist"), 2), 
  fval = c(anova(sfits1[[2]])[5:9,3], anova(sfits1[[5]])[5:9,3])
)

df.fs2 <- data.frame(
  hemi = "Right Hemi", 
  roi = rep(c("OFA", "FFA"), each=5), 
  measure = rep(c("Gender", "Exemplar: Prob", "Protoype: Dist", "Category: Prob", "Category: Dist"), 2), 
  fval = c(anova(sfits1[[12]])[5:9,3], anova(sfits1[[15]])[5:9,3])
)

df.fs <- rbind(df.fs1, df.fs2)
df.fs$measure <- factor(df.fs$measure, levels=c("Gender", "Exemplar: Prob", "Protoype: Dist", "Category: Prob", "Category: Dist"))

# plot
ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2^2
ggplot(df.fs, aes(x=measure, y=fval, fill=measure)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.5,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("F-Stats") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(hemi~roi) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1))






# OTHER

get.fits0 <- function(f, df.xmats) {
  sfits <- llply(1:ncol(ts.mats), function(ri) {
    cat(ri, "\n")
    
    fit = lme(f, random = ~ 1|subject/run, 
              data=cbind(y=ts.mats[,ri], df.xmats), 
              control = lmeControl(opt = "optim"))
    fit
  }, .parallel=T)
  names(sfits) <- rnames
  
  sfits
}




fits1 <- get.fits0(y ~ faces + quests + motion + gender.m, df.xmats)
fits1b <- get.fits0(y ~ faces + quests + motion + gender.m + kprobs, df.xmats)

ri <- 5
fit1 = lme(y ~ faces + quests + motion + gender.m, random = ~ 1|subject/run, 
          data=cbind(y=ts.mats[,ri], df.xmats), 
          control = lmeControl(opt = "optim"))
fit1a = lme(y ~ faces + quests + motion + gender.m + pdist, random = ~ 1|subject/run, 
           data=cbind(y=ts.mats[,ri], df.xmats), 
           control = lmeControl(opt = "optim"))
fit1b = lme(y ~ faces + quests + motion + gender.m + kprobs, random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"))
anova(fit1, fit1a)
anova(fit1, fit1b)

library(MASS)
fit1 = lme(y ~ faces + quests + motion + gender.m + pdist + kprobs + sprobs + sdists, random = ~ 1|subject/run, 
           data=cbind(y=ts.mats[,ri], df.xmats), 
           control = lmeControl(opt = "optim"), method="ML")
aic.fits <- stepAIC(fit1, scope=list(upper=~faces + quests + motion + gender.m + kprobs + pdist + sprobs + sdists, lower=~faces + quests + motion + gender.m))

fit2 = lme(y ~ faces + quests + motion + gender.m, random = ~ 1|subject/run, 
           data=cbind(y=ts.mats[,ri], df.xmats), 
           control = lmeControl(opt = "optim"), method="ML")
aic.fits <- stepAIC(fit2, scope=list(upper=~faces + quests + motion + gender.m + kprobs + pdist + sprobs + sdists, lower=~faces + quests + motion + gender.m), direction="forward")

# TO FINISH: this will get the order of adding everything
ri <- 15

base <- "faces + quests + motion + gender.m"
lst.addtl <- c("kprobs", "pdist", "sprobs", "sdists")
seq.addtl <- c()

for (i in 1:length(lst.addtl)) {
  scope.f <- sprintf("~ %s + %s", base, paste(lst.addtl, collapse=" + "))
  f <- sprintf("y ~ %s", base)
  cat(f, "\n")
  cat(scope.f, "\n")
  
  fit = lme(as.formula(f), random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"), method="ML")
  
  ret <- add1(fit, as.formula(scope.f))
  mret <- as.matrix(ret)
  ind <- which.min(mret[-1,2])+1
  add.term <- rownames(mret)[ind]
  
  base <- paste(base, add.term, sep=" + ")
  lst.addtl <- lst.addtl[lst.addtl != add.term]
  seq.addtl <- c(seq.addtl, add.term)
}

# Then a STEP AIC for going back
base <- "faces + quests + motion + gender.m"
f <- sprintf("y ~ %s + %s", base, paste(seq.addtl, collapse = " + "))
fit = lme(as.formula(f), random = ~ 1|subject/run, 
           data=cbind(y=ts.mats[,ri], df.xmats), 
           control = lmeControl(opt = "optim"), method="ML")
upper.f <- sprintf("~ %s + %s", base, paste(seq.addtl, collapse = " + "))
lower.f <- sprintf("~ %s", base)
aic.fits <- stepAIC(fit, scope = list(upper=as.formula(upper.f), lower=as.formula(lower.f)), direction="both")
# final model!
anova(fit)
anova(aic.fits)

# use stepAIC to determine the final set of models

print(ret)

# base model
aic.fits$terms[[3]]

# this will give the order to add the additional things to
sub("- ", "", as.character(aic.fits$anova$Step)[-1])

drop1(fit1, test="Chisq")


sfits1a <- get.fits(y ~ faces + quests + motion + gender.m + pdist + kprobs, df.xmats)
sfits1b <- get.fits(y ~ faces + quests + motion + gender.m + kprobs + pdist, df.xmats)



# what if we combine the prototype with knn
sfits1 <- get.fits(y ~ faces + quests + motion + gender.m + pdist + kprobs, df.xmats)
sfits1.2 <- get.fits(y ~ faces + quests + motion + gender.m + kprobs + pdist, df.xmats)
sfits2 <- get.fits(y ~ faces + quests + motion + gender.m + kprobs + sprobs, df.xmats)
sfits2.2 <- get.fits(y ~ faces + quests + motion + gender.m + sprobs + kprobs, df.xmats)
sfits3 <- get.fits(y ~ faces + quests + motion + gender.m + avg.female + kprobs, df.xmats)
sfits4 <- get.fits(y ~ faces + quests + motion + gender.m + pdist + kprobs + sprobs, df.xmats)

## k-probs better if include first (less need to include prototype effect)
anova(sfits1[[5]])
anova(sfits1.2[[5]])

## k-probs also better (no need to include svm)
anova(sfits2[[5]])
anova(sfits2.2[[5]])

sfits1[[5]]$tTable
sfits2[[5]]$tTable
sfits3[[5]]$tTable #nope....
sfits4[[5]]$tTable

# maybe plot the Fstats or tstats for different elements
