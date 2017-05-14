
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



# Regressors --------------------------------------------------------------

library(caret)
library(e1071)
registerDoMC(30)
X <- vfeats; y <- vdf$gender

###
# Prototype Distances
###

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

###
# Exemplar Probs
###

knnFit <- run_caret(X, y, "knn", tlen=20)
mean.preds <- extract_probs(knnFit)
mode.preds <- extract_preds(knnFit)
probs.male   <- lm(mean.preds$Male ~ vdf$gender)$residuals * (vdf$gender == "Male")
probs.female <- lm(mean.preds$Female*1 ~ vdf$gender)$residuals * (vdf$gender == "Female")

###
# Category Probs + Dists
###

svmFit <- run_caret(X, y, "svmLinear2", tlen=10)
# svm probs
mean.preds <- extract_probs(svmFit)
mode.preds <- extract_preds(svmFit)
probs.male2   <- lm(mean.preds$Male ~ vdf$gender)$residuals * (vdf$gender == "Male")
probs.female2 <- lm(mean.preds$Female ~ vdf$gender)$residuals * (vdf$gender == "Female")
# svm dists (higher value, the more of that category you should be)
dvs <- svmFit$finalModel$decision.values
dvs.male   <- lm(dvs ~ vdf$gender)$residuals * (vdf$gender == "Male")
dvs.female <- lm((dvs*-1) ~ vdf$gender)$residuals * (vdf$gender == "Female")

###
# Combine
###

design_mat <- cbind(gender.m=(vdf$gender=="Male")*1, gender.f=(vdf$gender=="Female")*1, 
                    gender.f_gt_m=((vdf$gender=="Female")*2 - 1), 
                    avg.male=dist.to.male3, avg.female=dist.to.female3,
                    knn.probs.m=probs.male, knn.probs.f=probs.female, 
                    svm.probs.m=probs.male2, svm.probs.f=probs.female2, 
                    svm.dmale=dvs.male, svm.dfemale=dvs.female)

design_mat2 <- cbind(gender.f=(vdf$gender=="Female")*1, 
                    avg=dist.to.male3 + dist.to.female3, 
                    knn.probs=probs.male + probs.female, 
                    svm.probs=probs.male2 + probs.female2, 
                    svm.dmale=dvs.male + dvs.female)


# Plot weights
library(corrplot)
cmat <- cor(design_mat)
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                             "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                             "#4393C3", "#2166AC", "#053061"))(200))
corrplot(cmat, tl.col='black', tl.cex=.75, diag=F, col=col, is.corr=T) 

# Convolve
df.xmats <- gen.xmat(design_mat)
df.xmats2 <- gen.xmat(design_mat2)




# Regressionz -------------------------------------------------------------

# stepwise regression
motion  <- as.matrix(df.xmats[,5:10])
pdist   <- as.matrix(df.xmats[,c("avg.male", "avg.female")])
kprobs  <- as.matrix(df.xmats[,c("knn.probs.m", "knn.probs.f")])
sprobs  <- as.matrix(df.xmats[,c("svm.probs.m", "svm.probs.f")])
sdists  <- as.matrix(df.xmats[,c("svm.dmale", "svm.dfemale")])

library(MASS)
lst.fits <- list()
for (ri in 1:ncol(ts.mats)) {
  cat("======\n")
  cat(rnames[ri], "\n")
  
  ## Get the right order of the full model
  base <- "faces + quests + motion + gender.m"
  lst.addtl <- c("pdist", "kprobs", "sprobs", "sdists")
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
  
  ## Then do backwards selection for best subset model
  cat("backwards\n")
  base <- "faces + quests + motion + gender.m"
  f <- sprintf("y ~ %s + %s", base, paste(seq.addtl, collapse = " + "))
  fit = lme(as.formula(f), random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats), 
            control = lmeControl(opt = "optim"), method="ML")
  upper.f <- sprintf("~ %s + %s", base, paste(seq.addtl, collapse = " + "))
  lower.f <- sprintf("~ %s", base)
  aic.fits <- stepAIC(fit, scope = list(upper=as.formula(upper.f), lower=as.formula(lower.f)), direction="both")
  
  lst.fits[[ri]] <- list(full=fit, subset=aic.fits)
}
names(lst.fits) <- rnames
  
anova(lst.fits$rh.OFA$full)
anova(lst.fits$`rh.mFFA-1`$full)

summary(lst.fits$rh.OFA$subset)$tTable
summary(lst.fits$`rh.mFFA-1`$subset)$tTable

anova(lst.fits$lh.OFA$full)
anova(lst.fits$`lh.mFFA-1`$full)


ret <- stepwise.wrapper(5)

###
# Plot the f-stats
###

selrows <- c("gender.m", "pdist", "kprobs", "sprobs", "sdists")
selnames<- c("Gender", "Prototype: Dist", "Exemplar: Prob", "Category: Prob", "Category: Dist")
selrinds<- c(2,5,12,15)
df.fs <- ldply(selrinds, function(ri) {
  ret <- anova(lst.fits[[ri]]$full)
  inds <- sapply(selrows, function(x) which(rownames(ret) == x))
  data.frame(hemi=sub("[.].*", "", rnames[ri]), 
             roi=sub("[rl]h[.]", "", rnames[ri]), 
             measure=selnames, fstat=ret[inds,3])
})
df.fs$measure <- factor(df.fs$measure, levels=selnames)

ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2^2
ggplot(df.fs, aes(x=measure, y=fstat, fill=measure)) + 
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


ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2^2
ggplot(df.fs, aes(x=measure, y=fstat, fill=measure, alpha=hemi)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.5,1)) + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') +
  geom_hline(yintercept=0, col="black") + 
  ylab("F-Stats") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  facet_grid(roi~.) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1))


###
# Plot the t-stats
###

selrows <- c("gender.m", "pdist", "kprobs", "sprobs", "sdists")
sorig <- c("gender.m", "pdistavg.male", "pdistavg.female", 
           "kprobsknn.probs.m", "kprobsknn.probs.f", 
           "sprobssvm.probs.m", "sprobssvm.probs.f",
           "sdistssvm.dmale", "sdistssvm.dfemale")
snew <- c("Gender", "Prototype", "Prototype", "Exemplar", "Exemplar", 
          "Category", "Category", "Categorical Boundary", "Categorical Boundary")
snew2 <- c("Male", rep(c("Male", "Female"), 4))
selnames<- c("Gender", "Prototype: Dist", "Exemplar: Prob", "Category: Prob", "Category: Dist")
selrinds<- c(2,5,12,15)
df.ts <- ldply(selrinds, function(ri) {
  ret <- summary(lst.fits[[ri]]$full)$tTable
  inds <- unlist(lapply(selrows, function(x) grep(x, rownames(ret))))
  vec  <- ret[inds,4]
  data.frame(hemi=sub("[.].*", "", rnames[ri]), 
             roi=sub("[rl]h[.]", "", rnames[ri]), 
             mname=names(vec), measure=factor(snew, levels=snew), gender=snew2, 
             tstat=vec)
})


ggthemr('pale', type='outer', text_size=16, layout='plain')
fdr.tthr1 <- 2.7
ggplot(df.ts, aes(x=measure, y=tstat, fill=measure, alpha=gender)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Spectral") + 
  scale_alpha_discrete(range=c(0.75,1)) + 
  geom_hline(yintercept=c(-1,1)*fdr.tthr1, linetype='dotted') +
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



zfit <- get.fits(y ~ faces + quests + motion + gender.m + kprobs, df.xmats)
zfit[[5]]$tTable

zfit <- get.fits(y ~ faces + quests + motion + gender.m + pdist, df.xmats)
zfit[[2]]$tTable
