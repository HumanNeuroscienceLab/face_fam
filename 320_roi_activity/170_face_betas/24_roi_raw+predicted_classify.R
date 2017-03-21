
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

## Time-Series Data
#load("/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda", verbose=T)

# Beta Data
load("/data1/famface01/analysis/misc/320_roi_task_activity/dgamma_betas.rda", verbose=T)

# Rnames
rnames2 <- sub("[.]", " ", rnames)
rnames2 <- sub("^r", "R", rnames2)
rnames2 <- sub("^l", "L", rnames2)

# Features
df.demos <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_df-demos.csv")[,-c(1:2)]
vnames <- read.table("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_demo-vnames.txt")
vnames <- as.character(vnames[,1])
all.mat <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-mat.csv")[,-1]

# Center all.mat columns
all.grps <- as.character(read.table("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-groups.txt")[,1])
all.cmat <- all.mat
all.cmat[,all.grps %in% c("traits","age","makeup")] <- scale(all.cmat[,all.grps %in% c("traits","age","makeup")], center=T, scale=F)
all.smat <- all.mat
all.smat[,all.grps %in% c("traits","age","makeup")] <- scale(all.cmat[,all.grps %in% c("traits","age","makeup")], center=T, scale=T)

# Low-level features
pca.feats <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_sym-shape-texture.csv")

## REORDER
# dataframe
traits <- as.matrix(all.mat[,all.grps=="traits"])
df <- cbind(traits, df.demos)
inds <- sapply(bvnames, function(vname) which(vnames==vname))
all.equal(vnames[inds], bvnames)
df <- df[inds,]
## scale
dfs <- df
dfs[,c(1:7,9)] <- scale(dfs[,c(1:7,9)])
##
all.mat2 <- all.mat[inds,]
all.cmat2 <- all.cmat[inds,]
all.smat2 <- all.smat[inds,]
##
all.equal(as.character(pca.feats$X)[inds], bvnames)
pca.feats <- pca.feats[,-1]
pca.feats2 <- pca.feats[inds,]

# Predictor/Residuals
load("/data1/famface01/analysis/misc/320_roi_task_activity/41_predictors.rda", verbose=T)
## combine pred stuff together
colnames(res.multinomial$predprobs) <- colnames(res.multinomial$residprobs)
preds1 <- cbind(res.gaussian$predvals, res.multinomial$predprobs, res.binomial$predprobs)
preds2 <- cbind(res.gaussian$predvals, res.multinomial$predclass, res.binomial$predclass)
preds2 <- model.matrix(~., data=preds2)
## combine resid stuff together
resids1 <- cbind(res.gaussian$residvals, res.multinomial$residprobs, res.binomial$residprobs)
resids2 <- cbind(res.gaussian$residvals, res.multinomial$residclass, res.binomial$residclass)
## reorder
inds <- sapply(bvnames, function(vname) which(vnames==vname))
preds1 <- preds1[inds,]
preds2 <- preds2[inds,]
resids1 <- resids1[inds,]
resids2 <- resids2[inds,]



# Forward Classify ----------------------------------------------------------

source("/data1/famface01/command/misc/face_representations/misc/cv_glmnet.R")

###
# See which set of features needed
# - using the raw trait values
# 

# Go through each of the ROIs
raw.leaveout <- llply(1:length(rnames), function(ri) {
  # Full model
  X <- as.matrix(all.cmat2)
  y <- grp.bs1[,ri,2]
  raw.full <- run_repeated_cvglmnet(X, y, parallel=T, k=10, nreps=10)
  
  # Remove one set of features
  ## combine the age/gender and glasses/makeup
  all.grps2 <- all.grps
  #all.grps2[all.grps %in% c("age", "gender")] <- "age/gender"
  #all.grps2[all.grps %in% c("glasses", "makeup")] <- "glasses/makeup"
  ugrps <- unique(all.grps2)
  ## loop through
  raw.leaveout <- llply(ugrps, function(ugrp) {
    X <- as.matrix(all.cmat2[,all.grps2!=ugrp])
    y <- grp.bs1[,ri,2]
    tmp1 <- run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
    tmp1
  }, .parallel=T)
  names(raw.leaveout) <- ugrps
  
  ret <- rlist::list.prepend(raw.leaveout, full=raw.full)
  ret
}, .progress="text")
names(raw.leaveout) <- rnames

# Average across regions
barplot(sapply(raw.leaveout, function(x) x$full$mean.max.res), names.arg=rnames2)

# Now get the max.res vals throughout
fit.rsqs <- laply(raw.leaveout, function(x) sapply(x, function(xx) xx$max.res))
dimnames(fit.rsqs)[[1]] <- rnames
diffs <- apply(fit.rsqs[,,-1] - fit.rsqs[,,rep(1,dim(fit.rsqs)[3]-1)], c(1,3), mean)
sign.diffs <- apply(sign(fit.rsqs[,,-1] - fit.rsqs[,,rep(1,dim(fit.rsqs)[3]-1)]), c(1,3), mean)

# Make heatmap of the sign
heatmap(sign.diffs * (abs(sign.diffs) > 0.4), Colv=NA, Rowv=NA, col=rev(brewer.pal(11, "RdBu")), scale="none")
heatmap(diffs * (abs(sign.diffs) > 0.4), Colv=NA, Rowv=NA, col=rev(brewer.pal(11, "RdBu")), scale="none", zlim=c(-0.067, 0.067))

# For just traits, should we get the weightings?
sapply(raw.leaveout, function(x) x$full$coefs)
# Try the lassoscore to get the sig vals?


###
# Forward with just the traits and age
#

# TODO



###
# Compare to low-level features (FORWARD MODEL)

X <- as.matrix(cbind(all.cmat[,1:7], pca.feats2))
raw.fulls <- llply(1:length(rnames), function(i) {
  y <- grp.bs1[,i,2]
  run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
}, .parallel=T)
names(raw.fulls) <- rnames

X <- as.matrix(cbind(all.cmat[,1:7]))
raw.traits <- llply(1:length(rnames), function(i) {
  y <- grp.bs1[,i,2]
  run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
}, .parallel=T)
names(raw.traits) <- rnames

X <- as.matrix(cbind(pca.feats2))
raw.lows <- llply(1:length(rnames), function(i) {
  y <- grp.bs1[,i,2]
  run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
}, .parallel=T)
names(raw.lows) <- rnames

# traits shows more of a drop, meaning removing the low-level features isn't good
# so here this suggests that the low-level features are in fact better!
colMeans(sign(sapply(raw.traits, function(x) x$max.res) - sapply(raw.fulls, function(x) x$max.res)))
colMeans(sign(sapply(raw.lows, function(x) x$max.res) - sapply(raw.fulls, function(x) x$max.res)))



# Predicted/Residuals Classify -------------------------------------------

# 1: has the probability values for the binomial classifiers
rm.inds <- c(grep("facial_hair", colnames(preds1)), grep("^hair", colnames(preds1)), 
             grep("race", colnames(preds1)), grep("eye", colnames(preds1)))
sel.preds1 <- preds1[,-rm.inds][,1:7]
sel.resids1 <- resids1[,-rm.inds][,1:7]

# 2: has the class labels for the binomial classifiers
rm.inds <- c(grep("facial_hair", colnames(preds2)), grep("^hair", colnames(preds2)), 
             grep("race", colnames(preds2)), grep("eye", colnames(preds2)))
sel.preds2 <- preds2[,-rm.inds][,1:7]
rm.inds <- c(grep("facial_hair", colnames(resids2)), grep("^hair", colnames(resids2)), 
             grep("race", colnames(resids2)), grep("eye", colnames(resids2)))
sel.resids2 <- resids2[,-rm.inds][,1:7]


###
# What is better? Low-level or the predicted/residual scores 
#

# FULL MODEL
X <- as.matrix(cbind(sel.preds1, sel.resids1, pca.feats2))
pred.fulls <- llply(1:length(rnames), function(i) {
  y <- grp.bs1[,i,2]
  run_repeated_cvglmnet(X, y, parallel=F)
}, .parallel=T)
names(pred.fulls) <- rnames

# LEFT OUT FEATURES
X0 <- as.matrix(cbind(sel.preds1, sel.resids1, pca.feats2))
Xgrps <- rep(c("pred", "resid", "low"), c(ncol(sel.preds1), ncol(sel.resids1), ncol(pca.feats2)))
ugrps <- unique(Xgrps)
pred.leftout <- llply(1:length(rnames), function(i) {
  ret <- lapply(ugrps, function(grp) {
    X <- X0[,Xgrps!=grp]
    y <- grp.bs1[,i,2]
    run_repeated_cvglmnet(X, y, parallel=F)
  })
  names(ret) <- ugrps
  ret
}, .parallel=T)
names(pred.leftout) <- rnames

# INDIVIDUAL SETS
pred.indivs <- llply(1:length(rnames), function(i) {
  ret <- lapply(ugrps, function(grp) {
    X <- X0[,Xgrps==grp]
    y <- grp.bs1[,i,2]
    run_repeated_cvglmnet(X, y, parallel=F)
  })
  names(ret) <- ugrps
  ret
}, .parallel=T)
names(pred.indivs) <- rnames

## COMPARE R2 with raw
indiv.rsq <- cbind(
  raw=sapply(raw.traits, function(x) x$mean.max.res), 
  pred=sapply(pred.indivs, function(x) x$pred$mean.max.res), 
  resid=sapply(pred.indivs, function(x) x$resid$mean.max.res), 
  low=sapply(pred.indivs, function(x) x$low$mean.max.res)
)
round(indiv.rsq, 3)
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
heatmap(indiv.rsq, col=cols, scale="none", Colv=NA, Rowv=NA)

## OTHER STUFF
tmp1 <- laply(pred.fulls, function(x) x$max.res)
tmp2 <- laply(pred.leftout, function(x) sapply(x, function(xx) xx$max.res))
tmp.diff <- laply(1:3, function(i) tmp2[,,i] - tmp1)
dimnames(tmp.diff)[[1]] <- ugrps
dimnames(tmp.diff)[[2]] <- rnames
apply(tmp.diff, 1:2, function(x) mean(sign(x))) # consistent drop in aFFA, pFFA, EBA
apply(tmp.diff, 1:2, function(x) mean(x))




## I WANT TO ONLY LOOK AT PRED + LOW
X <- as.matrix(cbind(sel.preds1, pca.feats2))
pred.fulls2 <- llply(1:length(rnames), function(i) {
  y <- grp.bs1[,i,2]
  run_repeated_cvglmnet(X, y, parallel=F)
}, .parallel=T)
names(pred.fulls2) <- rnames
sapply(pred.fulls2, function(x) x$mean.max.res)


###
# Now I want to classify the predicted and residual data
# and visualize the coefs
# 

X <- as.matrix(cbind(sel.preds1, sel.resids1))
pred.resid <- llply(1:length(rnames), function(i) {
  y <- grp.bs1[,i,2]
  
  fit <- run_repeated_cvglmnet(X, y, parallel=F)
  
  fit.score  <- lassoscore(y, X, lambda=fit$lambda, family="gaussian")
  fit$lscore <- fit.score
  fit$zvals  <- qt(fit.score$p.model, Inf, lower.tail=F)
  names(fit$zvals) <- names(fit$coefs)
  
  fit  
}, .parallel=T)
names(pred.resid) <- rnames

coefs1 <- sapply(pred.resid, function(x) abs(x$coefs))
zvals1 <- sapply(pred.resid, function(x) x$zvals * (x$coefs!=0))
zvals2 <- sapply(pred.resid, function(x) x$zvals * (x$zvals>1.96))

colnames(coefs1) <- rnames2
colnames(zvals1) <- rnames2
colnames(zvals2) <- rnames2

pinds <- 1:7
rinds <- 8:nrow(coefs1)


## HEATMAPS

quick.heat <- function(mat, ...) {
  cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
  seplines <- sum(grepl("^R", rnames2)) + 0.5
  
  tmat <- mat * (mat>0) # make sure no negatives
  zlim <- c(0,max(tmat))
  tmat[tmat==0] <- NA
  
  heatmap(tmat, scale="none", col=cols, 
          zlim=zlim, Rowv=NA, Colv=NA, 
          margins=c(5,0), cexCol = 1.25, cexRow = 1.25, 
          add.expr = abline(v=seplines, lwd=4, lty=3, col="grey70"), 
          ...)
}

quick.heat(coefs1[pinds,], main="Predicted Factors - abs(coefs)")
quick.heat(coefs1[rinds,], main="Residual Factors - abs(coefs)")

quick.heat(zvals1[pinds,], main="Predicted Factors - zvals*coefs")
quick.heat(zvals1[rinds,], main="Residual Factors - zvals*coefs")

quick.heat(zvals2[pinds,], main="Predicted Factors - abs(zvals2)")
quick.heat(zvals2[rinds,], main="Residual Factors - abs(zvals2)")

# TODO: regression on left-over features




###
# Invert the model to see the effects
# - using the raw trait values
# 

library(lassoscore)

raw.inverse <- llply(1:8, function(i) {
  X <- as.matrix(grp.bs1[,,2])
  y <- all.cmat2[,i]
  
  fit <- run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
  
  fit.score <- lassoscore(y, X, lambda=fit$lambda, family="gaussian")
  fit$lscore <- fit.score
  fit$zvals  <- qt(fit.score$p.model, Inf, lower.tail=F)
  
  fit
}, .parallel=T)
names(raw.inverse) <- colnames(all.cmat2)[1:8]

# fit rsqs
barplot(sapply(raw.inverse, function(x) x$mean.max.res))

# get the thresholded coefs
thr.zvals1 <- sapply(raw.inverse, function(x) x$zvals * (x$coefs!=0))
thr.zvals2 <- sapply(raw.inverse, function(x) x$zvals * (x$zvals>1.96))
rownames(thr.zvals1) <- rnames
rownames(thr.zvals2) <- rnames



###
# Run inverse using the predicted trait values
# 

library(corrplot)

# visualize the brain data [first]
cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
seplines <- sum(grepl("^R", rnames2)) + 0.5
heatmap(grp.bs1[,,2], Colv=NA, scale="none", labRow=NA, col=cols, labCol=rnames2, 
        add.expr = abline(v=seplines, lwd=4, lty=3, col="grey70"))

pred.inverse <- llply(1:8, function(i) {
  X <- as.matrix(grp.bs1[,,2])
  y <- preds1[,i]
  
  fit <- run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
  
  fit.score <- lassoscore(y, X, lambda=fit$lambda, family="gaussian")
  fit$lscore <- fit.score
  fit$zvals  <- qt(fit.score$p.model, Inf, lower.tail=F)
  
  fit
}, .parallel=T)
names(pred.inverse) <- colnames(preds1)[1:8]

resid.inverse <- llply(1:8, function(i) {
  X <- as.matrix(grp.bs1[,,2])
  y <- resids1[,i]
  
  fit <- run_repeated_cvglmnet(X, y, parallel=F, k=10, nreps=10)
  
  fit.score <- lassoscore(y, X, lambda=fit$lambda, family="gaussian")
  fit$lscore <- fit.score
  fit$zvals  <- qt(fit.score$p.model, Inf, lower.tail=F)
  
  fit
}, .parallel=T)
names(resid.inverse) <- colnames(resids1)[1:8]

# fit rsqs
barplot(sapply(pred.inverse, function(x) x$mean.max.res))

mat <- rbind(raw=sapply(raw.inverse, function(x) x$mean.max.res), 
             pred=sapply(pred.inverse, function(x) x$mean.max.res), 
             resid=sapply(resid.inverse, function(x) x$mean.max.res))
barplot(mat, beside=T, legend=T, col=brewer.pal(3, "Spectral"))

# get the thresholded coefs
thr.zvals1 <- sapply(raw.inverse, function(x) x$zvals * (x$coefs!=0))
thr.zvals2 <- sapply(raw.inverse, function(x) x$zvals * (x$zvals>1.96))
rownames(thr.zvals1) <- rnames2
rownames(thr.zvals2) <- rnames2

quick.heat(t(thr.zvals1), main="Raw Factors - zvals*coefs")
quick.heat(t(thr.zvals2), main="Raw Factors - zvals")

quick.heat(zvals2[pinds,], main="Predicted Factors - abs(zvals2)")
quick.heat(zvals2[rinds,], main="Residual Factors - abs(zvals2)")
