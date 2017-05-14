# This script will load the AFNI prototype/exemplar model as well as the ROI data
# and then apply an LME



# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(doMC)
registerDoMC(30)
library(plyr)
library(nlme)

subjects <- sprintf("sub%02i", 1:6)



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




# Load regressors ---------------------------------------------------------

# Function to load AFNI matrices
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

taskdir <- "/data1/famface01/analysis/task_activity"

# load the xfiles for each subject
df.xmats <- ldply(subjects, function(subj) {
  xfile <- sprintf("%s/%s/face_deviations_unfam/nnet2_pca_avg+knn+covars.reml/xmat.1D", taskdir, subj)
  xmat <- read.xmat(xfile, rm.nums=T)
  # change the runs to one column indicating the run
  rinds <- grep("^Run.*Pol", colnames(xmat))
  rnums <- rowSums(sweep(xmat[,rinds], 2, 1:length(rinds), FUN="*"))
  runs  <- sprintf("run%02i", rnums)
  # add the subject and runs
  xmat <- cbind(subject=subj, run=runs, xmat[,-rinds])
  xmat
}, .progress="text")

# traits
df.xmats2 <- ldply(subjects, function(subj) {
  xfile <- sprintf("%s/%s/face_basics_unfam/raw_demo_trait_feats_sm4.reml/xmat.1D", taskdir, subj)
  xmat <- read.xmat(xfile, rm.nums=T)
  # change the runs to one column indicating the run
  rinds <- grep("^Run.*Pol", colnames(xmat))
  rnums <- rowSums(sweep(xmat[,rinds], 2, 1:length(rinds), FUN="*"))
  runs  <- sprintf("run%02i", rnums)
  # add the subject and runs
  xmat <- cbind(subject=subj, run=runs, xmat[,-rinds])
  xmat
}, .progress="text")



# Run regression ----------------------------------------------------------

# get the lme results
sfits <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + avg_dist + knn_dist + pose + frame_diff + frame_pose + roll + pitch + yaw + dS + dL + dP,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

# put them in table form
df.sfit <- ldply(rnames, function(rname) {
  ts <- sfits[[rname]]$tTable[,4]
  zs <- qt(sfits[[rname]]$tTable[,5], Inf, lower.tail=F)
  ps <- sfits[[rname]]$tTable[,5]
  hemi <- sub("[.].*", "", rname)
  name <- sub("[lr]h.", "", rname)
  ord  <- rords[rnames==rname]
  sdf <- data.frame(hemi=hemi, roi=name, ord=ord, 
                    measure=rownames(sfits[[rname]]$tTable), 
                    tval=ts, pval=ps, zval=zs)
  sdf[-1,] # remove intercept
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

#ret <- fdrtool(c(0.05,df.sfit$pval), statistic="pvalue")
#tmp <- p.adjust(c(0.05,df.sfit$pval), method="fdr")



# Viz ---------------------------------------------------------------------

library(ggplot2)
library(ggthemr)
library(RColorBrewer)

ggthemr('pale', type='outer', text_size=14, layout='plain')

df.sfit2 <- subset(df.sfit, measure == "faces" & hemi == "rh")

ggplot(df.sfit2, aes(x=roi, y=abs(tval), fill=roi)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=0) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5), 
        legend.position = "none")

ggplot(df.sfit2, aes(x=ord, y=abs(tval))) + 
  geom_line() + 
  geom_point(aes(color=roi), size=3) + 
  scale_color_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5), 
        legend.position = "none")


df.sfit2 <- subset(df.sfit, measure %in% c("avg_dist", "knn_dist") & hemi == "rh")
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

df.sfit2 <- subset(df.sfit, measure %in% c("avg_dist", "knn_dist") & hemi == "lh")
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

df.sfit2 <- subset(df.sfit, measure %in% c("pose", "frame_pose", "frame_diff") & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(size=3, shape=21) + 
  #scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))


# Run regression for traits -------------------------------------------------

# get the lme results
sfits2 <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + raw_demo_gender + raw_demo_age + roll + pitch + yaw + dS + dL + dP,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats2))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits2) <- rnames


sfits3 <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + raw_t01_unemotional + raw_t02_competent + raw_t03_trustworthy + raw_t04_memorable + raw_t05_attractive + roll + pitch + yaw + dS + dL + dP,
            random = ~ 1|subject/run, 
            data=cbind(y=ts.mats[,ri], df.xmats2))
  sfit <- summary(fit)
  sfit
}, .parallel=F)
names(sfits3) <- rnames

# put them in table form
df.sfit <- ldply(rnames, function(rname) {
  sfits <- sfits2
  ts <- sfits[[rname]]$tTable[,4]
  zs <- qt(sfits[[rname]]$tTable[,5], Inf, lower.tail=F)
  ps <- sfits[[rname]]$tTable[,5]
  hemi <- sub("[.].*", "", rname)
  name <- sub("[lr]h.", "", rname)
  ord  <- rords[rnames==rname]
  sdf <- data.frame(hemi=hemi, roi=name, ord=ord, 
                    measure=rownames(sfits[[rname]]$tTable), 
                    tval=ts, pval=ps, zval=zs)
  sdf[-1,] # remove intercept
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

df.sfit2 <- subset(df.sfit, measure %in% c("raw_demo_age", "raw_demo_gender") & hemi == "rh")
ggplot(df.sfit2, aes(x=ord, y=abs(tval), group=measure)) + 
  geom_line(aes(color=measure)) + 
  scale_color_hue() + 
  geom_point(size=3, shape=21) + 
  #scale_fill_brewer(palette = "Spectral") + 
  geom_hline(yintercept=fdr.tthr1, linetype='dotted') + 
  ylab("Absolute T-Value") + 
  scale_x_continuous(breaks=df.sfit2$ord, labels=df.sfit2$roi) + 
  scale_y_continuous(expand=c(0,0)) + 
  expand_limits(y=c(0,max(abs(df.sfit2$tval))*1.05)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_line(color="black", size=0.5))





# One-Class SVM  ----------------------------------------------------------

library(e1071)

# go in 2% increments to find the outliers
thrs <- seq(0.01, 0.99, 0.01)
rank.outliers <- laply(thrs, function(nu) {
  fit <- svm(vfeats, type="one-classification", scale=T, kernel="radial", 
             nu=nu, cost=1)
  pred <- predict(fit, decision.values=F, probability=F)
  (!pred) * 1
}, .parallel=T)
rank.outliers <- rowSums(t(rank.outliers)) + 1
table(rank.outliers) # lowest value is most in the center of cluster (core)

# this doesn't work
fit <- svm(vfeats, type="one-classification", scale=T, kernel="radial", 
           nu=0.5, cost=1)
w <- t(fit$coefs) %*% fit$SV;  
b <- -1 * fit$rho;
dvals <- ((w %*% t(vfeats)) + b) / as.numeric(sqrt(w %*% t(w)))
dvals <- as.vector(dvals)


# relate to atypicality?
summary(lm(vdf$typical ~ dist.to.avg + rank.outliers)) # HAHA IT WORKS!
summary(lm(vdf$attractive ~ dist.to.avg + rank.outliers)) # nope but dist works
summary(lm(vdf$typical ~ dist.to.avg + dvals)) # dvals on own ok
summary(lm(vdf$typical ~ dist.to.avg + rank.outliers + dvals)) # but not as good...

# get design
design_mat <- cbind(avg=dist.to.avg, outlier=rank.outliers, dvals=dvals, 
                    typical=vdf$typical)
round(cor(design_mat), 3)
df.xmats <- gen.xmat(design_mat)

# get the fits
motion <- as.matrix(df.xmats[,5:10])
sfits1 <- get.fits(y ~ faces + quests + motion + avg + outlier, df.xmats)
df.sfits1 <- get.sfits(sfits1)
## decision vals
sfits2 <- get.fits(y ~ faces + quests + motion + avg + outlier + dvals, df.xmats)
df.sfits2 <- get.sfits(sfits2)
## typical
sfits3 <- get.fits(y ~ faces + quests + motion + avg + typical, df.xmats)
df.sfits3 <- get.sfits(sfits3)


# 1
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

# 2
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

# 3
fdr.tthr1 <- 2.7
sdf <- subset(df.sfits3, measure != "faces" & hemi == "rh")
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






## classify

library(caret)
library(glmnet)
registerDoMC(24)

# subs.folds
subs.folds <- as.numeric(df.xmats$subject)
nrepeats <- 1
ufolds  <- sort(unique(subs.folds))
nfolds  <- length(ufolds)
runFolds <- createMultiFolds(ufolds, k=nfolds, times = nrepeats)
folds    <- lapply(runFolds, function(x) {
  which(subs.folds %in% x)
})

fitControl <- trainControl(
  method = "repeatedcv", 
  number = nfolds, 
  repeats = nrepeats, 
  index = folds, 
  returnResamp = "final", 
  savePredictions = "final", 
  classProbs = T, 
  allowParallel = T
)

alphas <- c(0,0.5,1)
y.resid <- lm(ts.mats ~ subject*run, data=df.xmats)$residuals
x.resid <- lm(as.matrix(df.xmats[,-c(1:2)]) ~ subject*run, data=df.xmats)$residuals

cvfits <- llply(1:ncol(y.resid), function(i) {
  tuneGrid <- ldply(alphas, function(alpha) {
    tmp <- glmnet(x.resid, y.resid[,i], family="gaussian", nlambda=50, alpha=alpha)
    data.frame(alpha=alpha, lambda=tmp$lambda)
  }, .parallel=T)
  
  fit1 <- train(x.resid, y.resid[,i], 
                method = "glmnet", 
                trControl = fitControl, 
                preProcess = c("center", "scale"), 
                tuneGrid = tuneGrid)
  fit1
}, .progress="text")

resfits <- ldply(1:length(cvfits), function(i) {
  cvfit <- cvfits[[i]]
  ri <- as.numeric(rownames(cvfit$bestTune))
  cbind(roi=rnames[i], cvfit$results[ri,])
})
resfits
vifits <- sapply(cvfits, function(cvfit) {
  varImp(cvfit)$importance$Overall
})
rownames(vifits) <- colnames(x.resid)
colnames(vifits) <- rnames
vifits[3:7,]
tmp <- varImp(cvfits[[1]])
tmp$importance

i <- 1
fit2 <- train(x.resid, y.resid[,i], 
              method = "gbm", 
              trControl = fitControl, 
              preProcess = c("center", "scale"), verbose=F)
ri <- as.numeric(rownames(fit2$bestTune))
fit2$results[ri,]
varImp(fit2)












library(PMA)
?PMA

x <- lm(ts.mats ~ subject*run, data=df.xmats)$residuals
z <- lm(as.matrix(df.xmats[,-c(1:2)]) ~ subject*run, data=df.xmats)$residuals

perm.out <- CCA.permute(x, z, typex="standard", typez="standard", nperms=25)
#print(perm.out)
out <- CCA(x,z,K=10,typex="standard",typez="standard",penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)
rownames(out$u) <- colnames(x)
rownames(out$v) <- colnames(z)
round(out$u,3)
round(out$v,3)

out2 <- CCA(x[df.xmats$subject!="sub06",],z[df.xmats$subject!="sub06",],K=10,typex="standard",typez="standard",penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)
rownames(out2$u) <- colnames(x)
rownames(out2$v) <- colnames(z)
round(out2$u,3)
round(out2$v,3)

out3 <- CCA(x[df.xmats$subject!="sub04",],z[df.xmats$subject!="sub04",],K=10,typex="standard",typez="standard",penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)
rownames(out3$u) <- colnames(x)
rownames(out3$v) <- colnames(z)
round(out3$u,3)
round(out3$v,3)



geom_line(aes(linetype=sex), size=1) +     # Set linetype by sex
  geom_point(size=3, fill="white") +         # Use larger points, fill with white

#df.sfit2 <- subset(df.sfit, measure %in% c("faces", "pose", "avg_dist", "knn_dist", "frame_pose", "frame_diff"))

tmp <- cancor(x, z)
tmp$xcoef
tmp$ycoef
+ 
  expand_limits(y=max(df.sfit2$tval)*1.05)
+ 
  coord_cartesian(ylim=c(95,100.5)) + 
  
  
  #scale_fill_manual(values=rev(brewer.pal(7, "Spectral")[1:3])) + 
  facet_wrap(~measure, scales = "free_y") + 
  ylab("Percent Accuracy") + 
  scale_y_continuous(breaks=95:100, expand=c(0,0)) + 
  expand_limits(y=100.5) + 
  coord_cartesian(ylim=c(95,100.5)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.line.x = element_line(color="black", size=1), 
        legend.position = "none")







# "/data1/famface01/analysis/task_activity/sub01/face_deviations_unfam/nnet2_pca_avg+knn+covars.reml/xmat.1D"
"/data1/famface01/analysis/task_activity/sub01/face_basics_unfam/raw_demo_trait_feats_sm4.reml/xmat.1D"


df.xmats2 <- ldply(subjects, function(subj) {
  xfile <- sprintf("%s/%s/face_basics_unfam/raw_demo_trait_feats_sm4.reml/xmat.1D", taskdir, subj)
  xmat <- read.xmat(xfile, rm.nums=T)
  # change the runs to one column indicating the run
  rinds <- grep("^Run.*Pol", colnames(xmat))
  rnums <- rowSums(sweep(xmat[,rinds], 2, 1:length(rinds), FUN="*"))
  runs  <- sprintf("run%02i", rnums)
  # add the subject and runs
  xmat <- cbind(subject=subj, run=runs, xmat[,-rinds])
  xmat
}, .progress="text")
head(df.xmats2)




head(df.xmats)
round(cor(df.xmats[,-c(1:2)]), 3) # avg and knn not correlated...

# try the LME
require(lmerTest)


sfits <- llply(1:ncol(ts.mats), function(ri) {
  fit <- lmer(ts.mats[,ri] ~ faces + quests + avg_dist + knn_dist + pose + frame_diff + frame_pose + roll + pitch + yaw + dS + dL + dP + (1|subject/run), data=df.xmats)
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits) <- rnames

coef(sfits$rh.V1)[-c(9:14),]
coef(sfits$rh.OFA)[-c(9:14),]
coef(sfits$`rh.pFFA-1`)[-c(9:14),]
coef(sfits$`rh.pFFA-2`)[-c(9:14),]
coef(sfits$`rh.mFFA-1`)[-c(9:14),]
coef(sfits$`rh.mFFA-2`)[-c(9:14),]
coef(sfits$`rh.aFFA-1`)[-c(9:14),]
coef(sfits$`rh.aFFA-2`)[-c(9:14),]
coef(sfits$`rh.vATL`)[-c(9:14),]

fit <- lmer(ts.mats[,1] ~ faces + quests + avg_dist + knn_dist + pose + frame_diff + frame_pose + roll + pitch + yaw + dS + dL + dP + (1|subject/run), data=df.xmats)
sfit <- summary(fit)
coef(sfit)[-c(9:14),]

fit <- lmer(ts.mats[,2] ~ faces + quests + avg_dist + knn_dist + pose + frame_diff + frame_pose + roll + pitch + yaw + dS + dL + dP + (1|subject/run), data=df.xmats)
sfit2 <- summary(fit)
coef(sfit2)[-c(9:14),]



coef(summary(fit))

fit = lme(ts.mats[,1] ~ faces + quests + avg_dist + knn_dist + pose + frame_diff + frame_pose + roll + pitch + yaw + dS + dL + dP,
              random = ~ 1|subject/run,
              data = df.xmats)
sfit <- summary(fit)
anova(fit)


sfits2 <- llply(1:ncol(ts.mats), function(ri) {
  cat(ri, "\n")
  fit = lme(y ~ faces + quests + avg_dist + knn_dist + pose + frame_diff + frame_pose + memorable + roll + pitch + yaw + dS + dL + dP,
              random = ~ 1|subject/run, 
             data=cbind(y=ts.mats[,ri], df.xmats, memorable=df.xmats2$raw_t04_memorable))
  sfit <- summary(fit)
  sfit
}, .parallel=T)
names(sfits2) <- rnames

df.sfit <- ldply(rnames, function(rname) {
  ts <- sfits2[[rname]]$tTable[,4]
  zs <- qt(sfits2[[rname]]$tTable[,5], Inf, lower.tail=F)
  sdf <- data.frame(roi=rname, measure=rownames(sfits2[[rname]]$tTable), tval=ts, pval=sfits2[[rname]]$tTable[,5], zval=zs)
  sdf[-1,] # remove intercept
})

library(fdrtool)
ret <- fdrtool(c(0.05,df.sfit$pval), statistic="pvalue")
head(ret$pval < 0.05)*1
head(ret$qval < 0.05)*1
head(ret$lfdr < 0.05)*1
tmp <- p.adjust(c(0.05,df.sfit$pval), method="fdr")

plot.ts(subset(df.sfit, measure=="avg_dist")$tval[1:10])
plot.ts(subset(df.sfit, measure=="avg_dist")$tval[-c(1:10)])

plot.ts(subset(df.sfit, measure=="knn_dist")$tval[1:10])
plot.ts(subset(df.sfit, measure=="knn_dist")$tval[-c(1:10)])

plot.ts(subset(df.sfit, measure=="frame_pose")$tval[1:10])
plot.ts(subset(df.sfit, measure=="frame_pose")$tval[-c(1:10)])

plot.ts(abs(subset(df.sfit, measure=="faces")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="pose")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="avg_dist")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="knn_dist")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="frame_pose")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="frame_diff")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="memorable")$tval[1:10]))
plot.ts(abs(subset(df.sfit, measure=="memorable")$tval[-c(1:10)]))

head(df.sfit, 10)

df.sfit2 <- subset(df.sfit, measure %in% c("faces", "pose", "avg_dist", "knn_dist", "frame_pose", "frame_diff"))
df.sfit2$hemi  <- sub("[.].*", "", df.sfit2$roi)
df.sfit2$rname <- sub("[lr]h.", "", df.sfit2$roi)

library(ggplot2)
library(ggthemr)
library(RColorBrewer)
ggthemr('pale', type='outer', text_size=14, layout='plain')

ggplot(subset(df.sfit2, hemi=="rh"), aes(x=rname, y=abs(tval), fill=rname)) + 
  geom_bar(stat="identity", position="dodge") + 
  #scale_fill_manual(values=rev(brewer.pal(7, "Spectral")[1:3])) + 
  facet_wrap(~measure, scales = "free_y") + 
  ylab("Percent Accuracy") + 
  scale_y_continuous(breaks=95:100, expand=c(0,0)) + 
  expand_limits(y=100.5) + 
  coord_cartesian(ylim=c(95,100.5)) + 
  theme(axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.line.x = element_line(color="black", size=1), 
        legend.position = "none")




# Let's make a function that runs 3dDeconvolve given some timing information
3dDeconvolve -x1D_stop -polort A -nodata 300 2 -x1D stdout: | 1dplot -one -stdin

list.files("/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/sub01/")

tfile <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/sub01/stimam_nnet2_mean_diff.txt"
cmd <- sprintf("3dDeconvolve -global_times -nodata 5088 1 -polort -1 -x1D_stop -num_stimts 1 -stim_times_AM1 1 %s 'SPMG1(2)' -x1D stdout:", tfile)
tmp <- system(cmd, intern=T)
tmp <- tmp[!grepl("^#", tmp)]
tmp <- tmp[tmp!=""]
tmp <- as.numeric(tmp)
length(tmp)
tmp2 <- df.xmats$frame_pose[1:5088]
cor(tmp, tmp2)

tfile1 <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/sub01/stim_faces.txt"
tfile2 <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/sub01/stim_questions.txt"
infiles <- "/data1/famface01/analysis/preprocessed/sub01/func/unfam_vids/rois_asap_ventral_peaks/asap_ventral_peaks_run*.1D"
infiles <- "/data1/famface01/analysis/preprocessed/sub01/func/unfam_vids/rois_asap_ventral_peaks/asap_ventral_peaks_run*.nii.gz"
cmd <- sprintf("3dDeconvolve -global_times -input '%s' -force_TR 1 -polort -1 -x1D_stop -num_stimts 3 -stim_times 1 %s 'SPMG1(2)' -stim_times 2 %s 'SPMG1(4)' -stim_times_AM1 3 %s 'SPMG1(2)' -x1D tmp.txt", infiles, tfile1, tfile2, tfile)
system(cmd)
tmp3 <- read.xmat("tmp.txt.xmat.1D")
head(tmp3)
cor(tmp3[,ncol(tmp3)], tmp)
cor(tmp3[,ncol(tmp3)], tmp2)
cor(tmp3[,2], df.xmats$faces[1:5088])
cor(xmat$avg_dist, tmp)

cor()

xfile <- sprintf("%s/%s/face_deviations_unfam/nnet2_only_avgdist2.reml/xmat.1D", taskdir, subj)
xmat <- read.xmat(xfile, rm.nums=T)



subjects <- sprintf("sub%02i", 1:6)
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








