
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

# Beta Data
load("/data1/famface01/analysis/misc/320_roi_task_activity/dgamma_betas.rda", verbose=T)

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

# Rnames
rnames2 <- sub("[.]", " ", rnames)
rnames2 <- sub("^r", "R", rnames2)
rnames2 <- sub("^l", "L", rnames2)

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

all.mat2 <- all.mat[inds,]
all.cmat2 <- all.cmat[inds,]
all.smat2 <- all.smat[inds,]



# Regression --------------------------------------------------------------

# are those things o

###
# Regression
#

## fits
fits <- lm(grp.bs1[,,2] ~ ., data=all.mat2)
sfits <- summary(fits)
lm.tstats <- sapply(sfits, function(sfit) sfit$coefficients[-1,3])
colnames(lm.tstats) <- rnames
## heatmap
# colors
#cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
cols <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
all.grps2 <- all.grps
ccols <- factor(all.grps2, levels=levels(factor(all.grps2)), 
                labels=rev(brewer.pal(length(unique(all.grps2)), "Set3")))
ccols <- as.character(ccols)
# variables
mat <- t(lm.tstats)
tmat <- mat * (abs(mat) > 1.96)
zlim <- max(tmat[tmat!=0]) * c(-1,1)
tmat[tmat==0] <- NA
# separator lines
seplines <- c(0,cumsum(rle(all.grps2)$lengths)) + 0.5
# plot
heatmap(tmat, scale="none", col=cols, ColSideColors = ccols, 
        zlim=zlim, Rowv=NA, Colv=NA, 
        margins=c(9,5), cexCol = 1.25, cexRow = 1.25, 
        add.expr = abline(v=seplines, h=6.5, col="grey50"))


###
# ANOVA
# 

fits <- aov(grp.bs1[,,2] ~ ., data=df)
sfits <- summary(fits)
aov.tstats <- sapply(sfits, function(sfit) {
  fstats <- sfit$`F value`
  names(fstats) <- rownames(sfit)
  tstats <- sqrt(fstats[-length(fstats)])
  tstats
})
colnames(aov.tstats) <- rnames

# colors
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
all.grps2 <- all.grps
# variables
mat <- t(aov.tstats)
tmat <- mat * (abs(mat) > 1.65)
zlim <- max(tmat[tmat!=0]) * c(-1,1)
tmat[tmat==0] <- NA
# plot
heatmap(tmat, scale="none", col=cols, 
        zlim=zlim, Rowv=NA, Colv=NA, 
        margins=c(9,5), cexCol = 1.25, cexRow = 1.25)


# Forward Classify ----------------------------------------------------------

source("/data1/famface01/command/misc/face_representations/misc/cv_glmnet.R")


# Why does forward not work that well?

forward.mods0 <- llply(1:length(rnames), function(ri) {
  y <- sub.bs1[,,ri,2]
  dim(y) <- prod(dim(y))
  y <- as.vector(y)
  X <- as.matrix(all.mat)[rep(1:nrow(all.mat),6),]
  foldid <- rep(1:length(subjects), each=nrow(all.mat))
  ret1 <- run_cvglmnet(X=X, y=y, foldid=foldid, family="gaussian", type.measure="rsq", parallel=F, exclude.zero = T)
  ret1
}, .parallel=T)
names(forward.mods0) <- rnames
sapply(forward.mods0, function(x) x$bestfit$val)
forward.coefs0 <- sapply(forward.mods0, function(x) x$bestfit$coef[-1,])


forward.mods <- llply(1:length(rnames), function(ri) {
  ret1 <- run_repeated_cvglmnet(X=as.matrix(all.mat), y=grp.bs1[,ri,2], 
                                family="gaussian", type.measure="rsq", parallel=F, ezero = T)
  ret1
}, .parallel=T)
names(forward.mods) <- rnames
sort(sapply(forward.mods, function(x) x$mean.max.res))
forward.mods$l.EBA$coefs

forward.coefs <- sapply(forward.mods, function(x) x$coefs)

# Heatmap
## colors
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
all.grps2 <- all.grps
ccols <- factor(all.grps2, levels=levels(factor(all.grps2)), 
                labels=rev(brewer.pal(length(unique(all.grps2)), "Set3")))
ccols <- as.character(ccols)
## values
cmat <- abs(t(forward.coefs))
zlim <- range(abs(cmat[cmat!=0]))
cmat[cmat==0] <- NA
# separator lines
seplines <- c(0,cumsum(rle(all.grps2)$lengths)) + 0.5
# plot
heatmap(cmat, scale="none", col=cols, ColSideColors = ccols, 
        zlim=zlim, Rowv=NA, Colv=NA, 
        margins=c(9,5), cexCol = 1.25, cexRow = 1.25, 
        add.expr = abline(v=seplines, h=6.5, col="grey50"))

# get the lasso score sig
library(lassoscore)
forward.sigs <- laply(1:length(rnames), function(ri) {
  tmp <- run_repeated_cvglmnet(X=as.matrix(all.mat), y=grp.bs1[,ri,2], 
                                family="gaussian", type.measure="rsq", 
                                parallel=F, ezero=T)
  tmp2 <- lassoscore(grp.bs1[,ri,2], as.matrix(all.mat), lambda=tmp$lambda)
  sigs <- tmp2$p.model
  names(sigs) <- colnames(all.mat)
  sigs
}, .parallel=T)
forward.sigs # haha nothing is significant
lasso.sigz <- qt(forward.sigs, Inf, lower.tail=F)
rownames(lasso.sigz) <- rnames

cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
all.grps2 <- all.grps
ccols <- factor(all.grps2, levels=levels(factor(all.grps2)), 
                labels=rev(brewer.pal(length(unique(all.grps2)), "Set3")))
ccols <- as.character(ccols)


# HEATMAP
mat <- lasso.sigz
tmat <- mat * (mat > 1.96)
zlim <- range(tmat[tmat!=0])
tmat[tmat==0] <- NA
# separator lines
seplines <- c(0,cumsum(rle(all.grps2)$lengths)) + 0.5
# heatmap
heatmap(tmat, scale="none", col=cols, ColSideColors = ccols, 
        zlim=zlim, Rowv=NA, Colv=NA, 
        margins=c(9,5), cexCol = 1.25, cexRow = 1.25, 
        add.expr = abline(v=seplines, h=6.5, col="grey50"))
## so this isn't looking all that good



# Inverse Classify ----------------------------------------------------------

###
# Here I calculate the inverse models. 
# - I skip the amber eyes feature because it has too few observations with those eyes
# - I split the multinomial results into many binomial ones

all.mat2 <- all.mat[,-21] # remove amber eyes
inverse.mods <- llply(1:ncol(all.mat2), function(i) {
  if (all(all.mat2[,i] %in% c(0,1))) {
    family="binomial"
    tm="Kappa"
    y <- as.matrix(all.mat2)[,i]
    y <- factor(y, levels=c(0,1), labels=c("no","yes"))
  } else {
    family="gaussian"
    tm="rsq"
    y <- as.matrix(all.mat2)[,i]
  }
  ret1 <- run_repeated_cvglmnet(X=grp.bs1[,,2], y=y, 
                                family=family, type.measure=tm, parallel=F, 
                                ezero=T, alpha=0.5)
  ret1
}, .parallel=T)
names(inverse.mods) <- colnames(all.mat2)

#
# I can plot the mean R-squared or Kappa values
# 
mar <- par()$mar; par(mar=c(10.1, 4.1 ,4.1 ,2.1))
barplot(sort(sapply(inverse.mods, function(x) x$mean.max.res)), las=2)
par(mar=mar)

#
# Plot the coefficients with a heatmap
#
inverse.coefs <- sapply(inverse.mods, function(x) x$coefs)
## colors
cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
all.grps2 <- all.grps[-21]
ccols <- factor(all.grps2, levels=levels(factor(all.grps2)), 
                labels=rev(brewer.pal(length(unique(all.grps2)), "Set3")))
ccols <- as.character(ccols)
## values
cmat <- abs(inverse.coefs)
zlim <- range(abs(cmat[cmat!=0]))
cmat[cmat==0] <- NA
## separator lines
seplines <- c(0,cumsum(rle(all.grps2)$lengths)) + 0.5
# plot
heatmap(cmat, scale="none", col=cols, ColSideColors = ccols, 
        zlim=zlim, Rowv=NA, Colv=NA, 
        margins=c(9,5), cexCol = 1.25, cexRow = 1.25, 
        add.expr = abline(v=seplines, h=6.5, col="grey50"))



# Play with Inverse -------------------------------------------------------

# Trying caret and stuff on the data with inverse model

library(caret)

fitControl <- trainControl(
  method = "repeatedcv", 
  number = 10, # 10-fold
  repeats = 2, # with 10 repeats (you can start with 1 or 2 for speed)
  classProbs = TRUE, # allow you to get class probabilities (only for category classification)
  allowParallel = T, # run in parallel
  returnResamp = "final", 
  savePredictions = "all"
)

X <- grp.bs1[,,2]
i <- 10
y <- as.matrix(all.mat2)[,i]
y <- factor(y, levels=c(0,1), labels=c("no","yes"))

# no works
gn <- glmnet(X, y, family="binomial", alpha=1)
fit <- train(X, y, # X = observations x featurs, y = observations
             method = "glmnet", # can also use svmLinear or other method
             trControl = fitControl, # from above
             preProcess = c("center", "scale"), # if you want to center and scale the X data
             tuneGrid = data.frame(alpha=1, lambda=gn$lambda), metric="Kappa")
#tuneLength = 10) # this gives 10 decent values for each of the tuning parameters, can also do this manually by specifying tuneGrid
print(fit) # gives you summary

# no works
tg = expand.grid(n.trees=c(50,100,200), interaction.depth=1:3, shrinkage=c(0.1, 0.2), n.minobsinnode=10)
fit <- train(X, y, # X = observations x featurs, y = observations
             method = "gbm", # can also use svmLinear or other method
             trControl = fitControl, # from above
             preProcess = c("center", "scale"), # if you want to center and scale the X data
             tuneGrid = tg, 
             verbose=F, metric="Kappa")
#tuneLength = 10) # this gives 10 decent values for each of the tuning parameters, can also do this manually by specifying tuneGrid
print(fit) # gives you summary


# Facial Hair (always none)
y <- dfs$facial_hair
y <- revalue(y, c("Goatee"="Goatee and moustache"))
ret1 <- run_cvglmnet(X=grp.bs1[,,2], y=y, 
                     family="multinomial", type.measure="Kappa", parallel=T, 
                     exclude.zero=T, alpha=1)
ret1$bestfit$confMat

# Race (all white)
y <- dfs$race
#y <- revalue(y, c("Goatee"="Goatee and moustache"))
ret1 <- run_cvglmnet(X=grp.bs1[,,2], y=y, 
                     family="multinomial", type.measure="Kappa", parallel=T, 
                     exclude.zero=T, alpha=1)




