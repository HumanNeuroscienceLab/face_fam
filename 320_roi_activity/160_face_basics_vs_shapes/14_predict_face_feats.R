

# This script puts together the 18 or so factor scores collapsing the trait, 
# demographic and physical information. It then uses the PCA face features to 
# predict each of those factor scores


# Setup -------------------------------------------------------------------

# Setup
if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(doMC))
registerDoMC(30)
library(RColorBrewer)
library(dynamicTreeCut)

library(ggplot2)
source("/data1/famface01/command/encoding/FirstPass/plot_functions.R")
source("/data1/famface01/command/encoding/SubRois_Unfam/lm_functions.R")

library(glmnet)


# Load Data -----------------------------------------------------------

# Features
df.demos <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_df-demos.csv")[,-c(1:2)]
vnames <- read.table("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_demo-vnames.txt")
vnames <- as.character(vnames[,1])
all.mat <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-mat.csv")[,-1]

# Face Features
X <- read.csv("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_sym-shape-texture.csv")
all.equal(vnames, as.character(X$X)) # video order equal
X <- X[,-1]
X <- as.matrix(X)

# Center/scale all.mat columns
all.grps <- as.character(read.table("/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-groups.txt")[,1])
all.cmat <- all.mat
all.cmat[,all.grps %in% c("traits","age","makeup")] <- scale(all.cmat[,all.grps %in% c("traits","age","makeup")], center=T, scale=F)
all.smat <- all.mat
all.smat[,all.grps %in% c("traits","age","makeup")] <- scale(all.cmat[,all.grps %in% c("traits","age","makeup")], center=T, scale=T)

# Final data-frame with traits
traits <- as.matrix(all.mat[,all.grps=="traits"])
df <- cbind(traits, df.demos)
## scale
dfs <- df
dfs[,c(1:7,9)] <- scale(dfs[,c(1:7,9)])



# Classification Functions --------------------------------------------

source("/data1/famface01/command/misc/face_representations/misc/cv_glmnet.R")


# Classification: DF ----------------------------------------------------------

# We get the classification results from the data-frame values
# This means we will sometimes be running a multinomial classifier

head(df)
df2 <- df

#
# Gaussian:
# Run for all the trait factor scores and age + makeup
# 
cvrep.gaussian <- llply(c(1:7,9), function(i) {
  run_repeated_cvglmnet(X, df[,i], family="gaussian", type.measure="rsq", 
                        parallel=F, k=10, nreps=10)
}, .parallel=T)
names(cvrep.gaussian) <- colnames(df)[c(1:7,9)]


#
# Multinomial:
# Facial hair, race, hair
# 
## merge Goatee / Goatee and moustache / Moustache
df2$facial_hair <- revalue(df$facial_hair, c("Moustache"="Goatee or moustache", "Goatee"="Goatee or moustache", "Goatee and moustache"="Goatee or moustache"))
cvrep.fh <- run_repeated_cvglmnet(X, df2$facial_hair, family="multinomial", type.measure="Kappa", 
                                  parallel=T, k=10, nreps=10)
caret::confusionMatrix(cvrep.fh$fitted.class, df2$facial_hair)

## race
cvrep.race <- run_repeated_cvglmnet(X, df2$race, family="multinomial", type.measure="Kappa", 
                                    parallel=T, k=5, nreps=10)

## eye (not significant)
df2$eye <- revalue(df$eye, c("Amber"="Hazel"))
cvrep.eye <- run_repeated_cvglmnet(X, df2$eye, family="multinomial", type.measure="Kappa", 
                                   parallel=T, k=10, nreps=10)
caret::confusionMatrix(cvrep.eye$fitted.class, df2$eye)

## (for fun) hair
cvrep.hair <- run_repeated_cvglmnet(X, df2$hair, family="multinomial", type.measure="Kappa", 
                                    parallel=T, k=10, nreps=10)
caret::confusionMatrix(cvrep.hair$fitted.class, df2$hair)

cvrep.multinomial <- list(facial_hair=cvrep.fh, race=cvrep.race, eye=cvrep.eye, hair=cvrep.hair)


#
# Binomial:
# gender and glasses
#
cvrep.gender <- run_repeated_cvglmnet(X, df2$gender, family="binomial", type.measure="Kappa", 
                                      parallel=T, k=10, nreps=10)

cvrep.glasses <- run_repeated_cvglmnet(X, df2$glasses, family="binomial", type.measure="Kappa", 
                                       parallel=T, k=10, nreps=10)

cvrep.binomial <- list(gender=cvrep.gender, glasses=cvrep.glasses)


###
# So got the classifer results. Save the predicted results.
# For the multinomial and binomial, save both the classes and probabilities
#

save(vnames, df2, cvrep.gaussian, cvrep.multinomial, cvrep.binomial, file="/data1/famface01/analysis/misc/320_roi_task_activity/40_predict_face_feats_df.rda")




# Classify: Caret --------------------------------------------------------


library(caret)
fitControl <- trainControl(
  method = "repeatedcv",
  number = 10, 
  repeats = 10, 
  allowParallel = TRUE, 
  savePredictions = "final", 
  returnResamp = "final"
)

i <- 1
lambdas <- glmnet(X, df[,i], family="gaussian", nlambda=100, alpha=1)$lambda
grids <- expand.grid(alpha=1, lambda=lambdas)

fit <- train(X, df[,i], 
             method = "glmnet",
             trControl = fitControl, 
             tuneGrid = grids)

perf  <- getTrainPerf(fit)
tune  <- fit$bestTune
vi    <- varImp(fit, scale=F)

pmat <- fit$pred

# we will need to average by the repititions
reps    <- as.numeric(sub("Fold[0-9]{2}[.]Rep", "", pmat$Resample))
ureps   <- sort(unique(reps))

# order by row indices
rowInds <- pmat$rowIndex

# Get the predictions and average across reps
preds   <- pmat$pred[order(reps, rowInds)]
#reps    <- reps[order(reps, rowInds)]
preds   <- rowMeans(sapply(ureps, function(ri) preds[sort(reps)==ri]))

# Get the observed values and average across reps
obs     <- pmat$obs[order(reps, rowInds)]
all(obs[sort(reps)==1] == obs[sort(reps)==2])
all(obs[sort(reps)==1] == df[,1])
obs     <- rowMeans(sapply(ureps, function(ri) obs[sort(reps)==ri]))
all(obs == df[,1])

# Get the residuals of the obs
resid   <- lm(obs ~ preds)$residuals

cor(preds, cvrep.gaussian[[1]]$fitted) # very similar
tmp <- abs(cvrep.gaussian[[1]]$resids - resid)
mean(tmp)




# Classify: Mat -----------------------------------------------------------

## DONT DO IT

# Here we do each on their own
y <- factor(all.mat[,11], levels=0:1, labels=c("no","yes"))
cvrep2.fh1 <- run_repeated_cvglmnet(X, y, family="binomial", type.measure="Kappa", 
                                    parallel=F, k=10, nreps=10)
confusionMatrix(cvrep2.fh1$fitted.class, y)




# So we finally need to save this somehow, right
pca.face.feats <- X
outdir <- "/data1/famface01/analysis/misc/320_roi_task_activity"
dir.create(outdir)
save(fac.loads, pca.face.feats, fac.r2s, fac.preds, fac.resids, fac.coefs, fac.scores, demo.vnames, 
     file=file.path(outdir, "20_predict_face_feats.rda"))

