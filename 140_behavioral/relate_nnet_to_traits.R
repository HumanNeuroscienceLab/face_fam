

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

# Load neural network measure
load.nn.openface <- function(vnames=NULL) {
  base <- "/data1/famface01/analysis/encoding/12_Features"
  
  # Read in
  labels <- read.csv(sprintf('%s/labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$vid <- basename(as.character(labels.df$img))
  labels.df$vid <- sub(".png", "", labels.df$vid)
  labels.df$vid <- sub("_fr[0-9]{3}$", "", labels.df$vid)
  labels.df$vid <- factor(labels.df$vid)
  
  # Add the frame index
  frs <- sapply(1:nrow(labels.df), function(i) sub(as.character(labels.df$vid)[i], "", basename(as.character(labels.df$img)[i])))
  frs <- sub("_fr", "", frs)
  frs <- sub(".png", "", frs)
  frs <- as.numeric(frs)
  labels.df$frame <- frs
  
  # Remove vids of famous ppl
  inds <- sapply(as.character(labels.df$vid), function(vid) {
    any(vnames == vid)
  })
  labels.df <- labels.df[inds,]
  features  <- features[inds,]
  labels.df$vid <- factor(labels.df$vid)
  
  # Add index
  labels.df$X <- 1:nrow(labels.df)
  
  # reorder frames
  labels.df2 <- ddply(labels.df, .(vid), function(x) {
    x[order(x$frame),]
  })
  
  # Match up
  if (!is.null(vnames)) {
    oinds <- unlist(lapply(vnames, function(x) which(labels.df2$vid == x)))
    all.equal(as.character(labels.df2$vid[oinds[seq(1,length(oinds),by=8)]]), vnames)
    df.labels <- labels.df2[oinds,]
    mat.feats <- features[df.labels$X,]
  } else {
    df.labels <- labels.df2
    mat.feats <- features
  }
  
  return(list(labs=df.labels, feats=mat.feats))
}



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

# Load the neural network measure
openface <- load.nn.openface(vnames)
all.equal(as.character(openface$labs$vid[openface$labs$frame==3]), vnames)
openface$labs$X2 <- 1:nrow(openface$labs)
mat.nn.vid <- daply(openface$labs, .(vid), function(labs) {
  colMeans(openface$feats[labs$X2,])
}, .progress="text")
mat.nn.framediff <- daply(openface$labs, .(vid), function(labs) {
  openface$feats[labs$X2,]
  colMeans()
}, .progress="text")



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
  run_repeated_cvglmnet(mat.nn.vid, df[,i], family="gaussian", type.measure="rsq", 
                        parallel=F, k=10, nreps=10)
}, .parallel=T)
names(cvrep.gaussian) <- colnames(df)[c(1:7,9)]
rsqs <- sapply(cvrep.gaussian, function(x) x$mean.max.res)
rsqs


#
# Multinomial:
# Facial hair, race, hair
# 
## merge Goatee / Goatee and moustache / Moustache
df2$facial_hair <- revalue(df$facial_hair, c("Moustache"="Goatee or moustache", "Goatee"="Goatee or moustache", "Goatee and moustache"="Goatee or moustache"))
cvrep.fh <- run_repeated_cvglmnet(mat.nn.vid, df2$facial_hair, family="multinomial", type.measure="Kappa", 
                                  parallel=T, k=10, nreps=10)
caret::confusionMatrix(cvrep.fh$fitted.class, df2$facial_hair)

## race
cvrep.race <- run_repeated_cvglmnet(mat.nn.vid, df2$race, family="multinomial", type.measure="Kappa", 
                                    parallel=T, k=5, nreps=10)

## eye (not significant)
df2$eye <- revalue(df$eye, c("Amber"="Hazel"))
cvrep.eye <- run_repeated_cvglmnet(mat.nn.vid, df2$eye, family="multinomial", type.measure="Kappa", 
                                   parallel=T, k=10, nreps=10)
caret::confusionMatrix(cvrep.eye$fitted.class, df2$eye)

## (for fun) hair
cvrep.hair <- run_repeated_cvglmnet(mat.nn.vid, df2$hair, family="multinomial", type.measure="Kappa", 
                                    parallel=T, k=10, nreps=10)
caret::confusionMatrix(cvrep.hair$fitted.class, df2$hair)

cvrep.multinomial <- list(facial_hair=cvrep.fh, race=cvrep.race, eye=cvrep.eye, hair=cvrep.hair)

df.acc <- ldply(names(cvrep.multinomial), function(name) {
  confmat <- caret::confusionMatrix(cvrep.multinomial[[name]]$fitted.class, df2[[name]])
  c(name=name, confmat$overall)
})
df.acc1 <- df.acc
print(df.acc)


#
# Binomial:
# gender and glasses
#
cvrep.gender <- run_repeated_cvglmnet(mat.nn.vid, df2$gender, family="binomial", type.measure="Kappa", 
                                      parallel=T, k=10, nreps=10)

cvrep.glasses <- run_repeated_cvglmnet(mat.nn.vid, df2$glasses, family="binomial", type.measure="Kappa", 
                                       parallel=T, k=10, nreps=10)

cvrep.binomial <- list(gender=cvrep.gender, glasses=cvrep.glasses)

df.acc <- ldply(names(cvrep.binomial), function(name) {
  confmat <- caret::confusionMatrix(cvrep.binomial[[name]]$fitted.class, df2[[name]])
  c(name=name, confmat$overall)
})
df.acc2 <- df.acc
print(df.acc)


###
# So got the classifer results. Save the predicted results.
# For the multinomial and binomial, save both the classes and probabilities
#

save(vnames, df2, cvrep.gaussian, cvrep.multinomial, cvrep.binomial, file="/data1/famface01/analysis/misc/320_roi_task_activity/40_predict_face_feats_df.rda")

