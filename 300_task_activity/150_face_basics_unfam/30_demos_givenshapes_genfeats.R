# This script saves the classification fits for age and gender to be used as 
# voxelwise regressors in a regression analysis


###
# SETUP
###

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(plyr)

load("/data1/famface01/command/misc/stimuli/Classify/tmp_lst_fits.rda")


###
# Video Info
###

basedir <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes) # vnames to use


###
# Gender
# - we get the probability of each gender not the actual gender
###

# Get the data to fit
newx <- as.matrix(lst.fits$shape$gender$trainingData[,-ncol(lst.fits$shape$gender$trainingData)])
# Get the lambda
lam <- lst.fits$shape$gender$bestTune$lambda
probs <- predict(lst.fits$shape$gender$finalModel, newx, s=lam, type="response")
# visualize by gender
newy <- lst.fits$shape$gender$trainingData[,ncol(lst.fits$shape$gender$trainingData)]
plot.ts(probs[order(newy)])
# for saving
df.gender <- data.frame(shape.gender=probs[,1], gender=newy, row.names=shape.vnames)


###
# Age
# - this is the predicted value
###

pmat    <- lst.fits$shape$age$pred

# we will need to average by the repititions
reps    <- as.numeric(sub("Fold[0-9]{2}[.]Rep", "", pmat$Resample))

# order by row indices
rowInds <- pmat$rowIndex

# Get the predictions and average across reps
preds   <- pmat$pred[order(reps, rowInds)]
preds   <- (preds[sort(reps)==1] + preds[sort(reps)==2])/2

# Get the observed values and average across reps
obs     <- pmat$obs[order(reps, rowInds)]
all(obs[sort(reps)==1] == obs[sort(reps)==2])
obs     <- (obs[sort(reps)==1] + obs[sort(reps)==2])/2

# for saving
df.age  <- data.frame(shape.age=preds, age=obs, row.names=shape.vnames)


###
# Save
###

df <- cbind(df.gender, df.age)
write.csv(df, file="/data1/famface01/command/misc/face_representations/300_task_activity/150_face_basics_unfam/z_demos_givenshape.csv")
