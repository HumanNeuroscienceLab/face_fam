# This script saves the classification fits for age and gender to be used as 
# voxelwise regressors in a regression analysis


###
# SETUP
###

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(plyr)

load("/data1/famface01/command/misc/face_representations/120_features/tmp_combined_fits.rda")


###
# Video Info
###

basedir <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes) # vnames to use


###
# Traits - Factor Analysis
# - we get the probability of each gender not the actual gender
###

lst.dfs <- llply(1:length(fits$trait_factors), function(i) {
  pmat    <- fits$trait_factors[[i]]$pred
  
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
  obs     <- rowMeans(sapply(ureps, function(ri) obs[sort(reps)==ri]))
  
  # Get the residuals of the obs
  resid   <- lm(obs ~ preds)$residuals
  
  # for saving
  df.traitx  <- data.frame(preds, resid, obs, row.names=shape.vnames)
  colnames(df.traitx) <- c(sprintf("shape.trait%i", i), 
                           sprintf("resid.trait%i", i), 
                           sprintf("trait%i", i))
  
  return(df.traitx)
})
df.traitx <- do.call(cbind, lst.dfs)


###
# Save
###

write.csv(df.traitx, file="/data1/famface01/command/misc/face_representations/300_task_activity/150_face_basics_unfam/z_traitsfa_givenshape.csv")
