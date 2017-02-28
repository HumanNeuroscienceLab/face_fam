

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


# Load Trait Data -----------------------------------------------------------

# Load demographics information
base <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
demo.vnames <- demos$video # typo
demo.vnames <- sub("_fr[0-9]{3}", "", demo.vnames)
df.demos    <- demos[,-c(1:2)]
df.demos    <- df.demos[,-c(6:7)] # remove hair and eye color

# Load trait information
base         <- "/data1/famface01/analysis/encoding/12_Features"
traits       <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
trait.vnames <- as.character(traits[,1])
traits       <- traits[,-1]
# Reorder to match with demos
oinds       <- sapply(demo.vnames, function(x) which(trait.vnames == x))
all.equal(trait.vnames[oinds], demo.vnames)
df.traits   <- traits[oinds,]

# Save vnames
vnames <- demo.vnames

# Relevant functions
formula_to_mat <- function(formula, data) {
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  return(X)
}

compile.demo.trait <- function(df.demos, df.traits, ret.intercept=TRUE) {
  ## DEMOGRAPHICS
  # Remove the effect of gender from makeup
  df.demos2 <- df.demos
  df.demos2$makeup <- lm(makeup ~ gender, data=df.demos2)$residuals
  # Center age
  df.demos2$age <- scale(df.demos2$age, scale=F, center=T)
  # Set the reference for gender to Male
  df.demos2 <- within(df.demos2, gender <- relevel(gender, ref = "Male"))
  # Set the reference for facial hair to when none
  df.demos2 <- within(df.demos2, facial_hair <- relevel(facial_hair, ref = "None"))
  # Set the reference for race to when white
  df.demos2 <- within(df.demos2, race <- relevel(race, ref = "White"))
  # Get the matrix
  mat.demos1 <- formula_to_mat(~.-1, df.demos2)
  mat.demos2 <- formula_to_mat(~., df.demos2) # note this will have the intercept
  
  ## TRAITS
  mat.traits <- scale(df.traits, scale=F, center=T)
  
  if (ret.intercept) {
    return(cbind(mat.demos2, mat.traits))
  } else {
    return(cbind(mat.demos1, mat.traits))
  }
}

# Compile the demographics and traits
mat.all <- compile.demo.trait(df.demos, df.traits)


# Load Shape/Texture Data ---------------------------------------------------

# Finally load in the shape/texture data (this is )
basedir <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes) # vnames to use
sym.shapes.eigs <- read.table(file.path(basedir, "shape_sym_pca_eigs.txt"))
sym.textures    <- read.csv(file.path(basedir, "texture_sym_pca_scores.csv"), row.names = 1)
## combine
X1 <- sym.shapes[,1:65]
X2 <- sym.textures[,1:200]
colnames(X1) <- sprintf("shape_%s", colnames(X1))
colnames(X2) <- sprintf("texture_%s", colnames(X2))
X <- cbind(X1,X2)


# Rearrange rows for shape/texture ------------------------------------------

# rearrange of the shapes here
## get the new inds
inds <- sapply(vnames, function(vname) which(shape.vnames==vname))
all.equal(shape.vnames[inds], vnames)
## set new inds
X0 <- X # save
X <- as.matrix(X[inds,])


# OLD: Reduce Dimensionality -----------------------------------------------------

### I did this before but skip to next section for where I'm at now

# So from this analysis, it seems that maybe should split up traits from other things

# Use factor analysis on demographic/traits
library(psych)

# remove intercept
mat.all2 <- mat.all[,-1]
#mat.all2 <- mat.all2[,c(1,9,15,16:25)]
head(mat.all2)

# Determine the number of components needed
fa.parallel(mat.all2) # suggests 9 comps
vss(mat.all2, 25) # suggests 18 or 21 (going with 21)
# also went with 21 over 18 because it includes middle eastern participants
# also there are local minima at 16 and 11

# Get the factors
fac.res <- psych::fa(mat.all2, nfactors=18, residuals=T, rotate='varimax', 
                     fm='minres')
#loadings(fac.res)
#residuals(fac.res)

# Save the loadings and scores, rename the columns
fac.loads  <- unclass(fac.res$loadings)
fac.scores <- fac.res$scores
colnames(fac.loads) <- sprintf("factor%02i", 1:ncol(fac.loads))
colnames(fac.scores) <- sprintf("factor%02i", 1:ncol(fac.scores))

# Plot
library(corrplot)
hc1     <- hclust(dist(fac.loads))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- fac.loads[ord.hc1,]
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)


# Select Trait Factors ----------------------------------------------------

mat.all2 <- mat.all[,-1]
#mat.all2 <- mat.all2[,c(1,9,16:25)]
mat.all2 <- mat.all2[,c(16:25)]
head(mat.all2)

# Determine the number of components needed
fa.parallel(mat.all2) # suggests 5 factors
vss(mat.all2, ncol(mat.all2)) # also suggests 5, 6, 7 or 9 factors

# Get the factors (VSS complexity 2 maximizes at 6)
fac.res <- psych::fa(mat.all2, nfactors=6, residuals=T, rotate='varimax', 
                     fm='minres')

# Save the loadings and scores, rename the columns
fac.loads  <- unclass(fac.res$loadings)
fac.scores <- fac.res$scores
colnames(fac.loads) <- sprintf("factor%02i", 1:ncol(fac.loads))
colnames(fac.scores) <- sprintf("factor%02i", 1:ncol(fac.scores))

# Plot
library(corrplot)
hc1     <- hclust(dist(fac.loads))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- fac.loads[ord.hc1,]
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

# Assign Names
fac.names <- c("unemotional", "competent", "trustworthy", "typical", 
               "memorable", "attractive")
colnames(fac.scores) <- fac.names


# Other Features ----------------------------------------------------------

mat.all3 <- mat.all[,-1]
mat.all3 <- mat.all3[,-c(16:25)]

# split up
facial.hair <- mat.all3[,c(2:7)]
race <- mat.all3[,10:14]
misc <- mat.all3[,c(1,8,9,15)]

# Combine it all back
all.mat <- cbind(fac.scores, misc, facial.hair, race)
all.groups <- rep(c("traits", "age", "makeup", "gender", "glasses", "facial_hair", "race"), c(ncol(fac.scores), 1, 1, 1, 1, ncol(facial.hair), ncol(race)))



# Classification Functions --------------------------------------------

source("/data1/famface01/command/misc/face_representations/misc/cv_glmnet.R")


# Classification ----------------------------------------------------------




# Predict the scores
# for each score, we repeat the process 10 times
library(cvTools)
nreps <- 10
cfolds <- cvFolds(nrow(fac.scores), K=10, R=nreps)
cvrep.fits <- llply(1:ncol(fac.scores), function(i) {
  y <- fac.scores[,i]
  lapply(1:nreps, function(i) {
    foldid <- cfolds$which[cfolds$subsets[,i]]
    run_cvglmnet(X, y, foldid=foldid, parallel=F, exclude.zero=T)
  })
}, .parallel=T)
names(cvrep.fits) <- colnames(fac.scores)

# Calculate the average model fits across the 10 repeats
cvrep.r2s <- sapply(cvrep.fits, function(fits) {
  sapply(fits, function(x) x$bestfit$val)
})
fac.r2s <- colMeans(cvrep.r2s)
round(fac.r2s, 3)

# Combine the predicted values
fac.preds <- sapply(cvrep.fits, function(fits) {
  predvals <- sapply(fits, function(x) x$bestfit$preval)
  predvals <- rowMeans(predvals)
  return(predvals)
})

# Get the residuals
fac.resids <- sapply(1:ncol(fac.preds), function(i) {
  lm(fac.scores[,i] ~ fac.preds[,i])$residuals
})

# Combine the r2 values across the lambdas
mat.cvr2s <- lapply(cvrep.fits, function(fits) {
  lambda <- fits[[1]]$lambda
  cvr2   <- matrix(fits[[1]]$measures$rsq, nrow=1)
  cvnz   <- matrix(fits[[1]]$nzero, nrow=1)
  for (i in 2:nreps) {
    lambda0 <- fits[[i]]$lambda
    rsq0    <- fits[[i]]$measures$rsq
    nzeros0 <- fits[[i]]$nzero
    
    newlambda <- intersect(lambda, lambda0)
    cvr2      <- rbind(cvr2[,is.element(lambda, newlambda)], 
                       rsq0[is.element(lambda0, newlambda)])
    cvnz      <- rbind(cvnz[,is.element(lambda, newlambda)], 
                       nzeros0[is.element(lambda0, newlambda)])
    lambda    <- newlambda
  }
  # average across the repeats
  cbind(lambda=lambda, rsq=colMeans(cvr2), nzeros=colMeans(cvnz))
})
maxlambdas <- sapply(mat.cvr2s, function(x) {
  max.rsq <- max(x[x[,3] > 0, 2]) # not necessary as all are non-zero
  wmax.rsq<- x[,2] == max.rsq
  x[wmax.rsq,1]
})
## double check that not getting zero coefficients
sapply(mat.cvr2s, function(x) x[which.max(x[,2]),3])

# Determine the coefficients from the full model
# using the lambda with the best r2
## get full model fit first
full.fits <- lapply(1:ncol(fac.scores), function(i) {
  glmnet(X, fac.scores[,i], family="gaussian", lambda=maxlambdas[i])
})
names(full.fits) <- colnames(fac.scores)
## compile the coefficients
fac.coefs <- sapply(full.fits, function(x) as.vector(x$beta))
## note: two of the factors (9 and 16) have all non-zero values
round(colMeans(fac.coefs!=0), 2)

# # of coefs relative to # of r2
plot(round( cbind(colMeans(fac.coefs!=0), fac.r2s), 3 ))

# So we finally need to save this somehow, right
pca.face.feats <- X
outdir <- "/data1/famface01/analysis/misc/320_roi_task_activity"
dir.create(outdir)
save(fac.loads, pca.face.feats, fac.r2s, fac.preds, fac.resids, fac.coefs, fac.scores, demo.vnames, 
     file=file.path(outdir, "20_predict_face_feats.rda"))
