# This is all to create an initial dataset to classify age and gender
# (see Classify in misc/stimuli/Classify for more code)
# 

# Setup -------------------------------------------------------------------

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

library(caret)
library(glmnet)
library(gplots)



# Functions ---------------------------------------------------------------

get_bestfit <- function(cvfit, type.measure) {
  bestouts <- cvfit$measures[[type.measure]]
  extreme  <- ifelse(type.measure == "rmse", min, max)
  which.extreme <- ifelse(type.measure == "rmse", which.min, which.max)
  
  val <- extreme(bestouts)
  ind <- which.extreme(bestouts)
  
  bestfit  <- list(
    measure = type.measure, 
    val     = val, 
    ind     = ind, 
    lam     = cvfit$lambda[ind], 
    preval  = cvfit$fit.preval[,ind], 
    nzero   = cvfit$nzero[ind], 
    coef    = coef(cvfit, s=cvfit$lambda[ind])
  )
  bestfit
}

run_cvglmnet <- function(X, y, keep=T, parallel=T, type.measure="rsq", ...) 
{
  if (!(type.measure %in% c("rsq", "r", "rmse"))) stop("unknown type.measure: ", type.measure)
  
  rmse <- function(x1, x2) sqrt(mean((x1-x2)^2))
  
  #cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel)
  cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel, ...)
  
  rs     <- sapply(1:length(cvfit$lambda), function(i) cor(cvfit$fit.preval[,i], y))
  rsqs   <- rs^2
  rmses  <- sapply(1:length(cvfit$lambda), function(i) rmse(cvfit$fit.preval[,i], y))
  cvfit$measures <- list(
    r = rs, 
    rsq = rsqs, 
    rmse = rmses
  )
  
  cvfit$bestfit <- get_bestfit(cvfit, type.measure)
  
  return(cvfit)
}

# TODO
run_cvglmnet.binomial <- function(X, y, keep=T, parallel=T, type.measure="accuracy", ...) 
{
  if (!(type.measure %in% c("accuracy"))) stop("unknown type.measure: ", type.measure)
  
  cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel, ...)
  
  accuracies <- sapply(1:length(cvfit$lambda), function(i) cor(cvfit$fit.preval[,i], y))
  
  rs     <- sapply(1:length(cvfit$lambda), function(i) cor(cvfit$fit.preval[,i], y))
  rsqs   <- rs^2
  rmses  <- sapply(1:length(cvfit$lambda), function(i) rmse(cvfit$fit.preval[,i], y))
  cvfit$measures <- list(
    r = rs, 
    rsq = rsqs, 
    rmse = rmses
  )
  
  cvfit$bestfit <- get_bestfit(cvfit, type.measure)
  
  return(cvfit)
}

run_caret_glmnet <- function(X, y, nfolds=10, nrepeats=2, nlambda=25, alpha=1, 
                             metric="Rsquared", family="gaussian")
{
  fitControl <- trainControl(
    method = "repeatedcv",
    number = nfolds, 
    repeats = nrepeats, 
    allowParallel = TRUE, 
    savePredictions = "final"
  )
  
  # Tuning Grids
  gfit    <- glmnet(X, y, family=family, nlambda=nlambda, alpha=alpha)
  lambdas <- gfit$lambda
  grids   <- data.frame(alpha = alpha, lambda=lambdas)
  
  # Fit
  if (is.null(colnames(X))) colnames(X) <- 1:ncol(X)
  fit <- train(X, y, 
               method = "glmnet",
               trControl = fitControl, 
               tuneGrid = grids, 
               #tuneLength = 20, 
               preProcess = c("center", "scale"), 
               family=family, 
               metric=metric)
  
  return(fit)
}


# Load Shape/Texture Features ---------------------------------------------

basedir <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"

sym.shapes      <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames    <- rownames(sym.shapes)
sym.shapes.eigs <- read.table(file.path(basedir, "shape_sym_pca_eigs.txt"))
sym.textures    <- read.csv(file.path(basedir, "texture_sym_pca_scores.csv"), row.names = 1)

#nosym.shapes <- read.csv(file.path(basedir, "shape_pca_scores.csv"), row.names = 1)
#nosym.shapes.eigs <- read.table(file.path(basedir, "shape_pca_eigs.txt"))
#nosym.textures  <- read.csv(file.path(basedir, "texture_pca_scores.csv"), row.names = 1)


# Load Predictors ---------------------------------------------------------

# Load demographics information
base <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
demo.vnames <- demos$video # typo
demo.vnames <- sub("_fr[0-9]{3}", "", demo.vnames)
# Reorder rows to match features
oinds       <- sapply(shape.vnames, function(x) which(demo.vnames == x))
all.equal(demo.vnames[oinds], shape.vnames)
df.demos    <- demos[oinds,-c(1:2)]

# Load trait information
base         <- "/data1/famface01/analysis/encoding/12_Features"
traits       <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
trait.vnames <- as.character(traits[,1])
traits       <- traits[,-1]
# Reorder
oinds       <- sapply(shape.vnames, function(x) which(trait.vnames == x))
all.equal(trait.vnames[oinds], shape.vnames)
df.traits   <- traits[oinds,]

levels(df.demos$facial_hair)


# Distracted ---------------------------------------------------------

# Trying to make the df.demos into a thing to be used for regression

formula_to_mat <- function(formula, data) {
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  return(X)
}
df.demos <- within(df.demos, facial_hair <- relevel(facial_hair, ref = "None"))
tmp <- formula_to_mat(~facial_hair, data=df.demos)
head(tmp)


# Factor Analysis ---------------------------------------------------------

trait.fa <- factanal(df.traits, factors = 6, rotation = "varimax", 
                     na.action = na.omit, scores="regression")

# Reorder the rows of the rotation matrix for nicer viewing
hc1     <- hclust(dist(trait.fa$loadings))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- trait.fa$loadings[ord.hc1,]
round(loading, 3)
library(corrplot)
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)




# CCA ---------------------------------------------------------------------

library(PMA)

#' ### Centered Trait Scores
#' 
#+ cancor-analysis-traits

# Permute to get penalty
X1 <- as.matrix(sym.textures)[,1:200] # the first 200 comps is enough i think
X2 <- as.matrix(sym.shapes) # might not really need all the shapes
colnames(X2) <- sprintf("shape_%s", colnames(X2))
colnames(X1) <- sprintf("texture_%s", colnames(X1))
sym.combined  <- cbind(X2, X1)
perm.out <- CCA.permute(df.traits, sym.combined, typex="standard",typez="standard", nperms=25)
print(perm.out)
plot(perm.out)

# Get CCA, 9 is best we can do
out <- CCA(df.traits, sym.combined, typex="standard",  typez="standard", K=10,
           penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz, 
           v=perm.out$v.init)
print(out)

rownames(out$u) <- colnames(df.traits)
round(out$u, 2)
round(out$v, 2)

hc1     <- hclust(dist(out$u))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- out$u[ord.hc1,]
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

xscores <- scale(df.traits, scale=F) %*% out$u
yscores <- scale(sym.combined, scale=F) %*% out$v

tmp <- out$v %*% MASS::ginv(t(out$v))
head(tmp[,1:5])
head(scale(sym.combined, scale=F)[,1:5])

# Classify Shapes/Texture ---------------------------------------------------

# We use both to allow for rendering the full face

# Setup
fits <- list()
X1 <- as.matrix(sym.textures)[,1:200] # the first 200 comps is enough i think
X2 <- as.matrix(sym.shapes) # might not really need all the shapes
colnames(X2) <- sprintf("shape_%s", colnames(X2))
colnames(X1) <- sprintf("texture_%s", colnames(X1))
X  <- cbind(X2, X1)
## this is for later when we extract the relavant information
inds.shape   <- 1:ncol(X2)
inds.texture <- (1:ncol(X))[-inds.shape]

# Gender
y <- factor(df.demos$gender)
fit <- run_caret_glmnet(X, y, nfolds=10, nrepeats=5, nlambda=100, alpha=1, 
                        metric="Accuracy", family="binomial")
## get accuracy
confmat <- confusionMatrix(fit$pred$pred, fit$pred$obs)
confmat$overall[1]
## get coefficients
tmp <- coef(fit$finalModel, s=fit$bestTune$lambda)
tmp <- tmp[-1] # this will be used to weight the face in
tmp2<- varImp(fit)
## save
fits$gender <- fit

# Age
y <- as.numeric(df.demos$age) # cuz this is group average
fit <- run_caret_glmnet(X, y, nfolds=10, nrepeats=5, nlambda=100, alpha=1, 
                        metric="Rsquared", family="gaussian")
fit$results[rownames(fit$bestTune),]
#tmp <- fit$pred[grep("Rep1", fit$pred$Resample),]
#plot.ts(tmp[,c("pred","obs")], plot.type = "single", col=2:3)
## save
fits$age <- fit

# Attractiveness
y <- as.numeric(df.traits$attractive)
fit <- run_caret_glmnet(X, y, nfolds=10, nrepeats=5, nlambda=100, alpha=1, 
                        metric="Rsquared", family="gaussian")
fit$results[rownames(fit$bestTune),]
#tmp <- fit$pred[grep("Rep1", fit$pred$Resample),]
#plot.ts(tmp[,c("pred","obs")], plot.type = "single", col=2:3)
## save
fits$attractive <- fit

# Factor Scores (6 of them)
fitsFA <- llply(1:ncol(trait.fa$scores), function(i) {
  y <- as.numeric(trait.fa$scores[,i])
  fit <- run_caret_glmnet(X, y, nfolds=10, nrepeats=5, nlambda=100, alpha=1, 
                          metric="Rsquared", family="gaussian")
  fit
}, .progress="text")
ret <- ldply(1:length(fitsFA), function(i) {
  fit <- fits[[i]]
  data.frame(factor=sprintf("factor %i", i), fit$results[rownames(fit$bestTune),])
})
ret # not all traits have a great fit

fits$trait_factors <- fitsFA

# Save
save(fits, file="/data1/famface01/command/misc/face_representations/120_features/tmp_combined_fits.rda")


# Extract Coefs ------------------------------------------------------------

# the -1 removes the intercept (all 1s)
cfs.gender <- coef(fits$gender$finalModel, s=fits$gender$bestTune$lambda)[-1]
cfs.age    <- coef(fits$age$finalModel, s=fits$age$bestTune$lambda)[-1]
cfs.attractive <- coef(fits$attractive$finalModel, s=fits$attractive$bestTune$lambda)[-1]

# combine
cfs <- cbind(cfs.gender, cfs.age, cfs.attractive)
cfs.shape   <- cfs[inds.shape,]
cfs.texture <- cfs[inds.texture,]
head(cfs.texture)

# save
basedir <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
write.table(cfs.shape, file=file.path(basedir, "coefs_age+gender+attractive_shape.txt"))
write.table(cfs.texture, file=file.path(basedir, "coefs_age+gender+attractive_texture.txt"))




# Inverse Model ------------------------------------------------------------

library(investr)

# So I want to use the coefficients to constrain the forward model

# Select the coefs
X_age <- X[,cfs.age!=0]

# Do the inverse regression
y_age <- as.numeric(df.demos$age) # cuz this is group average
inv_fit_age <- lm(X ~ y_age)
predict(inv_fit_age, data.frame(y_age=0))

# Simply save this regression model to be used in python
cfs_age <- cfs.age
inds_shape <- inds.shape
inds_texture <- inds.texture
save(X, X_age, y_age, inv_fit_age, cfs_age, inds_shape, inds_texture, 
     file=file.path(basedir, "age+gender+attractive_inverse_model.rda"))









# Redo the regression based on those
y.age <- as.numeric(df.demos$age) # cuz this is group average
fit.age <- lm(y.age ~ X.age)
invest(fit.age, 0.5, interval="inversion", x0.name="shape_Comp002", newdata=as.data.frame(df.demos$age))

#b <- solve(t(X.age) %*% X.age) %*% t(X.age) %*% y.age
b <- fit.age$coefficients
invert.b <- t(b) %*% MASS::ginv(b %*% t(b))
invert.X <- y.age %*% invert.b # is that ok?
## fit
cor(invert.X[,-1], X.age)
mean(sapply(1:ncol(X.age), function(i) cor(X.age[,i], invert.X[,i+1])))


xx <- X.age[,1:2]
fit <- lm(y.age ~ xx)
fit$coefficients
y = 2.4*x0 + 0.04*x1 + -0.04*x2

rmse <- function(x,y) sqrt(sum((x-y)^2))

# Let's see how to get the best estimate
## calculate
xx <- X.age[,1:2]
fit <- lm(y.age ~ xx)
b <- fit$coefficients
## invert
invert.b <- t(b) %*% MASS::ginv(b %*% t(b))
invert.X <- y.age %*% invert.b # is that ok?
## test
cor(invert.X[,2],xx[,1])
cor(invert.X[,3],xx[,2])
rmse(invert.X[,2],xx[,1])
rmse(invert.X[,3],xx[,2])
rmse(invert.X[,2]+invert.X[,1]/2,xx[,1])
rmse(invert.X[,3]+invert.X[,1]/2,xx[,2])


# Finally let's try doing a simple multiple regression in reverse
xx <- X.age[,1:2]
inv.fit <- lm(xx ~ y.age)
invert.X <- predict(inv.fit, data.frame(y.age=y.age))
cor(invert.X[,1],xx[,1])
cor(invert.X[,2],xx[,2])
rmse(invert.X[,1],xx[,1])
rmse(invert.X[,2],xx[,2])


# Try on everything
fit <- lm(y.age ~ X.age)
b   <- fit$coefficients
invert.b <- t(b) %*% MASS::ginv(b %*% t(b))
invert.X1 <- y.age %*% invert.b # is that ok?
mean(sapply(1:ncol(X.age), function(i) rmse(invert.X1[,i+1], X.age[,i])))
mean(sapply(1:ncol(X.age), function(i) rmse(invert.X1[,i+1]+invert.X1[,1], X.age[,i])))

inv.fit <- lm(X.age ~ y.age)
invert.X2 <- predict(inv.fit, data.frame(y.age=y.age))
mean(sapply(1:ncol(X.age), function(i) rmse(invert.X2[,i], X.age[,i])))

plot(sapply(1:ncol(X.age), function(i) rmse(invert.X1[,i+1], X.age[,i])), 
     sapply(1:ncol(X.age), function(i) rmse(invert.X2[,i], X.age[,i])))
