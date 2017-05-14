
# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(doMC)
registerDoMC(30)
library(plyr)
library(bigmemory)

library(caret)
library(nlme)

library(ggplot2)
library(ggthemr)
library(RColorBrewer)

ggthemr('pale', type='outer', text_size=14, layout='plain')

subjects <- sprintf("sub%02i", 1:6)



# Functions ---------------------------------------------------------------

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

extract_probs <- function(fit, y) {
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], y)
  mean.preds <- ddply(preds, .(rowIndex), colwise(mean, levels(y)))
  return(mean.preds)
}

extract_preds <- function(fit, y) {
  preds <- fit$pred
  preds$Rep <- sub("Fold[0-9]+[.]", "", preds$Resample)
  preds <- preds[order(preds$Rep, preds$rowIndex),]
  all.equal(preds$obs[preds$Rep=="Rep01"], y)
  mode.preds <- ddply(preds, .(rowIndex), function(x) {
    c(pred=names(which.max(table(x$pred))))
  })
  mode.preds
}



# Load nnet + behav -------------------------------------------------------

indir <- "/data1/famface01/command/misc/face_representations/322_asap_rois/210_fampics_setup"

## nnets
face.norm <- faces.avg <- faces.exp <- faces.fullset <- faces.exp.targ <- faces.exp2 <- faces.avg.targ <- NULL
load(file.path(indir, "z_face_nnet_feats.rda"), verbose=T)

## behav
behav.df <- read.csv(file.path(indir, "z_likeness+rt+gender+set.csv"))[,-1]



# Prototype ---------------------------------------------------------------

###
# Collect Distances
###

# Get distances between everything
grps <- rep(1:4, c(nrow(face.norm), nrow(faces.avg.targ$feats), 
                   nrow(faces.exp2$feats), nrow(faces.fullset$feats)))
mat  <- rbind(face.norm, faces.avg.targ$feats, faces.exp2$feats, faces.fullset$feats)
dmat <- Rfast::Dist(mat)

# Save dists to norm
norm.dist <- dmat[grps==1,grps==3]

# Save dists to each ave (only keep the targets)
norm.avgs <- dmat[grps==2,grps==3]
rownames(norm.avgs) <- as.character(faces.avg.targ$df$person)
## only keep within identity average
norm.avgs2 <- norm.avgs
for (i in 1:nrow(norm.avgs)) {
  pname <- rownames(norm.avgs)[i]
  pinds <- faces.exp2$df$person == pname
  norm.avgs2[i,!pinds] <- 0
}
## remove the identity average
norm.avgs3 <- norm.avgs2
for (i in 1:nrow(norm.avgs)) {
  pname <- rownames(norm.avgs)[i]
  pinds <- faces.exp2$df$person == pname
  norm.avgs3[i,pinds] <- scale(norm.avgs2[i,pinds], center=T, scale=F)
}
## transpose
norm.avgs <- t(norm.avgs)
#norm.avgsX <- lm(norm.avgs ~ person, data=behav.df)$residuals
norm.avgs2 <- t(norm.avgs2)
norm.avgs3 <- t(norm.avgs3)
norm.avgs4 <- rowSums(norm.avgs3)


###
# Regress likeness
###

# First, let's see if the distances to the norm or to each average helps
# !!! Distance to norm is better than dists to each average, strange
fit <- lm(likeness ~ person + norm.dist, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + norm.avgs4, data=behav.df)
summary(aov(fit))

# Second, does it help to calculate the distances from an average to different identities (first one)?
# !!! No, second case is better (better looking only within the identity)
fit <- lm(likeness ~ person + norm.avgs, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + norm.avgs4, data=behav.df)
summary(aov(fit))

# Third, what if we combine the two (dist to norm and dist to averages). 
# Any added value?
## !!! Yes, improvement of fit with both.
fit1 <- lm(likeness ~ person + norm.dist, data=behav.df)
summary(aov(fit1))
fit2 <- lm(likeness ~ person + norm.dist + norm.avgs4, data=behav.df)
summary(aov(fit2))
anova(fit1, fit2)


###
# Regress RT
###

# First, let's see if the distances to the norm or to each average helps
# !!! Distance to norm is better than dists to each average, strange
fit <- lm(rt ~ person + norm.dist, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + norm.avgs4, data=behav.df)
summary(aov(fit))

# Second, does it help to calculate the distances from an average to different identities (first one)?
# !!! It doesn't seem to matter
fit <- lm(rt ~ person + norm.avgs, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + norm.avgs4, data=behav.df)
summary(aov(fit))

# Third, what if we combine the two (dist to norm and dist to averages). 
# Any added value?
## !!! No, adding the distances to each identity, doesn't help
fit1 <- lm(rt ~ person + norm.dist, data=behav.df)
summary(aov(fit1))
fit2 <- lm(rt ~ person + norm.dist + norm.avgs4, data=behav.df)
summary(aov(fit2))
anova(fit1, fit2)

# TODO: maybe having the distance to the second closest might help predict the RTs?


###
# Viz
###

# TODO: heat map of the distances btw averages (or maybe plotting them)
# Accuracy in the classification accuracy of the identities
# The model fits
# In some other section want to see if the RTs or likeness differ on average 
# betweenn the identities (seems yes, which means a familiarity effect)

# only 1 error
table(factor(faces.exp2$df$person), apply(norm.avgs, 1, which.min))



# kNN ---------------------------------------------------------------

###
# Classify
###

library(caret)
registerDoMC(30)

# restrict to eight IDs
targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
             "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
faces.full.targ <- list()
faces.full.targ$df <- faces.fullset$df[faces.fullset$df$person %in% targets,]
faces.full.targ$df$person <- factor(faces.full.targ$df$person)
faces.full.targ$feats <- faces.fullset$feats[faces.fullset$df$person %in% targets,]

# classify
X <- faces.full.targ$feats; 
y <- faces.full.targ$df$person
knnFit <- run_caret(X, y, "knn", tlen=20)

# preds/probs (not used just to see)
knn.probs <- extract_probs(knnFit, y)
knn.preds <- extract_preds(knnFit, y)
table(knn.preds$pred, y)

# test preds/probs (will save to run regression)
knn.probs <- predict(knnFit$finalModel, faces.exp2$feats, "prob")
knn.preds <- predict(knnFit$finalModel, faces.exp2$feats, "class")

faces.exp2$df$person <- factor(faces.exp2$df$person)
table(knn.preds, faces.exp2$df$person)

# simply knn measure
k <- knnFit$bestTune[[1]]
Xtest <- Rfast::Dist(faces.exp2$feats)
kdists <- sapply(1:nrow(Xtest), function(i) {
  mean(sort(as.numeric(Xtest[i,-i]))[1:k])
})


###
# Regress for likeness
###

# Combine probs together to be only that identity
knn.probs2 <- sapply(1:ncol(knn.probs), function(i) {
  name <- colnames(knn.probs)[i]
  inds <- faces.exp2$df$person == name
  x <- knn.probs[,i] * (inds*1)
  x[inds] <- scale(x[inds], center=T, scale=F)
  x
})
colnames(knn.probs2) <- colnames(knn.probs)
knn.probs3 <- rowSums(knn.probs2)

# First, let's see if modeling the probs together is helpful
# !!! Yes, #3 is combined and only within identity
fit <- lm(likeness ~ person + knn.probs, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + knn.probs2, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + knn.probs3, data=behav.df)
summary(aov(fit))

# Second, let's see if we looks at the knn overall dists
# This gives us how close our neighbors are or not
## !!! knn dist measure does very well, no need for class probs
fit1 <- lm(likeness ~ person + kdists, data=behav.df)
summary(aov(fit1))
fit2 <- lm(likeness ~ person + kdists + knn.probs3, data=behav.df)
summary(aov(fit2))
anova(fit1, fit2)


###
# Regress for RT
###

# First, let's see if modeling the probs together is helpful
# !!! None of them are good!
fit <- lm(rt ~ person + knn.probs, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + knn.probs2, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + knn.probs3, data=behav.df)
summary(aov(fit))

# Second, let's see if we looks at the knn overall dists
# This gives us how close our neighbors are or not
## !!! None of them are good!
fit1 <- lm(rt ~ person + kdists, data=behav.df)
summary(aov(fit1))
fit2 <- lm(rt ~ person + kdists + knn.probs3, data=behav.df)
summary(aov(fit2))
anova(fit1, fit2)




# SVM ---------------------------------------------------------------

###
# Classify
###

library(caret)
library(e1071)
registerDoMC(30)

# classify
X <- faces.full.targ$feats; 
y <- faces.full.targ$df$person
svmFit <- run_caret(X, y, "svmLinear2", tlen=20)

# preds/probs (not used just to see)
svm.probs <- extract_probs(svmFit, y)
svm.preds <- extract_preds(svmFit, y)
table(svm.preds$pred, y)

# test preds/probs (will save to run regression)
pred <- predict(svmFit$finalModel, faces.exp2$feats, 
                decision.values=T, probability=T)
head(attr(pred, "decision.values"))
head(attr(pred, "probabilities"))
svm.dvals <- attr(pred, "decision.values")
svm.probs <- attr(pred, "probabilities")

# combine dvals
# there are 7 dvals per identity
unames <- levels(faces.exp2$df$person)
svm.dvals.ave <- sapply(1:length(unames), function(i) {
  tmp <- svm.dvals[,grep(unames[i], colnames(svm.dvals))]
  rowMeans(tmp)
})
colnames(svm.dvals.ave) <- unames
svm.dvals.min <- sapply(1:length(unames), function(i) {
  tmp <- svm.dvals[,grep(unames[i], colnames(svm.dvals))]
  apply(tmp, 1, function(x) x[which.min(x)])
})
colnames(svm.dvals.min) <- unames
svm.dvals.in <- t(sapply(1:nrow(svm.dvals), function(i) {
  cname <- unames[faces.exp2$df$person[i]==unames]
  vec <- vector("numeric", 8)
  
  x <- svm.dvals[i,grep(cname, colnames(svm.dvals))]
  
  cur.names <- sub("/", "", sub(cname, "", names(x)))
  inds <- sapply(cur.names, function(cn) which(unames==cn))
  vec[inds] <- x
  names(vec) <- unames
  
  vec
}))
colnames(svm.dvals.in) <- unames

# collapse probs
svm.probs2 <- sapply(colnames(svm.probs), function(name) {
  x       <- svm.probs[,colnames(svm.probs)==name]
  inds    <- (faces.exp2$df$person==name)
  x       <- x * inds
  x[inds] <- scale(x[inds], center=T, scale=T)
  x
})
svm.probs2 <- rowSums(svm.probs2)

# collapse dvals
svm.dvals.ave2 <- sapply(colnames(svm.dvals.ave), function(name) {
  x       <- svm.dvals.ave[,colnames(svm.dvals.ave)==name]
  inds    <- (faces.exp2$df$person==name)
  x       <- x * inds
  x[inds] <- scale(x[inds], center=T, scale=T)
  x
})
svm.dvals.ave2 <- rowSums(svm.dvals.ave2)
svm.dvals.min2 <- sapply(colnames(svm.dvals.min), function(name) {
  x       <- svm.dvals.ave[,colnames(svm.dvals.min)==name]
  inds    <- (faces.exp2$df$person==name)
  x       <- x * inds
  x[inds] <- scale(x[inds], center=T, scale=T)
  x
})
svm.dvals.min2 <- rowSums(svm.dvals.min2)


###
# Regress for likeness
###

# !!! Not Significant
fit <- lm(likeness ~ person + svm.probs2, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + svm.dvals.ave2, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + svm.dvals.min2, data=behav.df)
summary(aov(fit))

# Want to have more than one regressor
# !!! The Dvals ave is significant
fit <- lm(likeness ~ person + svm.probs, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + svm.dvals.ave, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + svm.dvals.min, data=behav.df)
summary(aov(fit))
fit <- lm(likeness ~ person + svm.dvals.in, data=behav.df)
summary(aov(fit))

# !!! Dvals ave is good enough
fit1 <- lm(likeness ~ person + svm.dvals.ave, data=behav.df)
summary(aov(fit1))
fit2 <- lm(likeness ~ person + svm.probs, data=behav.df)
summary(aov(fit2))
fit3 <- lm(likeness ~ person + svm.probs + svm.dvals.ave, data=behav.df)
summary(aov(fit3))
anova(fit1,fit3) # not beneficial to add probs to dvals
anova(fit2,fit3) # good to add dvals to probs


###
# Regress for RT
###

fit <- lm(rt ~ person + svm.probs2, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + svm.dvals.ave2, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + svm.dvals.min2, data=behav.df)
summary(aov(fit))

fit <- lm(rt ~ person + svm.probs, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + svm.dvals.ave, data=behav.df)
summary(aov(fit))
fit <- lm(rt ~ person + svm.dvals.min, data=behav.df)
summary(aov(fit))



# Compare -----------------------------------------------------------------

###
# Likeness
###

# Let's now compare the different effects in a combined model
## !!! kdists is dropping out here
fit <- lm(likeness ~ person + norm.dist + norm.avgs4 + kdists + svm.dvals.ave, data=behav.df)
summary(aov(fit))

# Do it one by one
fit1 <- lm(likeness ~ person + norm.dist, data=behav.df)
summary(aov(fit1))
fit2 <- lm(likeness ~ person + norm.dist + norm.avgs4, data=behav.df)
summary(aov(fit2))
anova(fit1,fit2)
## not significat to add the knn
fit3 <- lm(likeness ~ person + norm.dist + norm.avgs4 + kdists, data=behav.df)
summary(aov(fit3))
anova(fit2,fit3)
## yes adds something with svm distances
fit4 <- lm(likeness ~ person + norm.dist + norm.avgs4 + svm.dvals.ave, data=behav.df)
summary(aov(fit4))
anova(fit2,fit4)


###
# RT
###

# Both dist to norm and dist to boundaries important for RT
fit1 <- lm(rt ~ person + norm.dist, data=behav.df)
summary(aov(fit1))
fit2 <- lm(rt ~ person + norm.dist + svm.dvals.ave, data=behav.df)
summary(aov(fit2))
anova(fit1,fit2)


ref <- summary(lm(rt ~ person, data=behav.df))$adj.r.squared
summary(lm(rt ~ person + norm.dist, data=behav.df))$adj.r.squared - ref
summary(lm(rt ~ person + svm.dvals.ave, data=behav.df))$adj.r.squared - ref
summary(lm(rt ~ person + norm.dist + svm.dvals.ave, data=behav.df))$adj.r.squared - 
# can we do some quick classification?
