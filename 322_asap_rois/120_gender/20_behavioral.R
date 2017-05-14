
# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(doMC)
registerDoMC(30)
library(plyr)
library(bigmemory)

subjects <- sprintf("sub%02i", 1:6)



# Load --------------------------------------------------------------------

# regressors etc
fdf <- read.csv("120_gender/z_df_info.csv")
fdf <- fdf[,-1]
head(fdf)

# face features
bdir   <- "/data1/famface01/analysis/misc/openface/concatenate_layers"
dfile  <- sprintf("clayer%02i.desc", 26)
bm     <- attach.big.matrix(file.path(bdir, dfile))
feats  <- as.matrix(bm)

# reduce to only vids
vdf    <- fdf[fdf$fr==3,]
vfeats <- apply(feats, 2, function(x) tapply(x, fdf$vid, mean))

# read in average
base <- "/data1/famface01/analysis/misc/openface"
labels <- read.csv(sprintf('%s/masked_labels_redo.csv', base), header=F)
features <- read.csv(sprintf('%s/masked_reps_redo.csv', base), header=F)
## find the average face
aind <- grep("average_face", labels$V2)
avg.all <- as.numeric(features[aind,])

# read in gender average
base <- "/data1/famface01/analysis/misc/openface/prototypes_gender"
labels <- read.csv(sprintf('%s/labels.csv', base), header=F)
features <- read.csv(sprintf('%s/reps.csv', base), header=F)
avg.male <- as.numeric(features[1,])
avg.female <- as.numeric(features[2,])

# cor the prototypes
round(cor(cbind(avg.all, avg.male, avg.female)), 3)



# Prototype 1 -------------------------------------------------------------

# This uses the above avgs as the prototypes

# calculate distance to averages
dist.to.avg    <- Rfast::Dist(rbind(avg.all, vfeats))[-1,1]
dist.to.male   <- Rfast::Dist(rbind(avg.male, vfeats))[-1,1]
dist.to.female <- Rfast::Dist(rbind(avg.female, vfeats))[-1,1]

# run a regression to see how well these things fit the data
fit <- glm(vdf2$gender ~ dist.to.avg + dist.to.male + dist.to.female, family=binomial(link='logit'))
summary(fit)




# Prototype 2 -------------------------------------------------------------

# This will use k-means clustering to get the centroids

library(MASS)

# Run a k-means on the data
cl <- kmeans(vfeats, 2, iter.max=200, nstart=20)

# calculate distance to averages
kdist.to.male   <- Rfast::Dist(rbind(cl$centers[2,], vfeats))[-1,1]
kdist.to.female <- Rfast::Dist(rbind(cl$centers[1,], vfeats))[-1,1]

# Calculate some distances?
# so how well do the distances to the center

# Visualize
dat  <- rbind(vfeats,cl$centers,avg.all,avg.male,avg.female)
#dat[2561:2568,] <- rnorm(ncol(dat)) # need to replace since identical to: 5497:5504
#dat[5473:5480,] <- rnorm(ncol(dat)) # 5473:5480 and 5505:5512
dmat <- Rfast::Dist(dat)
d    <- as.dist(dmat)
mds  <- isoMDS(d)
plot(mds$points[1:nrow(vfeats),])
n <- nrow(vfeats) + 1
points(mds$points[c(n,n+1),], col=2:3, cex=3)
points(mds$points[c(n+2,n+3,n+4),], col=1:3, cex=3, pch=6)




# Exemplar ----------------------------------------------------------------

nrepeats <- 5
nfolds <- 10
fitControl <- trainControl(
  method = "repeatedcv", 
  number = nfolds, 
  repeats = nrepeats, 
  returnResamp = "final", 
  savePredictions = "final", 
  classProbs = T, 
  allowParallel = T
)

y <- vdf2$gender
X <- vfeats
colnames(X) <- sprintf("feat%02i", 1:ncol(X))
#knnFit <- train(X, y, method = "glmnet", trControl = fitControl, preProcess = c("center","scale"), tuneLength = 10)
knnFit <- train(X, y, method = "knn", trControl = fitControl, preProcess = c("center","scale"), tuneLength = 20)
ri <- as.numeric(rownames(knnFit$bestTune))
k <- knnFit$results$k[ri]
varImp(knnFit)

dmat <- Rfast::Dist(vfeats)
ret <- t(sapply(1:nrow(dmat), function(i) {
  o <- order(dmat[i,-i])[1:k]
  tapply(vfeats[o,],  )
  prop.table(table(vdf2$gender[o]))
}))

# I'm trying a classifier to determine the most important features

library(caret)
library(glmnet)
registerDoMC(24)

# subs.folds
nrepeats <- 5
nfolds <- 10

fitControl <- trainControl(
  method = "repeatedcv", 
  number = nfolds, 
  repeats = nrepeats, 
  returnResamp = "final", 
  savePredictions = "final", 
  classProbs = T, 
  allowParallel = T
)

X <- cbind(avg=dist.to.avg, male=dist.to.male, female=dist.to.female, 
           kmale=kdist.to.male, kfemale=kdist.to.female)
alphas <- c(0,0.5,1)
tuneGrid <- ldply(alphas, function(alpha) {
  tmp <- glmnet(X, vdf2$gender, family="binomial", nlambda=50, alpha=alpha)
  data.frame(alpha=alpha, lambda=tmp$lambda)
}, .parallel=T)

fit <- train(X, vdf2$gender, 
             method = "glmnet", 
             trControl = fitControl, 
             preProcess = c("center", "scale"), 
             tuneGrid = tuneGrid)
ri <- as.numeric(rownames(fit$bestTune))
fit$results[ri,]

varImp(fit)
