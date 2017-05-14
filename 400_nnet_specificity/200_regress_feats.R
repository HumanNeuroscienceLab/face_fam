
# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(R.matlab)
library(doMC)
registerDoMC(30)
library(plyr)
library(RColorBrewer)
library(bigmemory)


# Load info and regressors ------------------------------------------------

# For each frame in the video, we extracted the features from each layer of the neural network
indir <- "/data1/famface01/analysis/misc/openface/layer_features"
infiles <- list.files(indir)
vidname <- sub("_fr[0-9]{3}.mat", "", infiles)
frame   <- as.integer(sub(".*fr", "", sub(".mat", "", infiles)))
vdf     <- data.frame(ind=1:length(infiles), vid=vidname, fr=frame, fn=infiles)
## select only the unfamiliar faces
base <- "/data1/famface01/analysis/encoding/12_Features"
tmp     <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
vdf     <- vdf[vdf$vid %in% as.character(tmp$X),]
vdf$ind <- 1:nrow(vdf)

# Load demographics information
base <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
demo.vnames <- sub("_fr[0-9]{3}", "", demos$video)
# Reorder rows to match features
oinds       <- sapply(vdf$vid, function(x) which(demo.vnames == x))
all.equal(demo.vnames[oinds], as.character(vdf$vid))
df.demos    <- demos[oinds,-c(1:2)]

# Load trait information
base         <- "/data1/famface01/analysis/encoding/12_Features"
traits       <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
trait.vnames <- as.character(traits[,1])
traits       <- traits[,-1]
# Reorder
oinds       <- sapply(vdf$vid, function(x) which(trait.vnames == x))
all.equal(trait.vnames[oinds], as.character(vdf$vid))
df.traits   <- traits[oinds,]

# Get pose information
pbase <- "/data1/famface01/data/stimuli/vids/eight_frames/three_dee"
df.pose <- ldply(unique(vdf$vid), function(vname) {
  fpath <- sprintf("%s/%s/%s_pose.txt", pbase, vname, vname)
  tab <- read.table(fpath)
  colnames(tab) <- c("pitch", "yaw", "roll")
  tab$fd <- Rfast::dista(c(0,0,0), t(tab))
  tab$vid <- vname
  tab
}, .progress="text")
all.equal(as.character(df.pose$vid), as.character(vdf$vid))
df.pose <- df.pose[,colnames(df.pose) != "vid"]

# Get luminance information
base <- "/data1/famface01/analysis/misc/400_nnet_specificity"
df.lum <- read.csv(file.path(base, "feats_ave_luminance.csv"))
lum.vnames <- as.character(df.lum$vid)
df.lum <- df.lum[,-1]
# Reorder
oinds       <- sapply(1:nrow(vdf), function(i) {
  which((lum.vnames == vdf$vid[i]) & (df.lum$frame == vdf$fr[i]))
})
all.equal(lum.vnames[oinds], as.character(vdf$vid))
df.lum   <- df.lum[oinds,"luminance",drop=F]

# Combine
df.info <- cbind(vdf, df.pose, df.lum, df.demos, df.traits)

### Load Raw Data (Images)
#library(png)
#imgdir <- "/data1/famface01/analysis/misc/openface/aligned_frames"
#imgfiles <- sub(".mat", ".png", df.info$fn)
#imgs <- laply(imgfiles, function(imgfile) {
#  img <- readPNG(file.path(imgdir, imgfile))
#  as.vector(0.2126*img[,,1] + 0.7152*img[,,2] + 0.0722*img[,,3])
#}, .progress="text")
## plot(1:2, type='n')
## rasterImage(as.raster(img), 1, 1, 2, 2, interpolate=T)

## old way
reps <- read.csv("/data1/famface01/analysis/misc/openface/reps.csv", header=F)
labs <- read.csv("/data1/famface01/analysis/misc/openface/labels.csv", header=F)
labs$fn <- sub("png", "mat", basename(as.character(labs$V2)))



# Load data ---------------------------------------------------------------

nlayers  <- 26
bdir <- "/data1/famface01/analysis/misc/openface/concatenate_layers"

# Subset if you want
vdf2     <- df.info
vdf2$vid <- factor(vdf2$vid)

# Load it all
if (length(list.files(bdir, pattern="bin$")) != nlayers) {
  feats <- llply(vdf2$fn, function(infile) {
    featlist <- readMat(file.path(indir, infile))
    if (length(featlist) != nlayers) stop("featlist len")
    featlist
  }, .progress="text")
} else {
  feats <- NULL
}

# Concatenate each layer across subjects and save
concat_layer <- function(num_layer, flatten=T) {
  lname  <- sprintf("layer%02i", num_layer)
  clayer <- laply(feats, function(featlist) featlist[[lname]])
  
  if (flatten) dim(clayer) <- c(dim(clayer)[1], prod(dim(clayer)[-1]))
  
  return(clayer)
}
bigconcat_layer <- function(num_layer, verbose="none") {
  lname  <- sprintf("layer%02i", num_layer)
  bfile  <- sprintf("c%s.bin", lname)
  dfile  <- sprintf("c%s.desc", lname)
  
  # get the end dimensions
  nr <- length(feats)
  indiv_dim <- dim(feats[[1]][[lname]])
  nc <- prod(indiv_dim)
  
  # create the file-backed big matrix
  bm <- big.matrix(nr, nc, backingpath=bdir, backingfile=bfile, 
                   descriptorfile=dfile)
  
  # loop through and add each row
  l_ply(1:nr, function(ri) {
    bm[ri,] <- as.vector(feats[[ri]][[lname]])
  }, .progress=verbose)
  
  # make sure all is written
  flush(bm)
  
  return(bm)
}

clayers <- llply(1:nlayers, function(li) {
  cat("layer", li, "\n")
  dfile <- sprintf("clayer%02i.desc", li)
  
  if (!file.exists(file.path(bdir, dfile))) {
    bm   <- bigconcat_layer(li, "text")
    flush(bm); rm(bm)
  }
  
  bm     <- attach.big.matrix(file.path(bdir, dfile))
  
  bm
})

rm(feats); gc(reset=T) # if feats is in the system, remove it to free memory


# Cluster Metrics for Each Layer ------------------------------------------

#registerDoMC(9)
accs <- laply(1:nlayers, function(li) {
  cat(li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  
  dmat   <- Rfast::Dist(X)
  nearest.vid <- sapply(1:nrow(dmat), function(i) vdf2$vid[-i][which.min(dmat[i,-i])])
  
  # percent accuracy
  id.acc <- sum(diag(table(nearest.vid, vdf2$vid)))/nrow(dmat) # accuracy
  # wb ratio
  within.dist <- daply(vdf2, .(vid), function(x) {
    d0 <- dmat[x$ind,x$ind]
    mean(d0[lower.tri(d0)])
  })
  between.dist <- daply(vdf2, .(vid), function(x) {
    d0 <- dmat[x$ind,-x$ind]
    mean(d0)
  }) 
  wb.ratio <- mean(within.dist)/mean(between.dist)
  # silhoutte width
  sii <- cluster::silhouette(as.numeric(vdf2$vid), dmatrix=dmat)
  sc <- summary(sii)
  avg.silwidth <- sc$avg.width
  # combine
  cstats <- c(accuracy=id.acc, win=within.dist, btw=between.dist, wb.ratio=wb.ratio, avg.silwidth=avg.silwidth)
  
  cstats
}, .parallel=T)


### Just for the raw input data
#dmat   <- Rfast::Dist(imgs)
#nearest.vid <- sapply(1:nrow(dmat), function(i) vdf2$vid[-i][which.min(dmat[i,-i])])
#
## percent accuracy
#id.acc <- sum(diag(table(nearest.vid, vdf2$vid)))/nrow(dmat) # accuracy
## wb ratio
#within.dist <- daply(vdf2, .(vid), function(x) {
#  d0 <- dmat[x$ind,x$ind]
#  mean(d0[lower.tri(d0)])
#})
#between.dist <- daply(vdf2, .(vid), function(x) {
#  d0 <- dmat[x$ind,-x$ind]
#  mean(d0)
#}) 
#wb.ratio <- mean(within.dist)/mean(between.dist)
## silhoutte width
#sii <- cluster::silhouette(as.numeric(vdf2$vid), dmatrix=dmat)
#sc <- summary(sii)
#avg.silwidth <- sc$avg.width
## combine
#cstats <- c(accuracy=id.acc, wb.ratio=wb.ratio, avg.silwidth=avg.silwidth)


### PLAY

fac.res <- psych::fa(vdf2[,18:ncol(vdf2)], nfactors=5, residuals=T, rotate='varimax', fm='minres')
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

knn.dmat <- sapply(1:nrow(dmat), function(i) mean(sort(dmat[i,vdf2$vid[i]!=vdf2$vid])[1:70]))

round(cor(cbind(avg.dmat, avg.dmat2, knn.dmat, fac.scores)), 3)
tmp <- cbind(avg.dmat, avg.dmat2, knn.dmat, fac.scores)
tmp <- apply(tmp, 2, function(xx) tapply(xx, vdf2$vid, mean))
round(cor(tmp), 3)

cor(cbind(avg.dmat, knn.dmat2, knn.dmat, knn.dmat3, knn.dmat4), fac.scores)

library(apcluster)
fac.red <- fac.scores[vdf2$fr==3,]
xred <- apply(X, 2, function(xx) tapply(xx, vdf2$vid, mean))
xred.dmat <- Rfast::Dist(xred)
clust <- apcluster(negDistMat(r=2), xred)
mds <- MASS::isoMDS(as.dist(xred.dmat))
plot(clust, mds$points)
plot(mds$points, col=as.numeric(vdf2$gender[vdf2$fr==3])) # it's all gender


# Regress with Pose/Demos/Trait -------------------------------------------

library(caret)
library(glmnet)
registerDoMC(20)

run_repeated_cvglmnet <- NULL
source("/data1/famface01/command/misc/face_representations/misc/cv_glmnet.R")

outdir <- "/data1/famface01/analysis/misc/400_nnet_specificity"

select.layers <- c(3,4,5,8,11,12,13,14:19,21:22,25,26) # we don't run everything
select.layers2 <- c(4, 16, 26)

gen.folds.vid <- function(y, k=8, nreps=5) {
  yuniq <- y[vdf2$fr==3] # get only one value in the frame
  
  # get folds for each video
  folds <- caret::createMultiFolds(yuniq, k=k, times=nreps)
  cfolds <- sapply(1:nreps, function(i) {
    names <- sprintf("Fold%i.Rep%i", 1:k, i)
    vec <- vector("numeric", length(yuniq))
    for (fi in 1:length(names)) {
      vec[-folds[[names[fi]]]] <- fi
    }
    vec
  })
  
  # expand to include each frame
  cfolds <- cfolds[as.numeric(vdf2$vid),]
  
  return(cfolds)
}

k     <- 8; nreps <- 5


###
# Pos: Gender
###

y     <- vdf2$gender
cfolds<- gen.folds.vid(y, k, nreps)

# run it all
fits.gender <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$gender
  cvrep.gender <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                        family="binomial", type.measure="Accuracy", 
                                        cfolds=cfolds, k=k, nreps=nreps)
  cvrep.gender
})
# load(file.path(outdir, "classify_layers_gender.rda"))
gender.accs <- sapply(fits.gender, function(x) max(x$mean.max.res))

mat <- cbind(accs[select.layers,-2], gender.accs, sapply(clayers[select.layers], dim)[2,])
colnames(mat) <- c("Clust Acc", "Sil Width", "Gender Acc", "# Features")
plot.ts(mat[,-4], xlab="Layer", main=NA)
plot.ts(scale(mat), xlab="Layer", main=NA, plot.type = "single", col=1:4)
select.layers[13]


###
# Pos: Age
###

y     <- vdf2$age
cfolds<- gen.folds.vid(y, k, nreps)

# run it all
fits.age <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$age
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})


###
# Pos: Typicality
###

y     <- vdf2$typical
cfolds<- gen.folds.vid(y, k, nreps)

# run it all
fits.typical <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$typical
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})
# load(file.path(outdir, "classify_layers_typical.rda"))
typical.accs <- sapply(fits.typical, function(x) max(x$mean.max.res))
plot.ts(typical.accs) # not very good


###
# Pos: Attractive
###

y     <- vdf2$attractive
cfolds<- gen.folds.vid(y, k, nreps)

# run it all
fits.attractive <- llply(select.layers2, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$attractive
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})


###
# Neg: Pose
###

y     <- vdf2$fd
cfolds<- gen.folds.vid(y, k, nreps)

fits.pose <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$fd
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})


###
# Neg: Luminance
###

y     <- vdf2$age
cfolds<- gen.folds.vid(y, k, nreps)

# run it all
fits.age <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$age
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})

y     <- vdf2$luminance
cfolds<- gen.folds.vid(y, k, nreps)

fits.luminance <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- vdf2$luminance
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})


###
# Neg: Makeup (independent of gender)
###

y     <- lm(makeup ~ gender, data=vdf2)$residuals
cfolds<- gen.folds.vid(y, k, nreps)

# run it all
fits.makeup <- llply(select.layers, function(li) {
  cat("layer", li, "\n")
  
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  
  cvrep <- run_repeated_cvglmnet(X, y, alpha=0, parallel=T, 
                                 family="gaussian", type.measure="rsq", 
                                 cfolds=cfolds, k=k, nreps=nreps)
  cvrep
})



# Extra Stuff -------------------------------------------------------------

# Don't use this anymore

pca_layer <- function(clayer) {
  clayer <- Rfast::standardise(clayer, center=T, scale=F)
  sk <- 2; ret <- list(d=NULL, u=NULL)
  for (k in seq(50,600,by=50)) {
    cat(k, "..", sep="")
    ret <- trlan.svd(clayer, k, lambda=ret$d, U=ret$u)
    ret$d2 <- ret$d/sqrt(max(1, nrow(clayer2) - 1))
    nfacs <- nFactors::nScree(ret$d2, cor=F)
    sk    <- max(nfacs$Components)
    if (sk < (k-50)) break
  }
  cat("\n")
  clayer <- ret$u[,1:sk,drop=F]
  clayer
}

run_caret_glmnet <- function(X, y, nfolds=10, nrepeats=10, nlambda=100, alpha=1, 
                             metric="Rsquared", family="gaussian", ...)
{
  fitControl <- trainControl(
    method = "repeatedcv",
    number = nfolds, 
    repeats = nrepeats, 
    allowParallel = TRUE, 
    savePredictions = "final", 
    classProbs = TRUE
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

fits.gender <- llply(1:nlayers, function(li) {
  X     <- as.matrix(clayers[[li]])
  Xvars <- Rfast::colVars(X)
  X     <- X[,Xvars!=0]
  y     <- df.info$gender
  fit   <- run_caret_glmnet(X, y, alpha=0.5, metric="Accuracy", family="binomial", 
                            nrepeats=5, nlambda=25)
  fit
})
