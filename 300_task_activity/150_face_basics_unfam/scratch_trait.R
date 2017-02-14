#' This script will load the beta data and regress the trait PCAs onto it

# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(doMC))
registerDoMC(24)
library(RColorBrewer)
library(dynamicTreeCut)

library(ggplot2)
source("/data1/famface01/command/encoding/FirstPass/plot_functions.R")
source("/data1/famface01/command/encoding/SubRois_Unfam/lm_functions.R")

library(caret)
library(glmnet)
library(gplots)

subjects <- sprintf("sub%02i", 1:6) # again skipping sub01
indir    <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
rnames   <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA",
              "l.vATL", "l.FFA", "l.OFA", "l.EBA")

source("/data1/famface01/command/encoding/ShapesAnalysis/000_funs_load_convolve.R")


# Functions ---------------------------------------------------------------

load.beta.data <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/betas_roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}

# Not needed but here as ref
load.traits <- function(beta.vids) {
  base      <- "/data1/famface01/analysis/encoding/12_Features"
  
  # Read in
  features  <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
  feat.vids <- as.character(features[,1])
  features  <- features[,-1]
  
  # match beta vid order with feat vid order
  inds <- sapply(beta.vids, function(vid) which(feat.vids == vid))
  if (!all.equal(feat.vids[inds], beta.vids)) stop("rearrange failure")
  
  features <- features[inds,]
  rownames(features) <- feat.vids[inds]
  
  return(features)
}

# loads the voxelwise data with relevant information like the raw features
load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}
# tmp <- load.cur("sub01")
# table(tmp$basics$timing$question)

load.trait.pca <- function(dat) {
  traits.df <- dat$features$traits
  trait.pca <- prcomp(traits.df, scale=F, retx=T)
  return(as.matrix(trait.pca$x))
}

getncomps.iter <- function(x, k, nsim=100, scale=FALSE, parallel=TRUE) {
  library(svd)
  library(MASS)
  r     <- nrow(x)
  c     <- ncol(x)
  if (scale) x <- scale(x)
  true.evals <- trlan.svd(x, k)$d
  perm.evals <- laply(1:nsim, function(i) {
    y <- mvrnorm(n=r, mu=rep(0,c), Sigma=diag(1,c), empirical=F)
    evals <- trlan.svd(y, k)$d
    evals
  }, .parallel=parallel)
  pvals <- sapply(1:length(true.evals), function(ci) {
    vec <- c(true.evals[ci], perm.evals[,ci]) # combine true with permuted
    sum(vec[1]<=vec)/length(vec)
  })
  ncomps  <- sum(pvals<0.05)
  return(ncomps)
}



# Load Data ---------------------------------------------------------------

dats <- lapply(subjects, load.beta.data)
names(dats) <- subjects


# Load trait data ---------------------------------------------------------

# The dataframe should be the same across participants, just re-arranged
traits.df <- dats$sub01$features$traits

# We can visualize it
library(corrplot)
col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                          "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                          "#4393C3", "#2166AC", "#053061"))(200)
corrplot(cor(traits.df), order = "hclust", tl.col='black', tl.cex=.75, diag=F, 
         col=rev(col)) 

# Perform PCA
trait.pca <- prcomp(traits.df, scale=F, retx=T)
dim(trait.pca$x)
## visualize the components
corrplot(trait.pca$rotation, tl.col='black', tl.cex=.75, diag=T, col=rev(col))


# Run for one -------------------------------------------------------------

ret1 <- sapply(subjects, function(subj) {
  rdat <- as.matrix(dats[[subj]]$betas.glms$r.aFFA)
  tdat <- load.trait.pca(dats[[subj]])
  fit <- simple_lm(rdat, tdat)
  rowMeans(fit$tvals)
})

ret2 <- sapply(subjects, function(subj) {
  rdat <- as.matrix(dats[[subj]]$betas.glms$r.pFFA)
  tdat <- load.trait.pca(dats[[subj]])
  fit <- simple_lm(rdat, tdat)
  rowMeans(fit$tvals)
})

ret3 <- sapply(subjects, function(subj) {
  rdat <- as.matrix(rowMeans(dats[[subj]]$betas.glms$r.pFFA))
  tdat <- load.trait.pca(dats[[subj]])
  fit <- simple_lm(rdat, tdat)
  fit$tvals
})
ret3


library(bigmemory)
dat <- load.cur("sub01")
convs <- convolve.features.byvid(dat$basics, dat$features$traits)
all.equal(convs, dat$convolve$traits)
dim(convs)

fit1 <- simple_lm(dat$fmri$dat$r.pFFA, convs)
fit2 <- simple_lm(dat$fmri$dat$r.pFFA, dat$convolve$traits)

rowMeans(fit1$tvals)
rowMeans(fit2$tvals)

# This is much better! So averaging the voxels first and then doing the stats is the way to go.
fit1 <- simple_lm(rowMeans(dat$fmri$dat$r.pFFA), convs)
fit2 <- simple_lm(rowMeans(dat$fmri$dat$r.pFFA), dat$convolve$traits)
fit1$tvals
fit2$tvals

# ok yes this looks like the results although seems I get other comps as well
trait.comps <- load.trait.pca(dat)
convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
fit3 <- simple_lm(rowMeans(dat$fmri$dat$r.pFFA), convs.comps)
fit3$tvals


# Now let's get this for all subjs
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects

## Look at the original trait measures
ret1 <- sapply(dat.vols, function(dat) {
  convs <- convolve.features.byvid(dat$basics, dat$features$traits)
  fit <- simple_lm(rowMeans(dat$fmri$dat$r.pFFA), convs)
  fit$tvals
})
rownames(ret1) <- colnames(dat.vols$sub01$features$traits)
colnames(ret1) <- subjects


## now look at the PCA
ret2 <- sapply(dat.vols, function(dat) {
  trait.comps <- load.trait.pca(dat)
  convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
  fit <- simple_lm(rowMeans(dat$fmri$dat$r.pFFA), convs.comps)
  fit$tvals
})
rownames(ret2) <- sprintf("Comp%02i", 1:10)
colnames(ret2) <- subjects

ret1
ret2

round(rowMeans(ret1), 4)
round(rowMeans(ret2), 4)



###
# CONVOLUTIONS
###


# Now let's get this for all subjs
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects

# Get the convolved features
conv.traits <- lapply(dat.vols, function(dat) {
  convs <- convolve.features.byvid(dat$basics, dat$features$traits)
  convs
})
conv.traits2 <- lapply(dat.vols, function(dat) {
  convs <- convolve.features.byvid(dat$basics, scale(dat$features$traits, scale=F))
  convs
})
conv.traits3 <- lapply(dat.vols, function(dat) {
  convs <- convolve.features.byvid(dat$basics, scale(dat$features$traits))
  convs
})

conv.traitspca <- lapply(dat.vols, function(dat) {
  trait.comps <- load.trait.pca(dat)
  convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
  convs.comps
})
conv.traitspca2 <- lapply(dat.vols, function(dat) {
  trait.comps <- load.trait.pca(dat)
  traits.df   <- dat$features$traits
  trait.pca   <- prcomp(traits.df, scale=T, retx=T)
  trait.comps <- trait.pca$x
  convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
  convs.comps
})

conv.traitsfa <- lapply(dat.vols, function(dat) {
  traits.df <- dat$features$traits
  res1 <- factanal(traits.df, factors = 6, rotation = "varimax", na.action = na.omit, 
                   scores="regression")
  convs.comps <- convolve.features.byvid(dat$basics, res1$scores)
  convs.comps
})
## sparse PCA
library(PMA)
traits.df <- dats$sub01$features$traits
cv.out <- SPC.cv(as.matrix(traits.df), sumabsvs=seq(1, sqrt(ncol(traits.df)), len=20), niter=10, orth = T)
print(cv.out)
conv.traitspma <- lapply(dat.vols, function(dat) {
  traits.df <- as.matrix(dat$features$traits)
  out <- SPC(traits.df, sumabsv=cv.out$bestsumabsv, K=10, orth = T)
  rownames(out$v) <- colnames(traits.df)
  convs.comps <- convolve.features.byvid(dat$basics, out$u)
  convs.comps
})

## sparse FA - NoCor
library(fanc)
fancs1 <- llply(dat.vols, function(dat) {
  traits.df <- as.matrix(dat$features$traits)
  ret <- fanc(scale(as.matrix(traits.df), scale=F), 5, cor.factor=F)
  ret2 <- select(ret, criterion="BIC", gamma=Inf, scores=T)
  ret2$scores
}, .parallel=T)
conv.traitsfanc1 <- lapply(1:length(subjects), function(i) {
  convs.comps <- convolve.features.byvid(dat.vols[[i]]$basics, fancs1[[i]])
  convs.comps
})
names(conv.traitsfanc1) <- subjects
## sparse FA - Cor
fancs2 <- llply(dat.vols, function(dat) {
  traits.df <- as.matrix(dat$features$traits)
  ret <- fanc(scale(as.matrix(traits.df), scale=F), 5, cor.factor=T)
  ret2 <- select(ret, criterion="BIC", gamma=Inf, scores=T)
  ret2$scores
}, .parallel=T)
conv.traitsfanc2 <- lapply(1:length(subjects), function(i) {
  convs.comps <- convolve.features.byvid(dat.vols[[i]]$basics, fancs2[[i]])
  convs.comps
})
names(conv.traitsfanc2) <- subjects

## Non-Matrix Factorization - 5
library(NMF)
conv.traitsnmf5 <- lapply(dat.vols, function(dat) {
  traits.df   <- dat$features$traits
  nres        <- nmf(as.matrix(traits.df), 5)
  convs.comps <- convolve.features.byvid(dat$basics, nres@fit@W)
  convs.comps
})
conv.traitsnmf7 <- lapply(dat.vols, function(dat) {
  traits.df   <- dat$features$traits
  nres        <- nmf(as.matrix(traits.df), 7)
  convs.comps <- convolve.features.byvid(dat$basics, nres@fit@W)
  convs.comps
})
conv.traitsnmf10 <- lapply(dat.vols, function(dat) {
  traits.df   <- dat$features$traits
  nres        <- nmf(as.matrix(traits.df), 10)
  convs.comps <- convolve.features.byvid(dat$basics, nres@fit@W)
  convs.comps
})
# ICA
library(fastICA)
conv.traitsica3 <- lapply(dat.vols, function(dat) {
  traits.df   <- scale(as.matrix(dat$features$traits), center=T, scale=F)
  trait.comps <- fastICA(traits.df,  3, method="C")$S
  convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
  convs.comps
})
conv.traitsica5 <- lapply(dat.vols, function(dat) {
  traits.df   <- scale(as.matrix(dat$features$traits), center=T, scale=F)
  trait.comps <- fastICA(traits.df,  5, method="C")$S
  convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
  convs.comps
})
conv.traitsica7 <- lapply(dat.vols, function(dat) {
  traits.df   <- scale(as.matrix(dat$features$traits), center=T, scale=F)
  trait.comps <- fastICA(traits.df,  7, method="C")$S
  convs.comps <- convolve.features.byvid(dat$basics, trait.comps)
  convs.comps
})


## to get loadings
#ret <- fanc(as.matrix(traits.df), 5, cor.factor=T)
#ret <- fanc(as.matrix(traits.df), 5, cor.factor=F) 
#select a model via model selection criterion
#ret2 <- select(ret, criterion="BIC", gamma=Inf, scores=T)
#ret2$loadings
## to get loadings
#out <- SPC(as.matrix(traits.df), sumabsv=cv.out$bestsumabsv, K=10, orth = T)
#rownames(out$v) <- colnames(traits.df)
#round(out$v, 4)
### to get the loadings
#traits.df <- dats$sub01$features$traits
#res1 <- factanal(traits.df, factors = 6, rotation = "varimax", na.action = na.omit, 
#                 scores="regression")  
#res1$loadings



###
# REGRESSIONS
###

# Run the traits
ret1 <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traits[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
tnames <- colnames(dat.vols$sub01$features$traits)
dimnames(ret1) <- list(subj=subjects, trait=c("face", tnames), roi=dimnames(ret1)[[3]])

# Quick bar plot
barplot(apply(ret1, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the traits with centering
ret1b <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traits2[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
tnames <- colnames(dat.vols$sub01$features$traits)
dimnames(ret1b) <- list(subj=subjects, trait=c("face", tnames), roi=dimnames(ret1b)[[3]])

# Quick bar plot
barplot(apply(ret1b, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the traits with centering and scaling
ret1c <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traits3[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
tnames <- colnames(dat.vols$sub01$features$traits)
dimnames(ret1c) <- list(subj=subjects, trait=c("face", tnames), roi=dimnames(ret1c)[[3]])

# Quick bar plot
barplot(apply(ret1c, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the trait PCA
ret2 <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitspca[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret2) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:10)), roi=dimnames(ret2)[[3]])

# Quick bar plot
barplot(apply(ret2, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)

# Run the trait PCA with scaling
ret2b <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitspca2[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret2b) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:10)), roi=dimnames(ret2b)[[3]])

# Quick bar plot
barplot(apply(ret2b, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)



# Run the trait FA
ret3 <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsfa[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret3) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:6)), roi=dimnames(ret3)[[3]])

# Quick bar plot
barplot(apply(ret3, 2:3, mean), beside=T, col=brewer.pal(7,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)

# Run the sparse PCA
ret4 <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitspma[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret4) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:10)), roi=dimnames(ret4)[[3]])

# Quick bar plot
barplot(apply(ret4, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the sparse FA-1
ret5a <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsfanc1[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret5a) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:5)), roi=dimnames(ret5a)[[3]])

# Quick bar plot
barplot(apply(ret5a, 2:3, mean), beside=T, col=brewer.pal(6,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the sparse FA-2
ret5b <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsfanc2[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret5b) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:5)), roi=dimnames(ret5b)[[3]])

# Quick bar plot
barplot(apply(ret5b, 2:3, mean), beside=T, col=brewer.pal(6,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the NMF 5
ret6a <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsnmf5[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret6a) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:5)), roi=dimnames(ret6a)[[3]])

# Quick bar plot
barplot(apply(ret6a, 2:3, mean), beside=T, col=brewer.pal(6,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the NMF 7
ret6b <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsnmf7[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret6b) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:7)), roi=dimnames(ret6b)[[3]])

# Quick bar plot
barplot(apply(ret6b, 2:3, mean), beside=T, col=brewer.pal(8,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# Run the NMF 10
ret6c <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsnmf10[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret6c) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:10)), roi=dimnames(ret6c)[[3]])

# Quick bar plot
barplot(apply(ret6c, 2:3, mean), beside=T, col=brewer.pal(11,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# ICA - 3 comps
ret7a <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsica3[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret7a) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:3)), roi=dimnames(ret7a)[[3]])

# Quick bar plot
barplot(apply(ret7a, 2:3, mean), beside=T, col=brewer.pal(4,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# ICA - 5 comps
ret7b <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsica5[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret7b) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:5)), roi=dimnames(ret7b)[[3]])

# Quick bar plot
barplot(apply(ret7b, 2:3, mean), beside=T, col=brewer.pal(6,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


# ICA - 7 comps
ret7c <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsica7[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
dimnames(ret7c) <- list(subj=subjects, trait=c("face",sprintf("Comp-%02i", 1:7)), roi=dimnames(ret7b)[[3]])

# Quick bar plot
barplot(apply(ret7c, 2:3, mean), beside=T, col=brewer.pal(8,"Spectral"))
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)



###
# Canonical Correlation
###

dat <- dat.vols$sub01

# Get the data together
lst.rdats <- lapply(dat.vols, function(dat) {
  rdats <- sapply(dat$fmri$dat, rowMeans)
  scale(rdats, scale=F)
})
rdats <- do.call(rbind, lst.rdats)

# Get the trait data (regular and factor analysis)
lst.convs <- lapply(conv.traits, scale, scale=F)
lst.convs <- lapply(1:length(lst.convs), function(i) cbind(face=dat.vols[[i]]$convolve$face, conv.traits[[i]]))
convs1 <- do.call(rbind, lst.convs)
colnames(convs1) <- colnames(dat.vols$sub01$features$traits)
lst.convs <- lapply(conv.traitsfa, scale, scale=F)
convs2 <- do.call(rbind, lst.convs)

cc.res <- cancor(rdats, convs1)
dim(cc.res$xcoef)
dim(cc.res$ycoef)
dim(convs1)
barplot(scale(cc.res$xcoef), beside = T, col=brewer.pal(9, "Spectral"), legend=T)
barplot(scale(cc.res$ycoef), beside = T, col=brewer.pal(11, "Spectral"), legend=T)


perm.out <- CCA.permute(rdats, convs1, typex="standard",typez="standard", nperms=25)
print(perm.out)
plot(perm.out)
out <- CCA(rdats,convs1,typex="standard",typez="standard",K=9,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)
print(out)

rownames(out$u) <- rnames
rownames(out$v) <- colnames(dat.vols$sub01$features$traits)
round(out$u, 4)
round(out$v, 4)
barplot(out$u, beside = T, col=brewer.pal(9, "Spectral"), legend=F)
barplot(out$v, beside = T, col=brewer.pal(9, "Spectral"), legend=F)
