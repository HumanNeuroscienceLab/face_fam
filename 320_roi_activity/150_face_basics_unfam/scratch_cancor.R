
#' We will be doing the canonical correlation analysis.
#' 
#' # Setup
#' 
#+ setup
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


#' ## Functions
#' 
#+ setup-functions

# NOTE: betas not working well
load.beta.data <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/betas_roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}

# Not needed but here as ref
load.traits <- function() {
  base      <- "/data1/famface01/analysis/encoding/12_Features"
  
  features  <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
  feat.vids <- as.character(features[,1])
  features  <- features[,-1]
  rownames(features) <- feat.vids
  
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

#' ## Load Data
#' 
#+ setup-load-data

dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects


#' # Convolutions
#' 
#' We prepare the convolved time-series first.
#' 
#+ roi-convolve
library(bigmemory)
base <- "/data1/famface01/command/misc/face_representations/300_task_activity/150_face_basics_unfam/"
rdafile <- file.path(base, "42_trait_n_rois_convolve.rda")
if (file.exists(rdafile)) {
  load(rdafile)
} else {
  conv.traits <- lapply(dat.vols, function(dat) {
    convs <- convolve.features.byvid(dat$basics, scale(dat$features$traits, scale=F))
    convs
  })
  conv.traitspca <- lapply(dat.vols, function(dat) {
    trait.comps <- load.trait.pca(dat)
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
  save(conv.traits, conv.traitspca, conv.traitsfa, file=rdafile)
}


#' # ROI Cananical Correlation Analysis
#' 
#' We now run a sparse canonical correlation analysis to see the correspondence
#' between the trait scores and the ROI activity.
#' 
#' ## Gather data
#' 
#' We average the time-series in each ROI and concatenate across subjects.
#' 
#+ cancor-data

# Timeseries Data
lst.rdats <- lapply(dat.vols, function(dat) {
  rdats <- sapply(dat$fmri$dat, rowMeans)
  scale(rdats, scale=F)
})
rdats <- do.call(rbind, lst.rdats)

# Gather centered and factor data convolution
lst.convs <- lapply(conv.traits, scale, scale=F)
lst.convs <- lapply(1:length(lst.convs), function(i) cbind(face=dat.vols[[i]]$convolve$face, lst.convs[[i]]))
convs1 <- do.call(rbind, lst.convs)
colnames(convs1) <- c("face", colnames(dat.vols$sub01$features$traits))

lst.convs <- lapply(conv.traitsfa, scale, scale=F)
lst.convs <- lapply(1:length(lst.convs), function(i) cbind(face=dat.vols[[i]]$convolve$face, lst.convs[[i]]))
convs2 <- do.call(rbind, lst.convs)
colnames(convs2) <- c("face", sprintf("factor %i", 1:6))


#' ## Canonical Correlation
#' 
#+ cancor-analysis
library(PMA)

#' ### Centered Trait Scores
#' 
#+ cancor-analysis-traits

# Permute to get penalty
perm.out <- CCA.permute(rdats, convs1, typex="standard",typez="standard", nperms=25)
print(perm.out)
plot(perm.out)

# Get CCA, 9 is best we can do
out <- CCA(rdats, convs1, typex="standard",  typez="standard", K=9,
           penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz, 
           v=perm.out$v.init)
print(out)

par(mfrow=c(1,1))
rownames(out$u) <- rnames
rownames(out$v) <- c("face", colnames(dat.vols$sub01$features$traits))
round(out$u, 2)
round(out$v, 2)

#barplot(out$u, beside = T, col=brewer.pal(9, "Spectral"), legend=F)
#barplot(out$v, beside = T, col=brewer.pal(9, "Spectral"), legend=F)


#' ### Factor Analytic Scores
#' 
#+ cancor-analysis-traits-fa

# Permute to get penalty
perm.out <- CCA.permute(rdats, convs2, typex="standard",typez="standard", nperms=25)
print(perm.out)
plot(perm.out)

# Get CCA, 9 is best we can do
out <- CCA(rdats, convs2, typex="standard",  typez="standard", K=7,
           penaltyx=perm.out$bestpenaltyx, penaltyz=perm.out$bestpenaltyz, 
           v=perm.out$v.init)
print(out)

par(mfrow=c(1,1))
rownames(out$u) <- rnames
rownames(out$v) <- c("face", sprintf("factor %i", 1:6))
round(out$u, 2)
round(out$v, 2)

#' Show the factor analysis loadings
#' 
#+ more
trait.df <- load.traits()
trait.fa <- factanal(trait.df, factors = 6, rotation = "varimax", 
                     na.action = na.omit, scores="regression")
round(unclass(trait.fa$loadings), 3)


# y <- rowMeans(rdat)
# x <- conv.traitsfa$sub02
# fit=hierNet.path(x,y)
# fitcv=hierNet.cv(fit,x,y)
# print(fitcv)
# plot(fitcv)
# ret <- hierNet(x, y, lam=fitcv$lamlist[which.min(fitcv$cv.err)])
# ret
