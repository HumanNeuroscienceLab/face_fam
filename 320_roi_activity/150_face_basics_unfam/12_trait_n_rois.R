#' ---
#' title: "R Notebook Stub"
#' output:
#'  html_document:
#'    fig_width: 10
#'    fig_height: 5
#' ---
#'
#' We will be showing some information on the trait structure.
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

suppressMessages(library(gplots))

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


#' # ROI Regression Analysis
#' 
#' I will now do ROI analysis using the centered trait scores, PCA trait scores, 
#' and Factor Analytic trait scores.
#' 
#' Note that I tried the beta-series and voxelwise data but both didn't perform
#' that well. Instead, here I use the average timeseries (voxelwise).
#' 
#' ## Convolutions
#' 
#' We prepare the convolved time-series first.
#' 
#+ roi-convolve
suppressMessages(library(bigmemory))
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


#' ## Centered Trait Measures
#'
#' 
#+ roi-trait
# Regress
ret <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traits[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
tnames <- colnames(dat.vols$sub01$features$traits)
dimnames(ret) <- list(subj=subjects, trait=c("face", tnames), roi=dimnames(ret)[[3]])

# Average
grp.ret <- apply(ret, 2:3, mean)
round(grp.ret, 3)

#' ### Trait Plot
#' 
#+ roi-trait-plot

# Quick bar plot
barplot(grp.ret, beside=T, col=brewer.pal(11,"Spectral"), 
        legend.text=T, xlim=c(0,length(rnames)*12+20), args.legend = list(
          x = length(rnames)*12 + 33, 
          y = max(grp.ret), 
          bty = 'n'
        ), ylab="T-Statistic", main="Centered Trait Scores")
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)

# Only right hemi
barplot(grp.ret[,1:5], beside=T, col=brewer.pal(11,"Spectral"), 
        legend.text=T, xlim=c(0,length(rnames)*7+10), args.legend = list(
          x = length(rnames)*7 + 17, 
          y = max(grp.ret[1:5,]), 
          bty = 'n'
        ), ylab="T-Statistic", main="Centered Trait Scores (RH)")
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


#' ## PCA
#'
#' 
#+ roi-trait-pca

ret <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitspca[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
tnames <- colnames(dat.vols$sub01$features$traits)
dimnames(ret) <- list(subj=subjects, trait=c("face", tnames), roi=dimnames(ret)[[3]])

# Average
grp.ret <- apply(ret, 2:3, mean)
round(grp.ret, 3)

#' ### Trait Factor Plot
#' 
#+ roi-trait-pca-plot

# Quick bar plot
barplot(grp.ret, beside=T, col=brewer.pal(11,"Spectral"), 
        legend.text=T, xlim=c(0,length(rnames)*12+20), args.legend = list(
          x = length(rnames)*12 + 34, 
          y = max(grp.ret), 
          bty = 'n'
        ), ylab="T-Statistic", main="Trait PCA Scores")
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)

# Only right hemi
barplot(grp.ret[,1:5], beside=T, col=brewer.pal(11,"Spectral"), 
        legend.text=T, xlim=c(0,length(rnames)*7+10), args.legend = list(
          x = length(rnames)*7 + 17, 
          y = max(grp.ret[1:5,]), 
          bty = 'n'
        ), ylab="T-Statistic", main="Trait PCA Scores (RH)")
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)


#' ## Factor Analysis
#'
#' 
#+ roi-trait-fa

ret <- laply(subjects, function(subj) {
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsfa[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  
  sapply(dat$fmri$dat, function(rdat) {
    fit <- simple_lm(rowMeans(rdat), convs)
    fit$tvals
  })
})
tnames <- sprintf("factor %i", 1:6)
dimnames(ret) <- list(subj=subjects, trait=c("face", tnames), roi=dimnames(ret)[[3]])

# Average
grp.ret <- apply(ret, 2:3, mean)
round(grp.ret, 3)

#' ### Trait FA Plot
#' 
#+ roi-trait-fa-plot

# Quick bar plot
barplot(grp.ret, beside=T, col=brewer.pal(7,"Spectral"), 
        legend.text=c("Face", sprintf("Factor %i", 1:6)), 
        xlim=c(0,length(rnames)*7+19), args.legend = list(
          x = length(rnames)*7 + 27, 
          y = max(grp.ret), 
          bty = 'n'
        ), ylab="T-Statistic", main="Trait Factor Scores")
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)

# Only right hemi
barplot(grp.ret[,1:5], beside=T, col=brewer.pal(7,"Spectral"), 
        legend.text=c("Face", sprintf("Factor %i", 1:6)), 
        xlim=c(0,length(rnames)*3+19), args.legend = list(
          x = length(rnames)*3 + 23, 
          y = max(grp.ret), 
          bty = 'n'
        ), ylab="T-Statistic", main="Trait Factor Scores (RH)")
abline(h=1.96, lty=3)
abline(h=-1.96, lty=3)





# THIS IS COOL STUFF
#ret <- laply(subjects, function(subj) {
  subj <- "sub02"
  dat   <- dat.vols[[subj]]
  convs <- conv.traitsfa[[subj]]
  convs <- cbind(face=dat$convolve$face, convs)
  colnames(convs) <- c("face", sprintf("factor_%02i", 1:6))
  rdats <- sapply(dat$fmri$dat, rowMeans)
  head(rdats)
#})


library(mvtboost)
out2 <- mvtb(Y=rdats[,1:5], X=convs, 
             n.trees=250, 
             shrinkage=.01,
             interaction.depth=1:3,
             
             bag.fraction=.5,      # fit each tree to a sub sample of this fraction
             train.fraction=.5,    # only fit the model to this fraction of the data set
             cv.folds=10,           # number of cross-validation folds
             mc.cores=20,           # run the cross-validation in parallel
             seednum=103)
summary(out2)
nonlin.out <- mvtb.nonlin(out2,X=convs,Y=rdats[,1:5])
nonlin.out$r.EBA$rank.list
covex      <- mvtb.covex(out2, X=convs, Y=rdats[,1:5])
round(covex,2)
cc <- mvtb.cluster(covex, clust.method = "ward.D", dist.method = "manhattan")
round(cc,2)
mvtb.heat(covex, cexRow=0.6)

round(mvtb.ri(out2, relative = "tot"), 2) 
numformat <- function(val){sub("^(-?)0.", "\\1.", sprintf("%.1f", val))}

par(mar=c(8,10,1,1))
#col <- colorRampPaletteAlpha(RColorBrewer::brewer.pal(9,"Greys"),100)
col <- brewer.pal(9, "YlOrRd")[2:9]
ri.out2 <- t(mvtb.ri(out2))
mvtb.heat(ri.out2, cexRow=1, cexCol=1, numformat=numformat, col=col)
par()
