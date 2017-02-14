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


#' # Trait Data
#' 
#' Here we visualize the trait data and its correlations.
#' 
#' ## Load
#' 
#+ trait-load

trait.df <- load.traits()
head(round(trait.df, 2))

#' ## Heatmap
#' 
#' This will re-order both rows and columns
#' 
#+ trait-heatmap

smat <- scale(trait.df)
zlims <- max(abs(smat)) * c(-1,1)
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                            "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                            "#4393C3", "#2166AC", "#053061"))(200))
heatmap(smat, col=col, scale="none", zlim=zlims, labRow = NA, margins=c(7,0))

#col2 <- rev(brewer.pal(11, "RdBu"))
#hc.cols <- hclust(dist(t(smat)), method="ward.D2")
#hc.rows <- hclust(dist(smat), method="ward.D2")
#heatmap(smat, Rowv=as.dendrogram(hc.rows), Colv=as.dendrogram(hc.cols), 
#        col=col2, scale="none", zlim=zlims, labRow = NA)

#' ## Correlations
#' 
#' We can see how correlated the different measures are.
#' 
#+ trait-correlation
library(corrplot)
corrplot(cor(trait.df), order = "hclust", tl.col='black', tl.cex=.75, diag=F, 
         col=col) 
corrplot(abs(cor(trait.df)), order = "hclust", tl.col='black', tl.cex=.75, diag=F, 
         col=col) 
round(cor(trait.df), 3)


#' # Data Reduction
#' 
#' ## PCA
#' 
#' Note that the first principal component is the one that explains everything
#' in the fMRI analysis.
#' 
#+ dr-pca

# Do PCA
trait.pca <- prcomp(trait.df, scale=F, retx=T)

#' ### Percent Variance Explained
#' 
#+ dr-pca-var-explained

prop.cvar <- cumsum((trait.pca$sdev)^2) / sum(trait.pca$sdev^2)
plot(x=1:10, y=prop.cvar*100, type='l', xlab="Components", 
     ylab="Percent Variance Explained", ylim=c(0,100), col='blue')
points(x=1:10, y=prop.cvar*100)

#' ### Rotation Matrix
#' 
#+ dr-pca-rotmat

# Reorder the rows of the rotation matrix for nicer viewing
hc1 <- hclust(dist(t(trait.df)))
hc1 <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
rotmat <- trait.pca$rotation[ord.hc1,]
round(rotmat, 3)

# Plot weights
corrplot(rotmat, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=F) 


#' ## Factor Analysis
#' 
#' Not totally sure what I was doing. Did the standard options. I couldn't go 
#' higher than 6 components, otherwise the program would crash. 
#' 
#+ dr-fa

# Do Factor Analysis
trait.fa <- factanal(trait.df, factors = 6, rotation = "varimax", 
                     na.action = na.omit, scores="regression")

#' ### Percent Variance Explained
#' 
#' We now plot the percent of variance explained
#' 
#+ dr-fa-var-explained 
prop.cvar <- cumsum(colSums(unclass(trait.fa$loadings)^2)/ncol(trait.df))
plot(x=1:6, y=prop.cvar*100, type='l', xlab="Factors", col="blue", 
     ylab="Percent Variance Explained", ylim=c(0,100))
points(x=1:6, y=prop.cvar*100)

#' ### Loadings
#' 
#+ dr-fa-loadings

# Reorder the rows of the rotation matrix for nicer viewing
hc1     <- hclust(dist(trait.fa$loadings))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- trait.fa$loadings[ord.hc1,]
round(loading, 3)

# Now plot
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

#smat <- trait.fa$scores
#zlims <- max(abs(smat)) * c(-1,1)
#heatmap(smat, col=col, scale="none", zlim=zlims, labRow = NA, margins=c(7,0))


