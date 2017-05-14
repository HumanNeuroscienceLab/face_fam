# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(xlsx)
library(plyr)
suppressMessages(library(doMC))
registerDoMC(20)

# output
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings"

# subjects
subjects <- sprintf("sub%02i", 1:6)



# 4dface Cor----------------------------------------------------------------

# Note this stuff is repeated in the function below...

# I want to deal with the strong correlations between the different 4dface
# measures
# 
# This for now entails orthogonalizing some of the variables
# 
# Remove pose from everything
# 

base <- "/data1/famface01/analysis/misc/120_features"

df <- read.csv(sprintf("%s/df_meandiffs_4dface.csv", base))
df <- df[,-1]
mat <- as.matrix(df[,-1])
round(cor(mat), 3) # see the correlations

#mat.orth <- mat
#mat.orth[,-1] <- lm(mat[,-1] ~ df$pose)$residuals
#round(cor(mat.orth), 3)

## FACTOR ANALYSIS
# suggests 3 but we will use 4
# basically getting the absolute and framewise results
library(psych)
fa.parallel(mat) # suggests 3 factors
vss(mat, 7)
# Get the factors
fac.res <- fa(mat, nfactors=4, residuals=T, rotate='varimax', fm='minres')
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
corrplot(cor(raw.mat, fa.scores), tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

## PCA
# this actually shows each measure on it's own
# only issue is that the direction will need to be flipped on some to make
# dominant value positive
pca     <- prcomp(mat, retx = T)
# 
rot   <- apply(pca$rotation, 2, function(x) x * sign(x[which.max(abs(x))]))
pca.x <- sapply(1:ncol(mat), function(i) {
  s <- sign(pca$rotation[which.max(abs(pca$rotation[,i])),i])
  pca$x[,i] * s
})
comp.names <- colnames(mat)[apply(cor(mat, pca.x), 2, which.max)]
colnames(pca.x) <- sprintf("pca_%s", comp.names)
head(pca.x)
corrplot(rot, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)
corrplot(cor(raw.mat, pca.x), tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

## Raw
loading <- cor(mat)[c(1,3,5,7,2,4,6),c(1,3,5,7,2,4,6)]
corrplot(loading, tl.col='black', tl.cex=.75, diag=F, col=col, is.corr=T)

# save going forward
vidnames <- as.character(df$vid)
raw.mat <- mat
fa.mat  <- fac.scores
pca.mat <- pca.x
colnames(fa.mat) <- c("fac_dynamic", "fac_static", "fac_expression", "fac_texture")

colnames(raw.mat) <- sub("framewise", "fw", colnames(raw.mat))
colnames(pca.mat) <- sub("framewise", "fw", colnames(pca.mat))



# PCA Shape Measures ------------------------------------------------------

# Above we looked at shape with the whole face (and all the mesh points)
# These are the shape PCA from 4dface and ones I generated from the data

# 4dface shape PCA
df2 <- read.csv(sprintf("%s/df_meandiffs2_4dface.csv", base))
df2 <- df2[,-1]
mat2 <- as.matrix(df2[,-1])

# moi shape PCA
df3 <- read.csv(sprintf("%s/df_meandiffs3_4dface.csv", base))
df3 <- df3[,-1]
mat3 <- as.matrix(df3[,-1])

# Relation
# - Note: my shape stuff (df3) relate more to the pca.mat (original)
cor(df2[,c("shape", "framewise_shape")], pca.mat)
cor(df3[,c("shape", "framewise_shape")], pca.mat)


# PCA - 2 (4dface)
pca2    <- prcomp(mat2, retx = T)
## make sure comps positively related
rot2    <- apply(pca2$rotation, 2, function(x) x * sign(x[which.max(abs(x))]))
pca.x2  <- sapply(1:ncol(mat2), function(i) {
  s <- sign(pca2$rotation[which.max(abs(pca2$rotation[,i])),i])
  pca2$x[,i] * s
})
comp.names <- colnames(mat2)[apply(cor(mat2, pca.x2), 2, which.max)]
colnames(pca.x2) <- sprintf("pca_%s", comp.names)
## shape and expression here very related
corrplot(rot2, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)
corrplot(cor(raw.mat, pca.x2), tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)


# PCA - 3 (mine)
pca3    <- prcomp(mat3, retx = T)
## make sure comps positively related
rot3    <- apply(pca3$rotation, 2, function(x) x * sign(x[which.max(abs(x))]))
pca.x3  <- sapply(1:ncol(mat3), function(i) {
  s <- sign(pca3$rotation[which.max(abs(pca3$rotation[,i])),i])
  pca3$x[,i] * s
})
comp.names <- colnames(mat3)[apply(cor(mat3, pca.x3), 2, which.max)]
colnames(pca.x3) <- sprintf("pca_%s", comp.names)
## shape and expression here very related
corrplot(rot3, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)
corrplot(cor(raw.mat, pca.x3), tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

pca.mat2 <- pca.x2
pca.mat3 <- pca.x3



# Functions ----------------------------------------------------------------

load.nn.openface <- function(vnames=NULL) {
  base <- "/data1/famface01/analysis/encoding/12_Features"
  
  # Read in
  labels <- read.csv(sprintf('%s/labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$vid <- basename(as.character(labels.df$img))
  labels.df$vid <- sub(".png", "", labels.df$vid)
  labels.df$vid <- sub("_fr[0-9]{3}$", "", labels.df$vid)
  labels.df$vid <- factor(labels.df$vid)
  
  # Add the frame index
  frs <- sapply(1:nrow(labels.df), function(i) sub(as.character(labels.df$vid)[i], "", basename(as.character(labels.df$img)[i])))
  frs <- sub("_fr", "", frs)
  frs <- sub(".png", "", frs)
  frs <- as.numeric(frs)
  labels.df$frame <- frs
  
  # Remove vids of famous ppl
  inds <- sapply(as.character(labels.df$vid), function(vid) {
    any(vnames == vid)
  })
  labels.df <- labels.df[inds,]
  features  <- features[inds,]
  labels.df$vid <- factor(labels.df$vid)
  
  # Add index
  labels.df$X <- 1:nrow(labels.df)
  
  # reorder frames
  labels.df2 <- ddply(labels.df, .(vid), function(x) {
    x[order(x$frame),]
  })
  
  # Match up
  if (!is.null(vnames)) {
    oinds <- unlist(lapply(vnames, function(x) which(labels.df2$vid == x)))
    all.equal(as.character(labels.df2$vid[oinds[seq(1,length(oinds),by=8)]]), vnames)
    df.labels <- labels.df2[oinds,]
    mat.feats <- features[df.labels$X,]
  } else {
    df.labels <- labels.df2
    mat.feats <- features
  }
  
  return(list(labs=df.labels, feats=mat.feats))
}

load.nn.openface.masked <- function(vnames=NULL) {
  base <- "/data1/famface01/analysis/misc/openface"
  
  # Read in
  labels <- read.csv(sprintf('%s/masked_labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/masked_reps.csv', base), header=F)
  
  # Find the average face
  aind <- grep("average_face", labels$V2)
  avg.lab <- labels$V2[aind]
  avg.features <- features[aind,]
  ## remove the average
  labels <- labels[-aind,]
  features <- features[-aind,]
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$vid <- basename(as.character(labels.df$img))
  labels.df$vid <- sub(".png", "", labels.df$vid)
  labels.df$vid <- sub("_fr[0-9]{3}$", "", labels.df$vid)
  labels.df$vid <- factor(labels.df$vid)
  
  # Add the frame index
  frs <- sapply(1:nrow(labels.df), function(i) sub(as.character(labels.df$vid)[i], "", basename(as.character(labels.df$img)[i])))
  frs <- sub("_fr", "", frs)
  frs <- sub(".png", "", frs)
  frs <- as.numeric(frs)
  labels.df$frame <- frs
  
  # Remove vids of famous ppl
  inds <- sapply(as.character(labels.df$vid), function(vid) {
    any(vnames == vid)
  })
  labels.df <- labels.df[inds,]
  features  <- features[inds,]
  labels.df$vid <- factor(labels.df$vid)
  
  # Add index
  labels.df$X <- 1:nrow(labels.df)
  
  # reorder frames
  labels.df2 <- ddply(labels.df, .(vid), function(x) {
    x[order(x$frame),]
  })
  
  # Match up
  if (!is.null(vnames)) {
    oinds <- unlist(lapply(vnames, function(x) which(labels.df2$vid == x)))
    all.equal(as.character(labels.df2$vid[oinds[seq(1,length(oinds),by=8)]]), vnames)
    df.labels <- labels.df2[oinds,]
    mat.feats <- features[df.labels$X,]
  } else {
    df.labels <- labels.df2
    mat.feats <- features
  }
  
  return(list(labs=df.labels, feats=mat.feats, avg.lab=avg.lab, avg.feat=avg.features))
}

afni.timing <- function(onsets, runs, nruns) {
  lines <- sapply(1:nruns, function(i) {
    inds <- runs == i
    if (any(inds)) {
      x    <- onsets[inds]
      line <- paste(as.character(round(x, 5)), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  as.character(lines)
}

afni.timing.amp <- function(onsets, amps, runs, nruns, center=F, scale=F) {
  amps <- scale(amps, center=center, scale=scale)
  lines <- sapply(1:nruns, function(i) {
    inds <- runs == i
    if (any(inds)) {
      amp  <- amps[inds]
      onset<- onsets[inds]
      line <- paste(sprintf("%.5f*%.8f", onset, amp), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  
  as.character(lines)
}

load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}


# Load ---------------------------------------------------------------------

# For some reason, I'm first matching the order to shape.vnames
# but later it will be changed to the order of the individual subject
basedir      <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)

#
openface0    <- load.nn.openface(vnames=shape.vnames)
openface     <- load.nn.openface.masked(vnames=shape.vnames)
istarts <- which(openface$labs$frame==3)
all.equal(shape.vnames, as.character(openface$labs$vid[istarts]))

# Load timing information
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
head(dat.vols$sub01$basics$timing)


# Distances ---------------------------------------------------------------

istarts <- which(openface$labs$frame==3)

system.time(avg.dmat <- Rfast::Dist(rbind(openface$avg.feat, openface$feats)))
system.time(dmat1 <- Rfast::Dist(openface0$feats))
system.time(dmat2 <- Rfast::Dist(openface$feats))

# Get the mean diff to each frame
frame.diffs <- avg.dmat[-1,1]
mean.diffs <- sapply(istarts, function(ii) {
  mean(frame.diffs[ii:(ii+8-1)])
})

# Get mean framewise distance
mean.frame.diffs2 <- sapply(istarts, function(ii) {
  mean(sapply(1:7, function(si) dmat2[ii+si-1,ii+si+1-1])) # ave framewise distances
})

cor(mean.diffs, mean.frame.diffs2)

# relation
rinds <- sapply(shape.vnames, function(x) which(vidnames==x))
all.equal(shape.vnames, vidnames[rinds])
round(cor(raw.mat[rinds,], cbind(mean.diffs, mean.frame.diffs2)),3) # frame measures related
round(cor(pca.mat[rinds,], cbind(mean.diffs, mean.frame.diffs2)),3) # a little with frame measures
round(cor(fa.mat[rinds,], cbind(mean.diffs, mean.frame.diffs2)),3) # dynamic and expression and texture

round(cor(pca.mat2[rinds,], cbind(mean.diffs, mean.frame.diffs2)),3)
round(cor(pca.mat3[rinds,], cbind(mean.diffs, mean.frame.diffs2)),3)

#round(cor(fac.scores[rinds,], cbind(mean.diffs, mean.frame.diffs2)),3)

system.time(tmp3 <- apply(sym.shapes^2, 1, function(x) sqrt(mean(x))))
round(cor(tmp3, cbind(mean.diffs, mean.frame.diffs2)),3)


# Subject Timings -----------------------------------------------------------

# For the simple stuff
for (subj in subjects) {
  cat(subj, "\n")
  
  # Output
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  # Order the videos
  timing <- dat.vols[[subj]]$basics$timing
  tvids  <- as.character(timing$video)
  inds   <- sapply(tvids, function(vname) which(vname==vidnames))
  if (!all.equal(vidnames[inds], tvids)) stop("reorder didn't work")
  
  # Save the timings for the factor analysis and PCA
  for (i in 1:ncol(pca.mat)) {
    cname <- colnames(pca.mat)[i]
    x <- pca.mat[inds,i]
    lines  <- afni.timing.amp(timing$onset, x, timing$run, max(timing$run), center=T)
    save.tofile(lines, sprintf("stimam_4dface_%s.txt", cname))
  }
  for (i in 1:ncol(fa.mat)) {
    cname  <- colnames(fa.mat)[i]
    x      <- fa.mat[inds,i]
    lines  <- afni.timing.amp(timing$onset, x, timing$run, max(timing$run), center=T)
    save.tofile(lines, sprintf("stimam_4dface_%s.txt", cname))
  }
}


# For the extra shape stuff
for (subj in subjects) {
  cat(subj, "\n")
  
  # Output
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  # Order the videos
  timing <- dat.vols[[subj]]$basics$timing
  tvids  <- as.character(timing$video)
  inds   <- sapply(tvids, function(vname) which(vname==vidnames))
  if (!all.equal(vidnames[inds], tvids)) stop("reorder didn't work")
  
  # Save the timings for the factor analysis and PCA
  for (i in 1:ncol(pca.mat2)) {
    cname <- colnames(pca.mat2)[i]
    x <- pca.mat2[inds,i]
    lines  <- afni.timing.amp(timing$onset, x, timing$run, max(timing$run), center=T)
    save.tofile(lines, sprintf("stimam_4dface2_%s.txt", cname))
  }
  for (i in 1:ncol(pca.mat3)) {
    cname  <- colnames(pca.mat3)[i]
    x      <- pca.mat3[inds,i]
    lines  <- afni.timing.amp(timing$onset, x, timing$run, max(timing$run), center=T)
    save.tofile(lines, sprintf("stimam_4dface3_%s.txt", cname))
  }
}
