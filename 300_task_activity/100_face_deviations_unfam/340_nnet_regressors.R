# This script will generate the timing information for
# - pose, luminance
# - avg dist, and knn dist
# - framewise dist, framewise pose, framewise luminance

# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(plyr)
library(bigmemory)
suppressMessages(library(doMC))
registerDoMC(20)



# Load Video Characteristics ----------------------------------------------

# Note: we won't actually use all this info but I load it all anyway in cases

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

# Subset if you want
vdf2     <- df.info
vdf2$vid <- factor(vdf2$vid)



# Load Face Data ----------------------------------------------------------

# Get average first
load.aveface <- function() {
  base <- "/data1/famface01/analysis/misc/openface"
  
  # Read in
  labels <- read.csv(sprintf('%s/masked_labels_redo.csv', base), header=F)
  features <- read.csv(sprintf('%s/masked_reps_redo.csv', base), header=F)
  
  # Find the average face
  aind <- grep("average_face", labels$V2)
  avg.lab <- labels$V2[aind]
  avg.features <- features[aind,]
  
  row.names(avg.features) <- avg.lab
  avg.features
}
avg.feats <- load.aveface()


base <- "/data1/famface01/analysis/misc/openface"
# Read in
labels <- read.csv(sprintf('%s/masked_labels_redo.csv', base), header=F)
features <- read.csv(sprintf('%s/masked_reps_redo.csv', base), header=F)
ref.vnames <- sub(".mat", "", basename(as.character(vdf2$fn)))
match.vnames <- sub(".png", "", basename(as.character(labels$V2)))
oinds       <- sapply(1:nrow(vdf2), function(i) {
  which((match.vnames == ref.vnames[i]))
})
all.equal(ref.vnames, match.vnames[oinds])
fr.feats <- features[oinds,]
vid.feats <- apply(fr.feats, 2, tapply, vdf2$vid, mean)

# Load the last layer features
bdir <- "/data1/famface01/analysis/misc/openface/concatenate_layers"
bm   <- attach.big.matrix(file.path(bdir, sprintf("clayer%02i.desc", 26)))

# Extract the frame and video features
fr.feats2  <- as.matrix(bm)
vid.feats2 <- apply(fr.feats2, 2, tapply, vdf2$vid, mean)

## double check
#base <- "/data1/famface01/analysis/misc/openface"
## Read in
#labels <- read.csv(sprintf('%s/labels.csv', base), header=F)
#features <- read.csv(sprintf('%s/reps.csv', base), header=F)
#ref.vnames <- sub(".mat", "", basename(as.character(vdf2$fn)))
#match.vnames <- sub(".png", "", basename(as.character(labels$V2)))
#oinds       <- sapply(1:nrow(vdf2), function(i) {
#  which((match.vnames == ref.vnames[i]))
#})
#all.equal(ref.vnames, match.vnames[oinds])
#diag(cor(t(features[oinds,][1:10,]), t(fr.feats2[1:10,])))
### extract dists
#fr.feats3 <- features[oinds,]
#vid.feats3 <- apply(fr.feats3, 2, tapply, vdf2$vid, mean)


# Distances ---------------------------------------------------------------

# everything (masked)
vid.dmat <- Rfast::Dist(vid.feats)
fr.dmat  <- Rfast::Dist(fr.feats)

# everything (unmasked)
vid.dmat2 <- Rfast::Dist(vid.feats2)
fr.dmat2  <- Rfast::Dist(fr.feats2)

# to average
avg.dmat <- Rfast::Dist(rbind(avg.feats, vid.feats))[-1,1]
avg.dmat2 <- Rfast::Dist(rbind(colMeans(vid.feats2), vid.feats2))[-1,1]
avg.dmat3 <- Rfast::Dist(rbind(avg.feats, vid.feats2))[-1,1]

# to knn
knn.dmat <- sapply(1:nrow(vid.dmat), function(i) mean(sort(vid.dmat[i,-i])[1:10]))

## the two are fairly correlated
cor(cbind(avg.dmat, avg.dmat2, avg.dmat3, knn.dmat))



# Compare Dists with Traits -----------------------------------------------

# factor analysis
fac.res <- psych::fa(vdf2[vdf2$fr==3,18:ncol(vdf2)], nfactors=5, residuals=T, rotate='varimax', fm='minres')
fac.loads  <- unclass(fac.res$loadings)
fac.scores <- fac.res$scores
colnames(fac.loads) <- sprintf("factor%02i", 1:ncol(fac.loads))
colnames(fac.scores) <- sprintf("factor%02i", 1:ncol(fac.scores))

# relate
smat <- cbind(avg.dmat, avg.dmat3, knn.dmat)
summary(lm(fac.scores ~ avg.dmat + knn.dmat))
summary(lm(smat ~ fac.scores)) # see that avg.dmat might be better
summary(lm(smat ~ gender + age + glasses + makeup, data=vdf2[vdf2$fr==3,]))

summary(aov(smat ~ gender + age + pose + luminance, data=vdf2[vdf2$fr==3,]))
summary(aov(smat[,-2] ~ gender + age + pose + luminance, data=vdf2[vdf2$fr==3,]))
t.test(smat[,1] ~ vdf2$gender[vdf2$fr==3])
t.test(smat[,3] ~ vdf2$gender[vdf2$fr==3])


# Summarize for Timing ----------------------------------------------------

istarts <- which(vdf2$fr==3)

# low-level
luminance <- tapply(vdf2$luminance, vdf2$vid, mean)
pose <- tapply(vdf2$fd, vdf2$vid, mean)

# face
mean.diff <- avg.dmat3
knn.diff  <- knn.dmat

# dynamics
frame.luminance <- sapply(istarts, function(ii) {
  mean(diff(vdf2$luminance[ii:(ii+7)]))
})
frame.pose <- laply(istarts, function(ii) {
  mean(sapply(1:7, function(si) {
    Rfast::dista(vdf2[ii+si-1,5:7], vdf2[ii+si,5:7])
  }))
}, .parallel=T)
frame.mean.diff <- laply(istarts, function(ii) {
  mean(sapply(1:7, function(si) {
    fr.dmat2[ii+si-1,ii+si]
  }))
}, .parallel=T)

# compare
#mat <- cbind(luminance, pose, mean.diff, knn.diff, frame.luminance, frame.pose, frame.mean.diff)
mat <- cbind( pose, mean.diff, knn.diff, frame.pose, frame.mean.diff)
round(cor(mat), 3)



# PCA ---------------------------------------------------------------------

pca <- prcomp(mat, retx=T)
library(corrplot)
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))
corrplot(pca$rotation, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

# use the pca loadings too, make sure sign of max loading is positive
rot <- apply(pca$rotation, 2, function(x) x * sign(x[which.max(abs(x))]))
corrplot(rot, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)
pca.mat <- sapply(1:ncol(rot), function(i) pca$x[,i] * sign(pca$rotation[which.max(abs(pca$rotation[,i])),i]))
colnames(pca.mat) <- paste("pca", rownames(rot)[apply(rot, 2, which.max)], sep=".")
corrplot(cor(mat, pca.mat), tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)



# Timing ------------------------------------------------------------------

afni.timing <- function(onsets, runs, nruns=16) {
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

afni.timing.amp <- function(onsets, amps, runs, nruns=16, center=T, scale=F) {
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

# Load timing information
subjects <- sprintf("sub%02i", 1:6)
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
head(dat.vols$sub01$basics$timing)

# For the basic timing
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings"
vnames <- rownames(mat)
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
  inds   <- sapply(tvids, function(vname) which(vname==vnames))
  if (!all.equal(vnames[inds], tvids)) stop("names unequal")
  
  # Save the timings
  for (i in 1:ncol(mat)) {
    x     <- mat[,i]
    xname <- gsub("[.]", "_", colnames(mat)[i])
    lines  <- afni.timing.amp(timing$onset, x[inds], timing$run)
    save.tofile(lines, sprintf("stimam_nnet2_%s.txt", xname))
  }
}

# For the PCA timing
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings"
vnames <- rownames(mat)
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
  inds   <- sapply(tvids, function(vname) which(vname==vnames))
  if (!all.equal(vnames[inds], tvids)) stop("names unequal")
  
  # Save the timings
  for (i in 1:ncol(pca.mat)) {
    x     <- pca.mat[,i]
    xname <- gsub("[.]", "_", colnames(pca.mat)[i])
    lines  <- afni.timing.amp(timing$onset, x[inds], timing$run)
    save.tofile(lines, sprintf("stimam_nnet2_%s.txt", xname))
  }
}
