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



# Functions ----------------------------------------------------------------

load.nn.openface <- function(vnames=NULL, parallel=T) {
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

basedir      <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)

openface     <- load.nn.openface(vnames=shape.vnames)
all.equal(shape.vnames, as.character(openface$labs$vid[istarts]))

system.time(d <- Rfast::Dist(openface$feats))
dmat <- as.matrix(d)

# Load timing information
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
dat.vols$sub01$basics$timing



# Main Timings --------------------------------------------------------------

# Get mean diff and frame diff
istarts <- which(openface$labs$frame==3)
mean.diffs <- sapply(istarts, function(ii) {
  mean(dmat[ii:(ii+8-1),-c(ii:(ii+8-1))]) # ave diff with everything else
})
mean.frame.diffs <- sapply(istarts, function(ii) {
  mean(sapply(1:7, function(si) dmat[ii+si-1,ii+si+1-1])) # ave framewise distances
})


# Subject Timings -----------------------------------------------------------

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
  inds   <- sapply(tvids, function(vname) which(vname==shape.vnames))
  all.equal(shape.vnames[inds], tvids)
  
  # Save the timings
  lines  <- afni.timing.amp(timing$onset, mean.diffs[inds], timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_mean_diff.txt")
  lines  <- afni.timing.amp(timing$onset, mean.frame.diffs[inds], timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_mean_frame_diff.txt")
}






# Play ----------------------------------------------------------------

# Cluster
library(ClusterR)
opt <- Optimal_Clusters_KMeans(openface$feats, 20, criterion = "silhouette", 
                               num_init = 20, max_iters = 200, 
                               initializer="kmeans++", threads=12)
km_rc <- KMeans_rcpp(openface$feats, clusters = 6, num_init = 20, max_iters = 200, 
                     initializer = 'kmeans++', threads = 12, verbose = F)
dim(km_rc$centroids)
table(km_rc$clusters, openface$labs$frame)
## get mean diff for clusters
system.time(d2 <- Rfast::Dist(rbind(km_rc$centroids, openface$feats)))
dmat2 <- as.matrix(d2)
tmp <- dmat2[-c(1:6),1:6]
tmp2 <- t(sapply(istarts, function(ii) colMeans(tmp[ii:(ii+8-1),])))
plot.ts(tmp2[1:10,], plot.type = "single", col=1:6)

# Some visuals
hist(mean.diffs)
hist(mean.frame.diffs)
cor(mean.diffs, mean.frame.diffs) # these are somewhat related
