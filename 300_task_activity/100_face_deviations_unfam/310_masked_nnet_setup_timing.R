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

basedir      <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)

openface0    <- load.nn.openface(vnames=shape.vnames)
openface     <- load.nn.openface.masked(vnames=shape.vnames)
istarts <- which(openface$labs$frame==3)
all.equal(shape.vnames, as.character(openface$labs$vid[istarts]))

# Load timing information
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
head(dat.vols$sub01$basics$timing)



# Distances ---------------------------------------------------------------

# would need to 
system.time(avg.dmat <- Rfast::Dist(rbind(openface$avg.feat, openface$feats)))
avg.dmat <- avg.dmat[-1,1]

system.time(dmat1 <- Rfast::Dist(openface0$feats))
system.time(dmat2 <- Rfast::Dist(openface$feats))



# Checks ---------------------------------------------------------------

# First point to make is that the correlation of the original feats 
# (whole face) and the newer feats masked is very good
cs <- sapply(1:nrow(dmat1), function(i) cor(dmat1[i,], dmat2[i,]))
hist(cs, main="Orig vs Masked (frame-by-frame)")

# Relate the distances with correlations
avg.cmat <- cor(t(openface$avg.feat), t(openface$feats))[1,]
plot.ts(cbind(avg.dmat, avg.cmat)[1:100,])
cor(avg.dmat, avg.cmat)
sort(avg.dmat)[1:10]
order(avg.cmat)[1:10]



# Main Timings --------------------------------------------------------------

istarts <- which(openface$labs$frame==3)

# Get the mean diff to each frame
frame.diffs <- avg.dmat

# Get the average of the frame diffs
mean.diffs <- sapply(istarts, function(ii) {
  mean(frame.diffs[ii:(ii+8-1)])
})
hist(mean.diffs)
mean.diffs2 <- sapply(istarts, function(ii) {
  mean(avg.cmat[ii:(ii+8-1)])
})
mean.diffs3 <- sapply(istarts, function(ii) {
  mean(mah3[ii:(ii+8-1)])
})


# Get residual frame diffs (note: this is pretty similar)
mean.frame.diffs <- sapply(istarts, function(ii) {
  fds <- frame.diffs[ii:(ii+8-1)] - mean(frame.diffs[ii:(ii+8-1)])
  mean(fds)
})
hist(mean.frame.diffs)
hist(scale(mean.frame.diffs, center=T, scale=T))

# Get mean framewise distance
mean.frame.diffs2 <- sapply(istarts, function(ii) {
  mean(sapply(1:7, function(si) dmat2[ii+si-1,ii+si+1-1])) # ave framewise distances
})
hist(mean.frame.diffs2)

# relation
cor(mean.diffs, mean.frame.diffs2)



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
  inds   <- sapply(tvids, function(vname) which(vname==shape.vnames))
  all.equal(shape.vnames[inds], tvids)
  
  # Save the timings
  lines  <- afni.timing.amp(timing$onset, mean.diffs[inds], timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_masked_mean_diff.txt")
  lines  <- afni.timing.amp(timing$onset, mean.frame.diffs2[inds], timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_masked_mean_frame_diff.txt")
}


# For the more complicated stuff (rolling average etc)
library(bigmemory)
library(doMC)
registerDoMC(12)
for (subj in subjects) {
  cat(subj, "\n")
  
  # Output
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  timing <- dat.vols[[subj]]$basics$timing
  
  ## get the different onset times with amps
  diff.scan.times <- ldply(1:max(timing$run), function(ri) {
    rtiming  <- timing[timing$run==ri,]
    run.vids <- as.character(rtiming$video)
    
    # find when time between two face videos < 2.5 and so consecutive
    succ.inds <- which(diff(rtiming$local.onset) < 2.5)
    
    scan.times <- ldply(succ.inds, function(si) {
      icur.vid  <- which(openface$labs$vid==run.vids[si] & openface$labs$frame==45)
      inext.vid <- which(openface$labs$vid==run.vids[si+1] & openface$labs$frame==3)
      c(onset=rtiming$onset[si+1] - 0.25, diff=dmat2[icur.vid,inext.vid])
    })
    
    data.frame(run=ri, scan.times)
  })
  ## save
  lines <- afni.timing.amp(diff.scan.times$onset, diff.scan.times$diff, 
                           timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_masked_btw_vid_diff.txt")
  
  ## get difference from rolling average
  roll.scan.times <- ldply(1:max(timing$run), function(ri) {
    rtiming  <- timing[timing$run==ri,]
    run.vids <- as.character(rtiming$video)
    
    roll.avg <- big.matrix(1, ncol(openface$feats), init=0)
    
    scan.times <- ldply(1:length(run.vids), function(si) {
      vid <- run.vids[si]
      sind <- which(openface$labs$vid == vid & openface$labs$frame==3)
      eind <- which(openface$labs$vid == vid & openface$labs$frame==45)
      
      mean.feats   <- colMeans(openface$feats[sind:eind,])
      avg.diff     <- Rfast::dista(roll.avg[1,], mean.feats)
      roll.avg[1,] <- (roll.avg[1,]*(si-1) + mean.feats)/si
      
      c(onset=rtiming$onset[si], diff=avg.diff)
    })
    
    data.frame(run=ri, scan.times)
  }, .parallel=T)
  ## save
  lines <- afni.timing.amp(roll.scan.times$onset, roll.scan.times$diff, 
                           timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_masked_rollavg_diff.txt")
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


# Typicality

# Load demographics information
base <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
demo.vnames <- demos$video # typo
demo.vnames <- sub("_fr[0-9]{3}", "", demo.vnames)
# Reorder rows to match features
oinds       <- sapply(shape.vnames, function(x) which(demo.vnames == x))
all.equal(demo.vnames[oinds], shape.vnames)
df.demos    <- demos[oinds,-c(1:2)]

# Load trait information
base         <- "/data1/famface01/analysis/encoding/12_Features"
traits       <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
trait.vnames <- as.character(traits[,1])
traits       <- traits[,-1]
# Reorder
oinds       <- sapply(shape.vnames, function(x) which(trait.vnames == x))
all.equal(trait.vnames[oinds], shape.vnames)
df.traits   <- traits[oinds,]

# Factor analysis
trait.fa <- factanal(df.traits, factors = 5, rotation = "varimax", 
                     na.action = na.omit, scores="regression")

# Load the mean diff measures with shape/texture
mean.vid.vals <- read.csv("100_face_deviations_unfam/measures/z_mean_vid_vals.csv")
all.equal(as.character(mean.vid.vals$vids), shape.vnames)

# Correlate
cor(df.traits$typical, mean.diffs)
cor(df.traits, mean.frame.diffs2)

cor(df.traits, mean.diffs)
cor(df.demos[,c("age")], mean.diffs)
cor(df.demos[,c("age")], mean.diffs)
cor(trait.fa$scores, mean.diffs)
cor(trait.fa$scores, mean.frame.diffs2)

cor(df.traits, mean.vid.vals[,colnames(mean.vid.vals)!="vids"])
cor(df.traits$typical, mean.vid.vals[,colnames(mean.vid.vals)!="vids"]) # not very correlated with typical

round(cor(cbind(df.traits, age=df.demos$age), cbind(mean.face=mean.vid.vals$mean_face, mean.nnet.face=mean.diffs)), 3)
tmp <- trait.fa$scores
colnames(tmp) <- c("fac.unemotional", "fac.competent", 'fac.trustworthy', "fac.memorable", "fac.attractive")
round(cor(cbind(tmp, age=df.demos$age), cbind(mean.face=mean.vid.vals$mean_face, mean.nnet.face=mean.diffs, texture=mean.vid.vals$pca_texture)), 3)

cor(mean.diffs, mean.vid.vals$mean_face)
cor(mean.frame.diffs2, mean.vid.vals$mean_fds)



# Outliers ----------------------------------------------------------------

# Using the nnet data
frame.diffs1 <- sapply(istarts, function(ii) {
  sapply(1:8, function(si) mean(dmat1[ii:(ii+8-1),ii:(ii+8-1)][si,-si])) # ave framewise distances
})
hist(apply(frame.diffs1, 2, max))
hist(apply(frame.diffs1, 2, mean))
bad.inds1 <- apply(frame.diffs1, 2, max) > quantile(apply(frame.diffs1, 2, max), p=0.95)
frame.diffs1[,bad.inds1]
shape.vnames[bad.inds1]


frame.diffs2 <- sapply(istarts, function(ii) {
  sapply(1:7, function(si) dmat1[ii+si-1,ii+si+1-1]) # ave framewise distances
})
frame.diffs2[,3]
hist(apply(frame.diffs2, 2, max))
hist(apply(frame.diffs2, 2, mean))
bad.inds <- apply(frame.diffs2, 2, max) > quantile(apply(frame.diffs2, 2, max), p=0.95)
frame.diffs2[,bad.inds]
shape.vnames[bad.inds]


# Use the face shapes




# Mahalanabois distance doesn't work well I think

cholMaha <- function(X) {
  dec <- chol( cov(X) )
  tmp <- forwardsolve(t(dec), t(X) )
  Rfast::Dist(t(tmp))
}

x <- as.matrix(openface$feats)

system.time(Sx <- cov(x))
#system.time(Sx <- Rfast::cova(x))

#system.time(mah1 <- mahalanobis(x, colMeans(x), Sx)) # slowest
#system.time(mah2 <- Rfast::mahala(x, colMeans(x), Sx)) # slower

system.time(mah1 <- mvnfast::maha(x, colMeans(x), Sx, ncores=4))
system.time(mah2 <- mvnfast::maha(x, as.numeric(openface$avg.feat), Sx, ncores=4))

system.time(mahpw <- cholMaha(rbind(openface$avg.feat,x)))
mah3 <- mahpw[1,-1]
head(tmp)
head(mah)
