
# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(xlsx)
library(plyr)
suppressMessages(library(doMC))
registerDoMC(20)

# output
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/300_fampics/timings"

# subjects
subjects <- sprintf("sub%02i", 1:6)


# Load --------------------------------------------------------------------


load.nn.openface <- function(vnames=NULL) {
  base <- "/data1/famface01/analysis/misc/openface"
  
  # Read in
  labels <- read.csv(sprintf('%s/fampics_labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/fampics_reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$vid <- basename(as.character(labels.df$img))
  labels.df$vid <- sub(".png", "", labels.df$vid)
  labels.df$vid <- sub("_fr[0-9]{3}$", "", labels.df$vid)
  labels.df$vid <- factor(labels.df$vid)
  
  # Add the identity
  ids <- gsub("_[0-9]+", "", labels.df$vid)
  labels.df$person <- ids
  
  # Add index
  labels.df$X <- 1:nrow(labels.df)
  
  # Match up
  if (!is.null(vnames)) {
    oinds <- unlist(lapply(vnames, function(x) which(labels.df$vid == x)))
    all.equal(as.character(labels.df$vid[oinds[seq(1,length(oinds),by=8)]]), vnames)
    df.labels <- labels.df[oinds,]
    mat.feats <- features[df.labels$X,]
  } else {
    df.labels <- labels.df
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

afni.timing <- function(tdf) {
  lines <- sapply(1:8, function(i) {
    inds <- tdf$tot.run == i
    if (any(inds)) {
      x    <- tdf[inds,,drop=F]
      line <- paste(as.character(round(x$onset2_raw, 5)), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  as.character(lines)
}

afni.timing2 <- function(onsets, runs) {
  lines <- sapply(1:8, function(i) {
    inds <- runs == i
    if (any(inds)) {
      onset <- onsets[inds]
      line <- paste(as.character(round(onset, 5)), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  as.character(lines)
}

afni.timing.amp <- function(tdf, amps) {
  if (nrow(tdf) != length(amps)) stop("tdf and amps lengths differ")
  
  lines <- sapply(1:8, function(i) {
    inds <- tdf$tot.run == i
    if (any(inds)) {
      amp  <- amps[inds]
      x    <- tdf[inds,,drop=F]
      line <- paste(sprintf("%.5f*%.8f", x$onset2_raw, amp), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  
  as.character(lines)
}

afni.timing.amp2 <- function(onsets, runs, amps) {
  lines <- sapply(1:8, function(i) {
    inds <- runs == i
    if (any(inds)) {
      onset<- onsets[inds]
      amp  <- amps[inds]
      line <- paste(sprintf("%.5f*%.8f", onset, amp), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  
  as.character(lines)
}

afni.timing.dur <- function(tdf, durs) {
  if (nrow(tdf) != length(durs)) stop("tdf and amps lengths differ")
  
  lines <- sapply(1:8, function(i) {
    inds <- tdf$tot.run == i
    if (any(inds)) {
      dur  <- durs[inds]
      x    <- tdf[inds,,drop=F]
      line <- paste(sprintf("%.5f:%.8f", x$onset2_raw, dur), collapse=" ")
    } else {
      line <- "*"
    }
    line
  })
  
  as.character(lines)
}

load.mc <- function(subj) {
  funcdir  <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  df.paths <- read.table(file.path(funcdir, "df_paths.txt"), header=T)
  inds     <- df.paths$inindex[df.paths$name == "fam_pics"]
  fpaths   <- sprintf("%s/mc/func_run%02i_dfile.1D", funcdir, inds)
  
  mcs <- ldply(fpaths, function(fpath) {
    x <- read.table(fpath)
    x <- as.matrix(x)
    x <- scale(x, scale=F, center=T)
    x
  })
  mcs <- as.matrix(mcs)
  
  mcs
}


load.timing <- function(subj) {
  cat(subj, "\n")
  infiles <- list.files(sprintf("/data1/famface01/data/behavior/%s/fam_pics", subj), 
                        full.names = T)
  timing.df <- ldply(infiles, function(infile) {
    timing <- read.xlsx(infile, 1)
    timing <- timing[!is.na(timing$onset2_raw),]
    timing$sess <- as.numeric(as.character(timing$sess))
    timing
  }, .progress="text")
  
  targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
               "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
  
  # Fix for sub01
  if (subj == "sub01") {
    timing.df$run[108:214] <- 2
  }
  
  # Add the actual (total) run number
  timing.df$tot.run <- (timing.df$sess-1)*2 + timing.df$run
  if (length(unique(timing.df$tot.run)) != 8) stop("run error")
  
  # Is the trial of a target face
  timing.df$target <- (timing.df$name %in% targets)*1
  
  # What are the button presses for each target
  button.presses <- with(subset(timing.df, name %in% targets), table(factor(name), response_raw))
  #button.presses <- with(subset(timing.df, name %in% targets), table(factor(name), response_raw, sess))
  print(button.presses)
  opts <- gsub("'", "", colnames(button.presses)[-1])
  right.buttons <- opts[apply(button.presses[,-1], 1, which.max)]
  names(right.buttons) <- rownames(button.presses)
  right.buttons
  
  # Get no responses
  inds <- timing.df$name %in% targets
  timing.df$noresp <- 0
  timing.df$noresp[inds] <- (as.character(timing.df$response_raw[inds]) == "'--'") * 1
  
  # Indicate if button press is correct
  inds <- timing.df$name %in% targets
  timing.df2 <- ddply(timing.df, .(name), function(sdf) {
    if (sdf$name[1] %in% targets) {
      butt <- right.buttons[names(right.buttons) == as.character(sdf$name[1])]
      bad  <- (gsub("'", "", as.character(sdf$response_raw)) != butt) & (sdf$noresp!=1)
      sdf$incorrect <- bad*1
    } else {
      sdf$incorrect <- 0
    }
    sdf
  })
  ## fix for 1st subject
  if (subj == "sub01") {
    timing.df2$incorrect[timing.df2$name == "Angelina_Jolie"] <- 0
    timing.df2$incorrect[timing.df2$name == "Justin_Timberlake"] <- 0
    timing.df2$incorrect[timing.df2$name == "Will_Smith"] <- 0
    timing.df2$incorrect[timing.df2$name == "Julia_Roberts"] <- 0
    timing.df2$incorrect[timing.df2$name == "Julia_Roberts" & timing.df2$response_raw %in% c("'b'", "'g'")] <- 1
  }
  timing.df2 <- timing.df2[order(timing.df2$tot.run, timing.df2$order),]
  
  # set RT relative to onset
  raw_rts  <- as.numeric(as.character(timing.df2$RT_raw))
  inds     <- !is.na(raw_rts)
  timing.df2$rt <- NA
  timing.df2$rt[inds] <- raw_rts[inds] - timing.df2$onset2_raw[inds]
  
  # Reformat everything
  timing.df3 <- subset(timing.df2, select=c("tot.run", "name", "onset", "onset2_raw", "trial", "num", "fname", "stimtype", "dur", "rt", "target", "noresp", "incorrect"))
  
  cat("...incorrect\n")
  print(table(timing.df3$incorrect))
  cat("...no responses\n")
  print(table(timing.df3$noresp))
  
  return(timing.df3)
}


# Load Nnet Features -------------------------------------------------------

# Get order
basedir      <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)

# Get the famous face pics
openface.fpics    <- load.nn.openface(vnames=NULL)
## test to see if can automatically detect the identity
library(dynamicTreeCut)
dmat <- Rfast::Dist(openface.fpics$feats)
d <- as.dist(dmat)
hc <- hclust(d, method="ward.D2")
cl <- cutreeDynamic(hc, minClusterSize=2, distM=dmat)
table(openface.fpics$labs$person, cl) # whoops showed same sandra bullock photo many times

# Get masked unfam faces (has the average face)
openface     <- load.nn.openface.masked(vnames=shape.vnames)
istarts <- which(openface$labs$frame==3)
all.equal(shape.vnames, as.character(openface$labs$vid[istarts]))

# Load timing information
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
head(dat.vols$sub01$basics$timing)



library(dynamicTreeCut)
dmat0 <- Rfast::Dist(t(openface.fpics$feats))
d <- as.dist(dmat0)
hc <- hclust(d, method="ward.D2")
cl <- cutreeDynamic(hc, minClusterSize=2, distM=dmat0)
table(cl)



# Distances ---------------------------------------------------------------

# Compare each famous picture to the average face from the videos
system.time(avg.dmat <- Rfast::Dist(rbind(openface$avg.feat, openface.fpics$feats)))
avg.dmat <- avg.dmat[-1,1]

system.time(dmat1 <- Rfast::Dist(openface.fpics$feats))
system.time(dmat2 <- Rfast::Dist(openface$feats))



# Timing Files --------------------------------------------------------------

# Collect all the timings first
timings <- llply(subjects, load.timing)
names(timings) <- subjects
# REMEMBER THESE TIMES ARE LOCAL

### BASICS
#for (subj in subjects) {
#  cat(subj, "\n")
#  
#  timing.df3 <- timings[[subj]]
#  
#  sub.outdir <- file.path(outdir, subj)
#  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
#  
#  save.tofile <- function(lines, ofile) {
#    writeLines(lines, file.path(sub.outdir, ofile))
#  }
#  
#  # Face 
#  inds <- timing.df3$stimtype=="face"
#  runs <- timing.df3$tot.run[inds]
#  onsets <- timing.df3$onset2_raw[inds]
#  lines <- afni.timing(timing.df3[inds,])
#  save.tofile(lines, "stim_faces.txt")
#  
#  # Incorrect
#  inds <- timing.df3$incorrect==1
#  lines <- afni.timing(timing.df3[inds,])
#  save.tofile(lines, "stim_incorrect.txt")
#  
#  # No Response
#  inds <- timing.df3$noresp==1
#  lines <- afni.timing(timing.df3[inds,])
#  save.tofile(lines, "stim_noresp.txt")
#  
#  # Button Press
#  rts  <- timing.df3$rt
#  inds <- !is.na(rts) & timing.df3$stimtype=="face"
#  bp.onset <- timing.df3$onset2_raw[inds] + rts[inds]
#  lines <- afni.timing2(timing.df3[inds,], bp.onset)
#  save.tofile(lines, "stim_buttonpress.txt")
#  
#  # RT (amp)
#  rts  <- timing.df3$rt
#  inds <- !is.na(rts)
#  amps <- scale(rts[inds], center=T, scale=F)
#  lines <- afni.timing.amp(timing.df3[inds,], amps)
#  save.tofile(lines, "stimam_rt.txt")
#  
#  # RT (dur)
#  rts  <- timing.df3$rt
#  inds <- !is.na(rts)
#  durs <- rts[inds]
#  lines <- afni.timing.dur(timing.df3[inds,], durs)
#  save.tofile(lines, "stimdur_rt.txt")
#  
#  # Motion
#  mcs <- load.mc(subj)
#  ofile <- file.path(sub.outdir, "motion.1D")
#  write.table(mcs, file=ofile, row.names=F, col.names=F, quote=F)
#}


## MORE ADVANCED STUFF
for (subj in subjects) {
  cat(subj, "\n")

  timing.df1 <- timings[[subj]]
  timing.df3 <- timing.df1
  
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  inds <- timing.df3$stimtype=="face"
  timing.vids   <- sub(".jpg", "", basename(as.character(timing.df3$fname[inds])))
  openface.vids <- as.character(openface.fpics$labs$vid)
  
  ds.inds <- sapply(timing.vids, function(vname) which(openface.vids==vname))
  all.equal(openface.vids[ds.inds], timing.vids)
  amps <- scale(avg.dmat[ds.inds], center=T, scale=F)
  lines <- afni.timing.amp2(timing.df3$onset2_raw[inds], timing.df3$tot.run[inds], amps)
  save.tofile(lines, "stimam_nnet_avg_deviation.txt")
}

