# This script will generate the timing files for comparing the data to the average
# the average will be both the person and the whole
# We also want to make regressors for the likeness rating
# (can run a separate analysis with them in the model)

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



# Functions ----------------------------------------------------------------

load.nn.full <- function() {
  base <- "/data1/famface01/analysis/misc/openface/fampics"
  
  # Read in
  labels <- read.csv(sprintf('%s/masked_fullset_labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/masked_fullset_reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$img <- basename(as.character(labels.df$img))
  labels.df$img <- sub(".png", "", labels.df$img)
  labels.df$img <- factor(labels.df$img)
  
  # Add the person index
  persons <- sub("_[a-z]+_[0-9]*", "", labels.df$img)
  labels.df$person <- persons
  
  # Add original index
  labels.df$X <- 1:nrow(labels.df)
  
  # reorder by persons
  mat.feats <- as.matrix(features[order(labels.df$person),])
  labels.df <- labels.df[order(labels.df$person),]
  labels.df$X2 <- 1:nrow(labels.df)
  
  return(list(labs=labels.df, feats=mat.feats))
}

load.nn.ave <- function() {
  base <- "/data1/famface01/analysis/misc/openface/fampics"
  
  # Read in
  labels <- read.csv(sprintf('%s/masked_ave_labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/masked_ave_reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$img <- basename(as.character(labels.df$img))
  labels.df$img <- sub(".png", "", labels.df$img)
  labels.df$img <- factor(labels.df$img)
  
  # Add original index
  labels.df$X <- 1:nrow(labels.df)
  
  # reorder by persons
  mat.feats <- as.matrix(features[order(labels.df$img),])
  labels.df <- labels.df[order(labels.df$img),]
  labels.df$X2 <- 1:nrow(labels.df)
  
  return(list(labs=labels.df, feats=mat.feats))
}

load.nn.exp <- function() {
  base <- "/data1/famface01/analysis/misc/openface/fampics"
  
  # Read in
  labels <- read.csv(sprintf('%s/masked_exp_labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/masked_exp_reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  labels.df <- data.frame(inds=1:nrow(labels), img=as.character(labels[,2]))
  labels.df$img <- basename(as.character(labels.df$img))
  labels.df$img <- sub(".png", "", labels.df$img)
  labels.df$img <- factor(labels.df$img)
  
  # Add the person index
  persons <- gsub("_[0-9]+", "", labels.df$img)
  labels.df$person <- persons
  
  # Add original index
  labels.df$X <- 1:nrow(labels.df)
  
  # reorder by persons
  mat.feats <- as.matrix(features[order(labels.df$img),])
  labels.df <- labels.df[order(labels.df$img),]
  labels.df$X2 <- 1:nrow(labels.df)
  
  return(list(labs=labels.df, feats=mat.feats))
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

cholMaha <- function(X) {
  dec <- chol( cov(X) )
  tmp <- forwardsolve(t(dec), t(X) )
  Rfast::Dist(t(tmp))
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
  timing.df2$subj <- subj
  timing.df3 <- subset(timing.df2, select=c("subj", "tot.run", "name", "onset", "onset2_raw", "trial", "num", "fname", "stimtype", "dur", "rt", "target", "noresp", "incorrect"))
  
  cat("...incorrect\n")
  print(table(timing.df3$incorrect))
  cat("...no responses\n")
  print(table(timing.df3$noresp))
  
  return(timing.df3)
}



# Fun AFNI ----------------------------------------------------------------

afni.timing2 <- function(onsets, runs, nruns=8) {
  lines <- sapply(1:nruns, function(i) {
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

afni.timing.amp2 <- function(onsets, runs, amps, nruns=8) {
  lines <- sapply(1:nruns, function(i) {
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

afni.timing.dur2 <- function(onsets, runs, durs, nruns=8) {
  lines <- sapply(1:nruns, function(i) {
    inds <- runs == i
    if (any(inds)) {
      onset<- onsets[inds]
      dur  <- durs[inds]
      line <- paste(sprintf("%.5f:%.8f", onset, dur), collapse=" ")
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



# Load Likeness -------------------------------------------------------------

fpaths <- list.files("/data1/famface01/data/behavior/ratings", full.names = T)

# for two participants, we had them repeat the likeness rating
rep.fpaths <- c("/data1/famface01/data/behavior/ratings/fampics_likeness_sub2_rep1_setA_t6ks60wm_datetime_116-8-4_15:3:50.csv", "/data1/famface01/data/behavior/ratings/fampics_likeness_sub2_rep1_setB_qhde68x2_datetime_116-8-5_11:43:9.csv", "/data1/famface01/data/behavior/ratings/fampics_likeness_sub3_rep1_setA_lbgo79jj_datetime_116-8-2_22:57:3.csv", "/data1/famface01/data/behavior/ratings/fampics_likeness_sub3_rep1_setB_rkg3jgly_datetime_116-9-0_12:58:17.csv")

dfs <- ldply(fpaths, function(fpath) {
  raw.df <- read.csv(fpath, 1)
  df <- subset(raw.df, trial_type=="likert-simple" & person != "")
  keep.cols <- c("subid", "set", "person", "button_pressed", "rt", "pic_path")
  df <- df[,keep.cols]
  df$pic <- sub(".jpg", "", basename(as.character(df$pic_path)))
  df$pic_path <- NULL
  if (fpath %in% rep.fpaths) {
    df$rep <- 2
  } else {
    df$rep <- 1
  }
  df
}, .progress="text")

# note range is 0-6 or really 1-7
cdata <- ddply(dfs, c("person", "pic"), summarise,
               N    = length(button_pressed), 
               mean = mean(button_pressed),
               sd   = sd(button_pressed),
               se   = sd / sqrt(N)
)
cdata


# Load RT ------------------------------------------------------------------

# Load RT information
timing.df <- ldply(subjects, load.timing)
timing.df$pic <- sub(".jpg", "", basename(as.character(timing.df$fname)))
head(timing.df)
timing.df2 <- subset(timing.df, stimtype=='face')
rt.df <- ddply(timing.df2, .(name, pic), summarise, mean = mean(rt, na.rm=T))

# Select cdata
cinds <- sapply(cdata$pic, function(x) x %in% rt.df$pic)
cdata2 <- cdata[cinds,]
cdata2$pic <- factor(cdata2$pic)
cdata2$person <- factor(cdata2$person)

# Now select/order rt data
rinds <- sapply(cdata2$pic, function(x) which(x==rt.df$pic))
all.equal(as.character(cdata2$pic), as.character(rt.df$pic[rinds]))
rt.df2 <- rt.df[rinds,]



# Load Neural Net ---------------------------------------------------------

basedir      <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)

nn.ave   <- load.nn.ave()
nn.exp   <- load.nn.exp()
nn.full  <- load.nn.full()
nn.unfam <- load.nn.openface.masked(shape.vnames)

ninds <- sapply(cdata$pic, function(x) which(x==nn.exp$labs$img))
feats <- nn.exp$feats[ninds,]
labs  <- nn.exp$labs[ninds,]

ninds2 <- sapply(as.character(cdata2$pic), function(x) which(x==nn.exp$labs$img))
feats2 <- nn.exp$feats[ninds2,]
labs2  <- nn.exp$labs[ninds2,]

feats3 <- nn.exp$feats[-ninds,]
labs3  <- nn.exp$labs[-ninds,]


# Measures --------------------------------------------------------------

## Includes 2 additional photos

# Get each photo to the global average
avg1.dmat1 <- Rfast::dista(as.numeric(nn.unfam$avg.feat), t(feats))

# Get each photo to each person's average
avg2.dmat1 <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- Rfast::dista(ave.feat, cur.feat)
  d
})

cor(avg1.dmat1, avg2.dmat1)
summary(aov(cdata$mean ~ cdata$person + avg1.dmat1*avg2.dmat1))


## Removes those 2 additional photos and what was done in scanner

avg1.dmat2 <- Rfast::dista(as.numeric(nn.unfam$avg.feat), t(feats2))
avg2.dmat2 <- sapply(1:nrow(feats2), function(i) {
  # distance of average to each of these pics
  ave.feat <- nn.ave$feats[labs2$person[i] == nn.ave$labs$img,]
  cur.feat <- feats2[i,]
  d <- Rfast::dista(ave.feat, cur.feat)
  d
})

cor(avg1.dmat2, avg2.dmat2)
summary(aov(cdata2$mean ~ cdata2$person + avg1.dmat2*avg2.dmat2))
summary(aov(cdata2$mean ~ cdata2$person + avg2.dmat2))


## Get distances of distractor photos

avg1.dmat3 <- Rfast::dista(as.numeric(nn.unfam$avg.feat), t(feats3))



# AFNI Timing -------------------------------------------------------------

# Remember the more basic timing stuff has been done in another file

for (subj in subjects) {
  cat(subj, "\n")
  
  # Subset
  timing.df3 <- timing.df[timing.df$subj == subj,]
  
  # Output Setup
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  ## Distractors
  
  # Subset
  inds0        <- timing.df3$stimtype=="face" & timing.df3$target==0
  timing.imgs0 <- sub(".jpg", "", basename(as.character(timing.df3$fname[inds0])))
  imgs0        <- as.character(labs3$img)
  
  ds.inds0 <- sapply(timing.imgs0, function(iname) which(imgs0==iname))
  if (!all.equal(imgs0[ds.inds0], timing.imgs0)) stop("does not add up")
  
  # Faces
  lines <- afni.timing2(timing.df3$onset2_raw[inds0], timing.df3$tot.run[inds0])
  save.tofile(lines, "stim_faces_distractor.txt")
  
  # Avg distance to norm
  amps <- scale(avg1.dmat3[ds.inds0], center=T, scale=F)
  lines <- afni.timing.amp2(timing.df3$onset2_raw[inds0], timing.df3$tot.run[inds0], amps)
  save.tofile(lines, "stimam_nnet_avg_deviation_distractor.txt")
  
  ## Targets
  
  # Subset
  inds        <- timing.df3$stimtype=="face" & timing.df3$target==1
  timing.imgs <- sub(".jpg", "", basename(as.character(timing.df3$fname[inds])))
  imgs        <- as.character(labs2$img)
  
  ds.inds <- sapply(timing.imgs, function(iname) which(imgs==iname))
  if (!all.equal(imgs[ds.inds], timing.imgs)) stop("does not add up")
  
  # Save the whole sample deviation
  amps1 <- scale(avg1.dmat2[ds.inds], center=T, scale=F)
  lines <- afni.timing.amp2(timing.df3$onset2_raw[inds], timing.df3$tot.run[inds], amps1)
  save.tofile(lines, "stimam_nnet_avg_deviation_all.txt")
  
  # Save the person deviation
  amps2 <- scale(avg2.dmat2[ds.inds], center=T, scale=F)
  lines <- afni.timing.amp2(timing.df3$onset2_raw[inds], timing.df3$tot.run[inds], amps2)
  save.tofile(lines, "stimam_nnet_avg_deviation_person.txt")
  
  # Save each identity
  for (person in unique(as.character(cdata2$person))) {
    sinds <- timing.df3$stimtype=="face" & timing.df3$target==1 & timing.df3$name==person
    lines <- afni.timing2(timing.df3$onset2_raw[sinds], timing.df3$tot.run[sinds])
    save.tofile(lines, sprintf("stim_faces_%s.txt", person))
  }
  
  # Save the person deviation removing effect of identity
  resid <- lm(avg2.dmat2 ~ cdata2$person)$residuals # note: resids are centered
  lines <- afni.timing.amp2(timing.df3$onset2_raw[inds], timing.df3$tot.run[inds], resid)
  save.tofile(lines, "stimam_nnet_avg_deviation_person_resids.txt")
  
  # Save the likeness ratings
  ls.inds <- sapply(timing.imgs, function(iname) which(cdata2$pic==iname))
  if (!all.equal(ls.inds, ds.inds)) stop("does not add up")
  amps3 <- scale(cdata2$mean[ds.inds], center=T, scale=F)
  lines <- afni.timing.amp2(timing.df3$onset2_raw[inds], timing.df3$tot.run[inds], amps3)
  save.tofile(lines, "stimam_likeness_group.txt")
  
  ## NOTE: RTs saved elsewhere
}




# Check Xmat --------------------------------------------------------------

xmat_labs <- function(fn) {
  str <- system(sprintf("grep ColumnLabels %s | sed s/'#  ColumnLabels = '//", fn), intern=T)
  str <- gsub("\"", "", str)
  cols <- strsplit(str, ' ; ')[[1]]
  cols
}

xfile <- "/data1/famface01/analysis/task_activity/sub02/fampics/fampics_nnet_deviation.reml/xmat.1D"
xmat <- read.table(xfile)
colnames(xmat) <- xmat_labs(xfile)
xmat <- xmat[,c(9:14)]
colnames(xmat) <- sub("#0", "", colnames(xmat))
head(xmat)
round(cor(xmat), 3)
