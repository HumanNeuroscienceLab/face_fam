# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(plyr)
suppressMessages(library(doMC))
registerDoMC(20)

# output
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings"

# subjects
subjects <- sprintf("sub%02i", 1:6)

# load the demographics



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
  labels <- read.csv(sprintf('%s/masked_labels_redo.csv', base), header=F)
  features <- read.csv(sprintf('%s/masked_reps_redo.csv', base), header=F)
  
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

load.nn.openface.gender <- function() {
  base <- "/data1/famface01/analysis/misc/openface/prototypes_gender"
  
  # Read in
  labels <- read.csv(sprintf('%s/labels.csv', base), header=F)
  features <- read.csv(sprintf('%s/reps.csv', base), header=F)
  
  # Create new labels df and add the video name
  vids <- basename(as.character(labels[,2]))
  vids <- sub("average_face_", "", sub(".png", "", vids))
  
  rownames(features) <- vids
  features <- t(features)
  features <- as.data.frame(features)
  
  return(features)
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

# Get the PCA shapes (mainly for the ordering but also as a comparison)
sym.shapes   <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes)

# Load the nnet features
openface0    <- load.nn.openface(vnames=shape.vnames)
openface     <- load.nn.openface.masked(vnames=shape.vnames)
ofgender     <- load.nn.openface.gender()
istarts <- which(openface$labs$frame==3)
all.equal(shape.vnames, as.character(openface$labs$vid[istarts]))

# Load demographics
base  <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
demos <- demos[,-1]
## video names
demo.vnames <- sub("_fr[0-9]{3}", "", demos$video)
## rearrange
oinds <- sapply(shape.vnames, function(x) which(demo.vnames == x))
all.equal(demo.vnames[oinds], shape.vnames)
demos <- demos[oinds,]
## fill out
fr.gender <- demos$gender[rep(1:nrow(demos), each=8)]

# Load timing information
dat.vols <- lapply(subjects, load.cur)
names(dat.vols) <- subjects
head(dat.vols$sub01$basics$timing)



# Distances ---------------------------------------------------------------

# note: some dmats are actually dvecs

# would need to 
system.time(avg.dmat <- Rfast::Dist(rbind(openface$avg.feat, openface$feats)))
avg.dmat <- avg.dmat[-1,1]

system.time(male.dmat <- Rfast::Dist(rbind(ofgender$male, openface$feats))[-1,1])
system.time(female.dmat <- Rfast::Dist(rbind(ofgender$female, openface$feats))[-1,1])

# make a summary 'female/maleness' measure
diff.dmat <- (female.dmat/male.dmat)/max(female.dmat/male.dmat)

# make a difference to gender measure
# this is the difference for each frame to it's gender prototype
gender.dmat <- vector("numeric", length(fr.gender))
gender.dmat[fr.gender=="Male"] <- male.dmat[fr.gender=="Male"] - mean(male.dmat[fr.gender=="Male"])
gender.dmat[fr.gender=="Female"] <- female.dmat[fr.gender=="Female"] - mean(female.dmat[fr.gender=="Female"])

# see how the prototypes are related
# and see how the different distances to the average are related
cor(t(rbind(average=openface$avg.feat, male=ofgender$male, female=ofgender$female)))
cor(cbind(avg.dmat, male.dmat, female.dmat, gender.dmat, diff.dmat))
## note: diff measure should obviously be related to both elements

# plot the relation between the female and male dmats
plot(female.dmat, male.dmat, col=scales::alpha(as.numeric(fr.gender), 0.5), 
     xlab="Distance to Female Prototype", ylab="Distance to Male Prototype", 
     pch=16)

# how well can we use the gender information to explain things...
fit <- glm(fr.gender ~ avg.dmat + male.dmat*female.dmat, family=binomial(link='logit'))
summary(fit)

system.time(dmat1 <- Rfast::Dist(openface0$feats))
system.time(dmat2 <- Rfast::Dist(openface$feats))


m <- lm(male.dmat ~ fr.gender)$resid
f <- lm(female.dmat ~ fr.gender)$resid
mat <- cbind(gender=as.numeric(fr.gender), avg.dmat, male.m=m*(fr.gender=='Male'), 
             male.f=m*(fr.gender=='Female'), female.m=f*(fr.gender=='Male'), 
             female.f=f*(fr.gender=='Female'))
round(cor(mat), 3)
pca <- prcomp(mat[,], retx = T)
rot <- pca$rotation
corrplot(rot, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)


# PCA/FactorAnalysis ------------------------------------------------------

library(corrplot)
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))

# Trying to see if some sort-of of dimensionality reduction can be helpful
# since the different measures are very correlated

mat <- cbind(gender=as.numeric(fr.gender), avg.dmat, male.dmat, female.dmat, gender.dmat, diff.dmat)
cmat <- scale(mat, center=T, scale=T)

## PCA
# no this confuses things

pca <- prcomp(mat[,], retx = T)
rot <- pca$rotation
corrplot(rot, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

pca <- prcomp(mat[,-c(5:6)], retx = T)
rot <- pca$rotation
corrplot(rot, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

pca <- prcomp(mat[,-c(1,3,4,6)], retx = T)
rot <- pca$rotation
corrplot(rot, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)
corrplot(cor(mat[,-c(3,4,6)], pca$x), tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

## Factor Analysis
# no this confuses things

library(psych)
fa.parallel(mat) # suggests 3 factors
vss(mat, ncol(mat))
# Get the factors
fac.res <- fa(mat, nfactors=2, residuals=T, rotate='varimax', fm='minres')
# Save the loadings and scores, rename the columns
fac.loads  <- unclass(fac.res$loadings)
fac.scores <- fac.res$scores
colnames(fac.loads) <- sprintf("factor%02i", 1:ncol(fac.loads))
colnames(fac.scores) <- sprintf("factor%02i", 1:ncol(fac.scores))
# Plot
hc1     <- hclust(dist(fac.loads))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- fac.loads[ord.hc1,]
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)



# Main Timings --------------------------------------------------------------

istarts <- which(openface$labs$frame==3)

# We need to now reduce the values for each individual frame
mean.diffs <- sapply(istarts, function(ii) {
  mean(avg.dmat[ii:(ii+8-1)])
})
mean.diffs.male <- sapply(istarts, function(ii) {
  mean(male.dmat[ii:(ii+8-1)])
})
mean.diffs.female <- sapply(istarts, function(ii) {
  mean(female.dmat[ii:(ii+8-1)])
})
mean.diffs.gender <- sapply(istarts, function(ii) {
  mean(gender.dmat[ii:(ii+8-1)])
})
onefr.gender <- sapply(istarts, function(ii) {
  fr.gender[ii]
})
onefr.gender <- as.numeric(onefr.gender) - 1.5 # male > female

cor(cbind(mean.diffs, mean.diffs.male, mean.diffs.female))
cor(cbind(mean.diffs, mean.diffs.gender))



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
  #lines  <- afni.timing.amp(timing$onset, mean.diffs[inds], timing$run, max(t#iming$run), center=T)
  #save.tofile(lines, "stimam_nnet_masked_mean_diff.txt")
  #
  #lines  <- afni.timing.amp(timing$onset, mean.diffs.male[inds], timing$run, max(t#iming$run), center=T)
  #save.tofile(lines, "stimam_nnet_masked_mean_diff_male.txt")
  #
  #lines  <- afni.timing.amp(timing$onset, mean.diffs.female[inds], timing$run, max(t#iming$run), center=T)
  #save.tofile(lines, "stimam_nnet_masked_mean_diff_female.txt")
  #
  #lines  <- afni.timing.amp(timing$onset, onefr.gender[inds], timing$run, max(t#iming$run), center=F)
  #save.tofile(lines, "stimam_gender_diff.txt")
  
  lines  <- afni.timing.amp(timing$onset, mean.diffs.gender[inds], timing$run, max(timing$run), center=T)
  save.tofile(lines, "stimam_nnet_masked_mean_diff_gender.txt")
}

