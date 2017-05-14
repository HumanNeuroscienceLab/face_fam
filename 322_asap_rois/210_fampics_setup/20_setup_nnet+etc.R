library(plyr)
library(bigmemory)

# Load nnet ----------------------------------------------------------

# NOTE: THE MASKED DATA IS ACTUALLY NOT MASKED (FOR NOW)

targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
             "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")

# Load the experiment faces
load.exp <- function() {
  ## read in
  labs <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_exp_labels.csv", header=F)
  reps <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_exp_reps.csv", header=F)
  
  ## format
  ifiles <- basename(as.character(labs$V2))
  inames <- sub(".png", "", ifiles)
  pnames <- sub("_[0-9].*", "", inames)
  idf <- data.frame(ind=1, oind=1:nrow(reps), 
                    target=pnames %in% targets, 
                    person=pnames, pic=inames, fn=ifiles)
  idf <- idf[order(idf$pic),]
  idf$ind <- 1:nrow(idf)
  
  ## reorder
  ifeats <- reps[idf$oind,]
  
  list(df=idf, feats=ifeats)
}

# Load the average faces
load.ave <- function() {
  ## read in
  labs <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_ave_labels.csv", header=F)
  reps <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_ave_reps.csv", header=F)

  ## remove the average
  rm.ind <- grep("average_face2", labs$V2)
  labs <- labs[-rm.ind,]
  reps <- reps[-rm.ind,]
  
  ## format
  ifiles <- basename(as.character(labs$V2))
  pnames <- sub(".png", "", ifiles)
  idf <- data.frame(ind=1, oind=1:nrow(reps), 
                    target=pnames %in% targets, 
                    person=pnames, fn=ifiles)
  idf <- idf[order(idf$person),]
  idf$ind <- 1:nrow(idf)
  
  ## reorder
  ifeats <- reps[idf$oind,]
  
  list(df=idf, feats=ifeats)
}

load.ave2 <- function() {
  ## read in
  labs <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_ave_labels.csv", header=F)
  reps <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_ave_reps.csv", header=F)
  
  ## only keep the average
  keep.ind <- grep("average_face2", labs$V2)
  reps <- reps[keep.ind,]
  
  reps
}

# Load the fullset of faces
load.fullset <- function() {
  ## read in
  labs <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_fullset_labels.csv", header=F)
  reps <- read.csv("/data1/famface01/analysis/misc/openface/fampics/masked_fullset_reps.csv", header=F)
  
  ## format
  ifiles <- basename(as.character(labs$V2))
  inames <- sub(".png", "", ifiles)
  pnames <- sub("_large", "", sub("_medium", "", sub("_[0-9].*", "", inames)))
  idf <- data.frame(ind=1, oind=1:nrow(reps), 
                    target=pnames %in% targets, 
                    person=pnames, pic=inames, fn=ifiles)
  idf <- idf[order(idf$pic),]
  idf$ind <- 1:nrow(idf)
  
  ## reorder
  ifeats <- reps[idf$oind,]
  
  list(df=idf, feats=ifeats)
}

# Now actually load everything
face.norm     <- load.ave2()
faces.avg     <- load.ave()
faces.exp     <- load.exp()
faces.fullset <- load.fullset()

# save ordering of exp faces
exp.pics <- as.character(faces.exp$df$pic)
exp.pics2 <- as.character(faces.exp$df$pic[faces.exp$df$target])

# get only targets for exp
faces.exp.targ <- list()
faces.exp.targ$df <- faces.exp$df[faces.exp$df$target,]
faces.exp.targ$df$pic <- factor(faces.exp.targ$df$pic)
faces.exp.targ$df$person <- factor(faces.exp.targ$df$person)
faces.exp.targ$feats <- faces.exp$feats[faces.exp$df$target,]

faces.avg.targ <- list()
faces.avg.targ$df <- faces.avg$df[faces.avg$df$target,]
faces.avg.targ$df$person <- factor(faces.avg.targ$df$person)
faces.avg.targ$feats <- faces.avg$feats[faces.avg$df$target,]


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
ldata <- ddply(dfs, c("person", "pic"), summarise,
               N    = length(button_pressed), 
               mean = mean(button_pressed),
               sd   = sd(button_pressed),
               se   = sd / sqrt(N)
)
head(ldata)

# reorder to match...
oinds <- sapply(exp.pics2, function(epic) which(ldata$pic == epic))
all.equal(exp.pics2, ldata$pic[oinds])
ldata <- ldata[oinds,]
head(ldata)



# Load RT ------------------------------------------------------------------

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

library(xlsx)
subjects <- sprintf("sub%02i", 1:6)

# Load RT information
timing.df <- ldply(subjects, load.timing)
timing.df$pic <- sub(".jpg", "", basename(as.character(timing.df$fname)))
head(timing.df)
timing.df2 <- subset(timing.df, stimtype=='face' & target == 1)
rt.df <- ddply(timing.df2, .(name, pic), summarise, mean = mean(rt, na.rm=T))

# Select cdata and pictures cuz only show 18/20 in scanner
oinds <- sapply(rt.df$pic, function(x) which(ldata$pic == x))
all.equal(rt.df$pic, ldata$pic[oinds])
ldata2 <- ldata[oinds,]
ldata2$pic <- factor(ldata2$pic)

# Similarly select the features
oinds <- sapply(rt.df$pic, function(x) which(exp.pics2 == x))
faces.exp2 <- list()
faces.exp2$df <- faces.exp$df[faces.exp$df$target,][oinds,]
faces.exp2$df$pic <- factor(faces.exp2$df$pic)
faces.exp2$feats <- faces.exp$feats[faces.exp$df$target,][oinds,]
all.equal(as.character(faces.exp2$df$pic), rt.df$pic)



# Add gender+set+combine ---------------------------------------------------

## add gender label
female <- c("Angelina_Jolie", "Julia_Roberts", "Jennifer_Aniston", "Oprah_Winfrey")
male   <- c("Justin_Timberlake", "Will_Smith", "Brad_Pitt", "Johnny_Depp")
lst.gender <- factor((ldata2$person %in% male)*1, 
                     levels=c(0,1), labels=c("female", "male"))

## add set label
setA <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith")
setB <- c("Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
lst.set <- factor((ldata2$person %in% setB)*1, 
                  levels=c(0,1), labels=c("setA", "setB"))

## combine
behav.df2 <- data.frame(person=ldata2$person, pic=ldata2$pic, 
                        gender=lst.gender, set=lst.set, 
                        rt=rt.df$mean, likeness=ldata2$mean)


# Brief regressions to check --------------------------------------------

# Almost significant
X <- as.matrix(faces.exp2$feats)
all.equal(as.character(ldata2$pic), as.character(faces.exp2$df$pic))
summary(aov(lm(ldata2$mean ~ person + X, data=faces.exp2$df)))

# Now try the rt (oh...this is now significant)
X <- as.matrix(faces.exp2$feats)
all.equal(as.character(rt.df$pic), as.character(faces.exp2$df$pic))
summary(aov(lm(rt.df$mean ~ person + X, data=faces.exp2$df)))



# Save --------------------------------------------------------------------

setwd("/data1/famface01/command/misc/face_representations/322_asap_rois/210_fampics_setup")

## nnet
save(face.norm, faces.avg, faces.exp, faces.fullset, faces.exp.targ, faces.exp2, 
     faces.avg.targ, file="z_face_nnet_feats.rda")

## likeness and rt
write.csv(ldata, file="z_likeness_all.csv")
write.csv(behav.df2, file="z_likeness+rt+gender+set.csv")
