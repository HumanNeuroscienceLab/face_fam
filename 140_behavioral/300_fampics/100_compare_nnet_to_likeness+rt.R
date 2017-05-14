




# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(plyr)
suppressMessages(library(doMC))
registerDoMC(20)



# Load --------------------------------------------------------------------

fpaths <- list.files("/data1/famface01/data/behavior/ratings", full.names = T)

fpath <- fpaths[1]

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

cholMaha <- function(X) {
  dec <- chol( cov(X) )
  tmp <- forwardsolve(t(dec), t(X) )
  Rfast::Dist(t(tmp))
}


# Load nnet
nn.ave <- load.nn.ave()
nn.exp <- load.nn.exp()
nn.full <- load.nn.full()

# Load RT information
subjects <- sprintf("sub%02i", 1:6)
timing.df <- ldply(subjects, load.timing)
timing.df$pic <- sub(".jpg", "", basename(as.character(timing.df$fname)))
head(timing.df)
timing.df2 <- subset(timing.df, stimtype=='face')
rt.df <- ddply(timing.df2, .(name, pic), summarise, mean = mean(rt, na.rm=T))
#rt.df$pic <- factor(rt.df$pic)
cinds <- sapply(cdata$pic, function(x) x %in% rt.df$pic)
cdata2 <- cdata[cinds,]
cdata2$pic <- factor(cdata2$pic)
cdata2$person <- factor(cdata2$person)
rinds <- sapply(cdata2$pic, function(x) which(x==rt.df$pic))
all.equal(as.character(cdata2$pic), as.character(rt.df$pic[rinds]))
rt.df2 <- rt.df[rinds,]

# So this isn't significant
summary(lm(rt.df2$mean ~ cdata2$mean))
summary(lm(cdata2$mean ~ rt.df2$mean))
summary(lm(cdata2$mean ~ log(rt.df2$mean)))

ninds <- sapply(as.character(cdata2$pic), function(x) which(x==nn.exp$labs$img))
feats <- nn.exp$feats[ninds,]
labs  <- nn.exp$labs[ninds,]
all.equal(as.character(cdata2$pic), as.character(labs$img))
ds <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- Rfast::dista(ave.feat, cur.feat)
  d
})
ds2 <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  cur.feat <- feats[i,]
  comp.feat <- nn.full$feats[nn.full$labs$person==labs$person[i],]
  d <- Rfast::dista(cur.feat, t(comp.feat))
  mean(d)
})
sets <- rep(c("SetA", "SetB"), each=4)
targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
             "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
ds.compare <- t(sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  iis <- nn.ave$labs$img %in% targets
  ave.feat <- nn.ave$feats[iis,]
  cur.feat <- feats[i,]
  d <- Rfast::dista(cur.feat, t(ave.feat))
  names(d) <- nn.ave$labs$img[iis]
  #d[sets != sets[labs$person[i] == targets]] <- 0
  d
}))

summary(aov(cdata2$mean ~ ds))
summary(aov(cdata2$mean ~ lm(ds ~ ds.compare)$residuals))
summary(aov(cdata2$mean ~ cdata2$person + ds))
summary(aov(cdata2$mean ~ cdata2$person + ds + ds.compare))
summary(lm(cdata2$mean ~ cdata2$person + ds + ds.compare))
summary(aov(rt.df2$mean ~ cdata2$person + ds + ds.compare))
summary(lm(rt.df2$mean ~ cdata2$person + ds + ds.compare))


summary(aov(cdata2$mean ~ cdata2$person + ds + ds.compare))
summary(aov(rt.df2$mean ~ cdata2$person + ds + ds.compare))

cs <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- cor(ave.feat, cur.feat)
  d
})
ms1 <- unlist(lapply(unique(labs$person), function(name) {
  x <- rbind(feats[labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  #x <- rbind(nn.exp$feats[nn.exp$labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  Sx <- cov(x)
  d <- Rfast::mahala(x, nn.ave$feats[nn.ave$labs$img==name,], Sx)
  d[1:sum(labs$person==name)]
}))
ms2 <- unlist(lapply(unique(labs$person), function(name) {
  x  <- rbind(feats[labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  x2 <- rbind(nn.ave$feats[nn.ave$labs$img==name,], x)
  pd <- cholMaha(x2)
  pd[1,-1][1:sum(labs$person==name)]
}))

# Try
summary(lm(cdata2$mean ~ ds + cs + ms1 + ms2))
summary(lm(cdata2$mean ~ rt.df2$mean + ds + cs + ms1 + ms2))
summary(lm(cdata2$mean ~ ds + cs)) # not this though!
summary(lm(rt.df2$mean ~ ds + cs)) # significantly explain together!
summary(lm(rt.df2$mean ~ ds))
summary(lm(rt.df2$mean ~ cs))


library(MASS)

ninds <- sapply(as.character(cdata2$pic), function(x) which(x==nn.exp$labs$img))
feats <- nn.exp$feats[ninds,]
labs  <- nn.exp$labs[ninds,]
all.equal(as.character(cdata2$pic), as.character(labs$img))

isopts <- ldply(unique(labs$person), function(name) {
  cat(name, "\n")
  mat <- rbind(nn.ave$feats[nn.ave$labs$img==name,], 
               feats[labs$person==name,], 
               nn.full$feats[nn.full$labs$person==name,])
  d <- Rfast::Dist(mat)
  if (name == "Oprah_Winfrey") {
    d <- d[-19,][,-19]
  }
  isomat <- isoMDS(as.dist(d))
  pts <- isomat$points
  if (name == "Oprah_Winfrey") {
    pts <- rbind(pts[1:18,], pts[4,], pts[19:nrow(pts),])
  }
  #head(pts)
  cpts0 <- sweep(pts, 2, pts[1,])
  cpts  <- cpts0[-1,][1:sum(labs$person==name),]
  cpts
})
dim(isopts)
colnames(isopts) <- c("x", "y")

summary(aov(cdata2$mean ~ cdata2$person + ds))
summary(aov(cdata2$mean ~ cdata2$person + ds + isopts$x*isopts$y))








vi.df <- as.data.frame(vi$importance)
ds.classify <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  sel.dims <- vi.df[[labs$person[i]]] != 0
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- Rfast::dista(ave.feat[sel.dims], cur.feat[sel.dims])
  d
})
cs.classify <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  sel.dims <- vi.df[[labs$person[i]]] != 0
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- cor(ave.feat[sel.dims], cur.feat[sel.dims])
  d
})
ms1.classify <- unlist(lapply(unique(labs$person), function(name) {
  sel.dims <- vi.df[[name]] != 0
  x <- rbind(feats[labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  #x <- rbind(nn.exp$feats[nn.exp$labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  Sx <- cov(x[,sel.dims])
  d <- Rfast::mahala(x[,sel.dims], nn.ave$feats[nn.ave$labs$img==name,sel.dims], Sx)
  d[1:sum(labs$person==name)]
}))
ms2.classify <- unlist(lapply(unique(labs$person), function(name) {
  sel.dims <- vi.df[[name]] != 0
  x  <- rbind(feats[labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  x2 <- rbind(nn.ave$feats[nn.ave$labs$img==name,], x)
  pd <- cholMaha(x2[,sel.dims])
  pd[1,-1][1:sum(labs$person==name)]
}))

cor(cbind(ds, cs, dsc=ds.classify, csc=cs.classify))

summary(lm(cdata2$mean ~ ds + cs))
summary(lm(rt.df2$mean ~ ds + cs))

summary(lm(cdata2$mean ~ ds.classify + cs.classify))
summary(lm(rt.df2$mean ~ ds.classify + cs.classify))

summary(lm(cdata2$mean ~ ds + cs + ds.classify + cs.classify)) # not this though!
summary(lm(rt.df2$mean ~ ds + cs + ds.classify + cs.classify)) # significantly explain together!
summary(lm(cdata2$mean ~ ds + ds.classify))
summary(lm(rt.df2$mean ~ ds + ds.classify)) # significantly 

summary(lm(cdata2$mean ~ ms1 + ms1.classify))
summary(lm(rt.df2$mean ~ ms1 + ms1.classify))

fit <- lm(cdata2$mean ~ cdata2$person + ds)
summary(fit) # sig
plot(cdata2$mean, ds)
x1 <- lm(cdata2$mean ~ cdata2$person)$resid + mean(cdata2$mean)
y1 <- lm(ds ~ cdata2$person)$resid + mean(ds)
plot(x1, y1, 
     ylab="Deviance from Ave Famous Template", xlab="Likeness Rating")
fit <- lm(y1 ~ x1)
abline(a=fit$coefficients[1], b=fit$coefficients[2], col=2)

summary(aov(cdata2$mean ~ cdata2$person + ds)) # sig
summary(aov(cdata2$mean ~ cdata2$person + ds2)) # sig
summary(aov(cdata2$mean ~ cdata2$person + cs)) # sig
summary(aov(cdata2$mean ~ cdata2$person + ms1))
summary(aov(cdata2$mean ~ cdata2$person + ms2)) # wow really sig
summary(aov(cdata2$mean ~ cdata2$person + ds.classify))
summary(aov(cdata2$mean ~ cdata2$person + cs.classify))
summary(aov(cdata2$mean ~ cdata2$person + ms1.classify))
summary(aov(cdata2$mean ~ cdata2$person + ms2.classify))

summary(aov(rt.df2$mean ~ cdata2$person + ds))
summary(aov(rt.df2$mean ~ cdata2$person + ds2))
summary(aov(rt.df2$mean ~ cdata2$person + cs))
summary(aov(rt.df2$mean ~ cdata2$person + ms1))
summary(aov(rt.df2$mean ~ cdata2$person + ms2))
summary(aov(rt.df2$mean ~ cdata2$person + ds.classify))
summary(aov(rt.df2$mean ~ cdata2$person + cs.classify))
summary(aov(rt.df2$mean ~ cdata2$person + ms1.classify))
summary(aov(rt.df2$mean ~ cdata2$person + ms2.classify))


summary(aov(rt.df2$mean ~ cdata2$person + ds.classify))
summary(aov(rt.df2$mean ~ cdata2$person + cs.classify))
summary(aov(rt.df2$mean ~ cdata2$person + ms1.classify))


summary(lm(cdata2$mean ~ cdata2$person + ds.classify + cs.classify + ms1.classify))
summary(lm(rt.df2$mean ~ ds.classify + cs.classify + ms1.classify))
summary(lm(rt.df2$mean ~ ds.classify))
summary(lm(cdata2$mean ~ ds.classify))


head(cdata)
inds <- sapply(cdata$pic, function(x) which(x==nn.exp$labs$img))
all.equal(as.character(cdata$pic), as.character(nn.exp$labs$img[inds]))
nn.exp$labs
feats <- nn.exp$feats[inds,]
labs  <- nn.exp$labs[inds,]
ds <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- Rfast::dista(ave.feat, cur.feat)
  d
})
cs <- sapply(1:nrow(feats), function(i) {
  # distance of average to each of these pics
  ave.feat <- nn.ave$feats[labs$person[i] == nn.ave$labs$img,]
  cur.feat <- feats[i,]
  d <- cor(ave.feat, cur.feat)
  d
})

ms1 <- lapply(unique(labs$person), function(name) {
  x <- rbind(feats[labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  Sx <- cov(x)
  d <- Rfast::mahala(x, nn.ave$feats[i,], Sx)
  d[1:sum(labs$person==name)]
})
ms1 <- unlist(ms1)



ms2 <- lapply(unique(labs$person), function(name) {
  x  <- rbind(feats[labs$person==name,], nn.full$feats[nn.full$labs$person==name,])
  x2 <- rbind(nn.ave$feats[nn.ave$labs$img==name,], x)
  pd <- cholMaha(x2)
  pd[1,-1][1:sum(labs$person==name)]
})
ms2 <- unlist(ms2)

summary(lm(cdata$mean ~ ds + cs))
summary(lm(cdata$mean ~ ds + cs + ms1 + ms2))
summary(lm(cdata$mean ~ ms1))
summary(lm(cdata$mean ~ ms2))
summary(lm(cdata$mean ~ cs + ms1 + ms2))
summary(lm(cdata$mean ~ ms1 + ms2))



# TODO: create regressors relative to each of the average templates






# Classify ----------------------------------------------------------------

# Packages
library(caret)
library(glmnet)
library(doMC) # for parallel processing
registerDoMC(24) # set number of threads
library(plyr)

# To setup options for training model
fitControl <- trainControl(
  method = "repeatedcv", 
  number = 2, 
  repeats = 10, 
  summaryFunction = multiClassSummary, 
  classProbs = TRUE, 
  savePredictions = TRUE, 
  allowParallel = T
)

finds <- nn.full$labs$person %in% levels(cdata2$person)
X <- nn.full$feats[finds,]
y <- factor(nn.full$labs$person[finds])

#alphas <- c(0, 0.5, 1)
alphas <- c(1)
nlambda <- 100
tuneGrid <- ldply(alphas, function(alpha) {
  tmp <- glmnet(X, y, family="multinomial", nlambda=nlambda, alpha=alpha)
  data.frame(alpha=alpha, lambda=tmp$lambda)
}, .parallel=T)

# Run glmnet
fit <- train(X, y, # X = observations x featurs, y = observations
             method = "glmnet",
             trControl = fitControl, 
             tuneGrid = tuneGrid, 
             family="multinomial", 
             preProcess = c("center", "scale"), 
             metric="Kappa")
print(fit) # gives you summary
fit$bestTune
getTrainPerf(fit) # really really accurate
confusionMatrix(fit)
vi <- varImp(fit, scale=F)
apply(vi$importance, 2, function(x) mean(x==0))

summary_caret <- function(fit) {
  perf  <- getTrainPerf(fit)
  tune  <- fit$bestTune
  vi    <- varImp(fit, scale=F)
  pfeats<- mean(vi$importance!=0)*100
  
  df <- data.frame(roi=name, 
                   accuracy=perf$TrainAccuracy, kappa=perf$TrainKappa, 
                   alpha=tune$alpha, lambda=tune$lambda, 
                   percent.features=pfeats)
  
  df
}
