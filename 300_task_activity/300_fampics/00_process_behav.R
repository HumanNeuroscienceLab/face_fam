#' There are a few things here:
#' - face vs nonface
#' - 4 ids vs distractor ids + accuracy in identifying
#' - RT



# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(xlsx)
library(plyr)
suppressMessages(library(doMC))
registerDoMC(20)

# output
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/300_fampics/timings"


# Load --------------------------------------------------------------------

targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
             "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")

subj <- "sub02"
subjects <- sprintf("sub%02i", 1:6)


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


for (subj in subjects) {
  cat(subj, "\n")
  
  infiles <- list.files(sprintf("/data1/famface01/data/behavior/%s/fam_pics", subj), 
                        full.names = T)
  timing.df <- ldply(infiles, function(infile) {
    timing <- read.xlsx(infile, 1)
    timing <- timing[!is.na(timing$onset2_raw),]
    timing$sess <- as.numeric(as.character(timing$sess))
    timing
  }, .progress="text")
  
  # Add the actual (total) run number
  timing.df$tot.run <- (timing.df$sess-1)*2 + timing.df$run
  
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
  
  # Reformat everything
  timing.df3 <- subset(timing.df2, select=c("tot.run", "name", "onset", "onset2_raw", "trial", "num", "fname", "stimtype", "dur", "RT_raw", "target", "noresp", "incorrect"))
  
  cat("...incorrect\n")
  print(table(timing.df3$incorrect))
  cat("...no responses\n")
  print(table(timing.df3$noresp))
  
  
  
  # Timings --------------------------------------------------------------------
  
  # REMEMBER THESE TIMES ARE LOCAL
  
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  
  # Face 
  inds <- timing.df3$stimtype=="face"
  lines <- afni.timing(timing.df3[inds,])
  save.tofile(lines, "stim_faces.txt")
  
  # Incorrect
  inds <- timing.df3$incorrect==1
  lines <- afni.timing(timing.df3[inds,])
  save.tofile(lines, "stim_incorrect.txt")
  
  # No Response
  inds <- timing.df3$noresp==1
  lines <- afni.timing(timing.df3[inds,])
  save.tofile(lines, "stim_noresp.txt")
  
  # RT (amp)
  rts  <- as.numeric(as.character(timing.df3$RT_raw))
  inds <- !is.na(rts)
  amps <- scale(rts[inds], center=T, scale=F)
  lines <- afni.timing.amp(timing.df3[inds,], amps)
  save.tofile(lines, "stimam_rt.txt")
  
  # RT (dur)
  rts  <- as.numeric(as.character(timing.df3$RT_raw))
  inds <- !is.na(rts)
  durs <- rts[inds]
  lines <- afni.timing.dur(timing.df3[inds,], durs)
  save.tofile(lines, "stimdur_rt.txt")
  
  # Motion
  mcs <- load.mc(subj)
  ofile <- file.path(sub.outdir, "motion.1D")
  write.table(mcs, file=ofile, row.names=F, col.names=F, quote=F)
}
