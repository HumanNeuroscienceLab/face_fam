# Setup -------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(xlsx)
library(plyr)
suppressMessages(library(doMC))
registerDoMC(20)

# output
outdir <- "/data1/famface01/command/misc/face_representations/300_task_activity/400_imagine/timings"



# Load --------------------------------------------------------------------

targets <- c("Angelina_Jolie", "Julia_Roberts", "Justin_Timberlake", "Will_Smith", 
             "Jennifer_Aniston", "Oprah_Winfrey", "Brad_Pitt", "Johnny_Depp")
#subj <- "sub02"
subjects <- sprintf("sub%02i", 1:6)

afni.timing <- function(tdf, nruns=4) {
  lines <- sapply(1:nruns, function(i) {
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

afni.timing.amp <- function(tdf, amps, nruns=4) {
  if (nrow(tdf) != length(amps)) stop("tdf and amps lengths differ")
  
  lines <- sapply(1:nruns, function(i) {
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

afni.timing.dur <- function(tdf, durs, nruns=4) {
  if (nrow(tdf) != length(durs)) stop("tdf and amps lengths differ")
  
  lines <- sapply(1:nruns, function(i) {
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

load.mc <- function(subj, task) {
  funcdir  <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  df.paths <- read.table(file.path(funcdir, "df_paths.txt"), header=T)
  inds     <- df.paths$inindex[df.paths$name == task]
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
  
  # For sub01 we are missing the timing files but the data is still there
  # So for runs1-2, we will use subj02's data
  
  infiles <- list.files(sprintf("/data1/famface01/data/behavior/%s/imagine", subj), 
                        full.names = T)
  if (subj == "sub01") {
    infiles0 <- list.files(sprintf("/data1/famface01/data/behavior/%s/imagine", "sub02"), 
                           full.names = T)
    infiles <- c(infiles0[1:2], infiles)
  }
  timing.df <- ldply(1:length(infiles), function(i) {
    infile <- infiles[i]
    timing <- read.xlsx(infile, 1)
    timing <- timing[!is.na(timing$onset2_raw),]
    timing$sess <- i
    timing$run <- 1
    timing$tot.run <- i
    timing
  }, .progress="text")
  timing.df$onset <- as.numeric(as.character(timing.df$onset))
  timing.df$celeb <- factor(as.character(timing.df$celeb))
  print(table(timing.df$celeb))
  
  
  # Timings --------------------------------------------------------------------
  
  # REMEMBER THESE TIMES ARE LOCAL
  
  sub.outdir <- file.path(outdir, subj)
  if (!file.exists(sub.outdir)) dir.create(sub.outdir)
  
  save.tofile <- function(lines, ofile) {
    writeLines(lines, file.path(sub.outdir, ofile))
  }
  
  # Imagine each celeb
  #celebs <- names(table(timing.df$celeb))
  celebs <- targets
  for (celeb in celebs) {
    inds <- timing.df$celeb == celeb
    lines <- afni.timing(timing.df[inds,,drop=F], nruns=max(timing.df$tot.run))
    save.tofile(lines, sprintf("stim_celeb_%s.txt", celeb))
  }
  
  # Motion
  mcs   <- load.mc(subj, "imagine")
  ofile <- file.path(sub.outdir, "motion.1D")
  write.table(mcs, file=ofile, row.names=F, col.names=F, quote=F)
}

# Manually add 
# Note: subj 5 has fewer
