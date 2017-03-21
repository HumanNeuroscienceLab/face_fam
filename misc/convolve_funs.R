spmt_ref <- function(times, peak_delay=6, under_delay=16, 
                     peak_disp=1, under_disp=1, 
                     p_u_ratio=6, normalize=TRUE) 
{
  hrf   <- vector("numeric", length(times))
  pos_t <- times[times > 0]
  peak  <- dgamma(pos_t, shape=peak_delay/peak_disp, scale=peak_disp)
  undershoot <- dgamma(pos_t, shape=under_delay/under_disp, scale=under_disp)
  hrf[times > 0] <- peak - undershoot / p_u_ratio
  
  if (normalize) hrf <- hrf / max(hrf)
  
  return(hrf)
}

spmt <- function(times) {
  spmt_ref(times)
}

dspmt <- function(times) {
  spmt(times) - spmt(times - 1)
}

ddspmt <- function(times) {
  (spmt(times) - spmt_ref(times, peak_disp=1.01))/0.01
}

# returns HRFs by voxels (rows) and time (columns)
spm_hrf <- function(hrfs, times) {
  hrfs[1,] %*% t(spmt(times)) + hrfs[2,] %*% t(dspmt(times)) + hrfs[3,] %*% t(ddspmt(times))
}

# This one uses the custom HRF (needs to know the subject id as well)
convolve.features.worker2 <- function(subj, vid.timing, features, tr=1, frate=1/24, 
                                      frames.sec=1/frate, to.convolve=TRUE) {
  # We upsample to the number of frames in the video
  # and do the convolution with that
  # and then downsample to the TRs we have
  
  lapply(roi.names, function(rname) {
    cat("====\n")
    cat(rname, "\n")
    
    # Get the canonical HRF to do the convolution
    basedir <- sprintf("/data1/famface01/analysis/task_activity/%s/unfam_vids.betas", subj)
    fname   <- file.path(basedir, sprintf("hrfs_%s.txt", rname))
    hrfs    <- as.matrix(read.table(fname))
    
    ntpts    <- 318
    nruns    <- max(vid.timing$run)
    tot.ntpts<- ntpts*nruns
    tpts     <- seq(frate, ntpts, frate) # 318 * 16runs = 2288 * 24frames
    chrf     <- spm_hrf(hrfs, tpts)
    
    # Upsample the features (faster way to do this?)
    cat("upsample\n")
    up.features <- big.matrix(tot.ntpts*frames.sec, ncol(features), shared=T, init=0)
    l_ply(1:nrow(vid.timing), function(ri) {
      #cat(ri, "\n")
      frame.row <- vid.timing[ri,]
      onset <- frame.row$onset * frames.sec + 1
      dur   <- frame.row$duration * frames.sec
      offset<- onset + dur - 1
      l_ply(1:ncol(features), function(fi) {
        up.features[onset:offset,fi] <- features[ri,fi]
      })
    }, .parallel=T)
    ## get runs (do convolution within runs)
    up.runs <- rep(1:nruns, each=ntpts * frames.sec)
    
    # Convolve
    cat("convolve\n")
    conv.up.features <- laply(1:ncol(features), function(fi) {
      lst.convs <- llply(1:nruns, function(irun) {
        vec <- up.features[up.runs==irun,fi]
        s.convs <- apply(chrf, 1, function(x) {
          s.conv <- convolve(x, rev(vec))
          s.conv <- s.conv/max(s.conv)
          s.conv[is.na(s.conv)] <- 0
          s.conv
        })
      })
      do.call(rbind, lst.convs)
    }, .parallel=T, .drop=F)
    
    # Downsample
    cat("downsample\n")
    ix <- seq(1, dim(conv.up.features)[2], tr/frate)
    conv.features <- conv.up.features[,ix,,drop=F]
    if (dim(conv.features)[2] != tot.ntpts) {
      stop("length of nrow(conv.features) not same as tot.tpts")
    }
    
    return(conv.features)
  })
}


convolve.fftw <- function(x, y) {
  library(fftw)
  n  <- length(x)
  ny <- length(y)
  if (n != ny) stop("x and y must be the same lengths")
  
  # pad
  n1 <- ny - 1
  x <- c(rep.int(0,n1), x)
  n <- length(y <- c(y, rep.int(0,n1)))
  
  # convolve
  p <- fftw::planFFT(x)
  x1 <- IFFT(FFT(x, plan=p) * Conj(FFT(y, plan=p)), plan=p, scale=F)
  x1 <- Re(x1)/n
  x1 <- x1[1:ny]
  
  return(x1)
}

convolve.mfftw <- function(x, ys) {
  library(fftw)
  n  <- length(x)
  ny <- nrow(ys)
  if (n != ny) stop("x and y must be the same lengths")
  
  # pad
  n1 <- ny - 1
  x <- c(rep.int(0,n1), x)
  n <- length(x)
  ys0 <- ys
  ys <- matrix(0, n, ncol(ys))
  ys[1:ny,] <- ys0
  rm(ys0)
  
  # convolve
  p <- fftw::planFFT(x)
  xx <- FFT(x, plan=p)
  x1s <- sapply(1:ncol(ys), function(i) {
    Re(IFFT(xx * Conj(FFT(ys[,i], plan=p)), plan=p, scale=F))
  })
  x1s <- x1s/n
  x1s <- x1s[1:ny,,drop=F]
  
  return(x1s)
}

convolve.features.worker <- function(vid.timing, features, tr=1, frate=1/24, 
                                     frames.sec=1/frate, to.convolve=TRUE, 
                                     parallel=T, verbose=T) {
  suppressMessages(library(neuRosim))
  # We upsample to the number of frames in the video
  # and do the convolution with that
  # and then downsample to the TRs we have
  
  # Get the canonical HRF to do the convolution
  ntpts    <- 318
  nruns    <- max(vid.timing$run)
  tot.ntpts<- ntpts*nruns
  tpts     <- seq(frate, ntpts, frate) # 318 * 16runs = 2288 * 24frames
  #chrf     <- canonicalHRF(tpts, verbose = FALSE)
  chrf     <- spmt(tpts)
  
  # Upsample the features (faster way to do this?)
  if (verbose) cat("upsample\n")
  ## setup upsampled data
  up.features <- big.matrix(tot.ntpts*frames.sec, ncol(features), 
                            shared=T, init=0)
  ## find new indices to put things
  onsets <- round(vid.timing$onset * frames.sec + 1)
  durs   <- round(vid.timing$duration * frames.sec)
  offsets<- sapply(1:length(onsets), function(i) max(onsets[i] + durs[i] - 1, onsets[i]))
  ns     <- offsets - onsets + 1
  up.inds <- unlist(lapply(1:length(onsets), function(i) onsets[i]:offsets[i]))
  ## corresponding indices in regularly sampled data
  reg.inds<- rep(1:length(onsets), ns)
  ## transfer data
  up.features[up.inds,] <- features[reg.inds,]
  ## get runs (do convolution within runs)
  up.runs <- rep(1:nruns, each=ntpts * frames.sec)
  
  # Convolve
  if (to.convolve) {
    if (verbose) cat("convolve\n")
    conv.up.features <- llply(1:nruns, function(irun) {
      ys <- up.features[up.runs==irun,,drop=F]
      ys <- ys[nrow(ys):1,,drop=F] # same as rev to each column
      s.convs <- convolve.mfftw(chrf, ys)
      for (i in 1:ncol(s.convs)) {
        s.convs[,i] <- s.convs[,i]/max(s.convs[,i]) # normalize
        s.convs[is.na(s.convs[,i]),i] <- 0
      }
      s.convs
    }, .parallel=parallel) # TODO: test if parallel is actually faster!
    conv.up.features <- do.call(rbind, conv.up.features)
    conv.up.features <- as.matrix(conv.up.features)
  } else {
    conv.up.features <- as.matrix(up.features)
  }
  
  # Downsample
  if (verbose) cat("downsample\n")
  ix <- seq(1, nrow(conv.up.features), tr/frate)
  conv.features <- conv.up.features[ix,,drop=F]
  if (nrow(conv.features) != tot.ntpts) {
    stop("length of nrow(conv.features) not same as tot.tpts")
  }
  
  return(conv.features)
}

convolve.features.worker0 <- function(vid.timing, features, tr=1, frate=1/24, 
                                      frames.sec=1/frate, to.convolve=TRUE) {
  suppressMessages(library(neuRosim))
  # We upsample to the number of frames in the video
  # and do the convolution with that
  # and then downsample to the TRs we have
  
  # Get the canonical HRF to do the convolution
  ntpts    <- 318
  nruns    <- max(vid.timing$run)
  tot.ntpts<- ntpts*nruns
  tpts     <- seq(frate, ntpts, frate) # 318 * 16runs = 2288 * 24frames
  chrf     <- canonicalHRF(tpts, verbose = FALSE)
  
  # Upsample the features (faster way to do this?)
  cat("upsample\n")
  up.features <- big.matrix(tot.ntpts*frames.sec, ncol(features), shared=T, init=0)
  l_ply(1:nrow(vid.timing), function(ri) {
    #cat(ri, "\n")
    frame.row <- vid.timing[ri,]
    onset <- round(frame.row$onset * frames.sec + 1)
    dur   <- round(frame.row$duration * frames.sec)
    offset<- max(onset + dur - 1, onset)
    n     <- offset - onset + 1
    for (fi in 1:ncol(features)) {
      up.features[onset:offset,fi] <- features[rep(ri,n),fi]
    }
  }, .parallel=T)
  ## get runs (do convolution within runs)
  up.runs <- rep(1:nruns, each=ntpts * frames.sec)
  
  # Convolve
  if (to.convolve) {
    cat("convolve\n")
    conv.up.features <- laply(1:ncol(features), function(fi) {
      lst.convs <- llply(1:nruns, function(irun) {
        vec <- up.features[up.runs==irun,fi]
        s.conv <- convolve(chrf, rev(vec), type="open")[1:length(vec)]
        s.conv <- s.conv/max(s.conv)
        s.conv[is.na(s.conv)] <- 0
        s.conv
      })
      unlist(lst.convs)
    }, .parallel=T)
    if (ncol(features) > 1) {
      conv.up.features <- t(conv.up.features)
    } else {
      conv.up.features <- as.matrix(conv.up.features)
    }
  } else {
    conv.up.features <- as.matrix(up.features)
  }
  
  # Downsample
  cat("downsample\n")
  ix <- seq(1, nrow(conv.up.features), tr/frate)
  conv.features <- conv.up.features[ix,,drop=F]
  if (nrow(conv.features) != tot.ntpts) {
    stop("length of nrow(conv.features) not same as tot.tpts")
  }
  
  return(conv.features)
}


convolve.features.byframe <- function(basics, features, inds=1:nrow(features), 
                                      hrf.sub=NULL, ...) 
{
  frame.timing <- basics$frame.timing[inds,,drop=F]
  features     <- features[inds,,drop=F]
  
  if (is.null(hrf.sub)) {
    conv.features <- convolve.features.worker(frame.timing, features, ...)
  } else {
    conv.features <- convolve.features.worker2(sub, frame.timing, features, ...)
  }
  
  conv.features
}

convolve.features.byvid <- function(basics, features, inds=1:nrow(features), 
                                    hrf.sub=NULL, ...) 
{
  func.timing  <- basics$timing[inds,,drop=F]
  features     <- features[inds,,drop=F]
  
  if (is.null(hrf.sub)) {
    conv.features <- convolve.features.worker(func.timing, features, ...)
  } else {
    conv.features <- convolve.features.worker2(sub, func.timing, features, ...)
  }
  
  conv.features
}

# only for demographics for now
convolve.features.demos <- function(basics, features, inds=1:nrow(features), 
                                    hrf.sub=NULL) 
{
  func.timing  <- basics$timing[inds,,drop=F]
  features     <- features[inds,,drop=F]
  
  cnames <- colnames(features)
  lst.rhs <- llply(cnames, function(cname) {
    cat("==", cname, "\n")
    f <- as.formula(sprintf("~ %s - 1", cname))
    rhs.frame <- model.frame(f, features, drop.unused.levels = TRUE)
    rhs       <- model.matrix(f, rhs.frame)
    rhs       <- rhs[inds,,drop=F]
    if (is.null(hrf.sub)) {
      conv.rhs  <- convolve.features.worker(func.timing, rhs)
    } else {
      conv.rhs  <- convolve.features.worker2(hrf.sub, func.timing, rhs)
    }
    
    conv.rhs
  })
  names(lst.rhs) <- cnames
  lst.rhs
}

