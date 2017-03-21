#' Apply the hyperalignment for each ROI to the functional unfam vids data
#' 
#' This time we will remove the covariates: compcor, motion, and FD (motion)
#' 

#--- Setup ---#

if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
library(methods)
library(plyr)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(20)
library(biganalytics)


#--- Variables/Paths ---#

subs <- sprintf("sub%02i", 1:6)
tbs  <- c( "tb3429", "tb3390", "tb3510", "tb3392", "tb3663", "tb3580" )
faceloc.subs <- as.character(read.table("/data1/faceloc02/notes/goodsubs.txt")[,1])
which.subs <- sapply(tbs, function(tb) which(faceloc.subs==tb))

# ROIs
roi.names <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", "r.Amyg", 
               "l.vATL", "l.FFA", "l.OFA", "l.EBA", "l.Amyg")
nrois <- length(roi.names)
roi.names1 <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", 
                "l.vATL", "l.FFA", "l.OFA", "l.EBA")
nrois1 <- length(roi.names1)
roi.names2 <- c("r.Amyg", "l.Amyg")
nrois2 <- length(roi.names2)

mrois1 <- read.mask("/data1/famface01/analysis/roi/Functional_v3/group_rois.nii.gz", NULL)
mask1 <- mrois1 != 0

mrois2 <- read.mask("/data1/famface01/analysis/roi/Functional_v3/amygdala/group_rois_select.nii.gz", NULL)
mask2 <- mrois2 != 0

mrois <- mrois1 * 0
for (i in 1:nrois) {
  if (any(roi.names[i] == roi.names1)) {
    ri <- which(roi.names[i] == roi.names1)
    mrois[mrois1 == ri] <- i
  } else if (any(roi.names[i] == roi.names2)) {
    ri <- which(roi.names[i] == roi.names2)
    mrois[mrois2 == ri] <- i
  } else {
    stop("error")
  }
}
rnames <- roi.names
grp.mask <- mrois!=0
rois <- mrois
roi.nums <- 1:length(roi.names)


pre.base  <- "/data1/famface01/analysis/preprocessed"


#--- Functions ---#

load.hyperalign <- function(rname) {
  load(sprintf("/data1/famface01/analysis/roi/Functional_v3/hyperalignments/hyper_%s.rda", rname))
  roi.mask <- mask[grp.mask]
  list(mask=roi.mask, transforms=hypers)
}

load.mc <- function(subj) {
  funcdir  <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  df.paths <- read.table(file.path(funcdir, "df_paths.txt"), header=T)
  inds     <- df.paths$inindex[df.paths$name == "unfam_vids"]
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

load.fds <- function(subj) {
  funcdir  <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  df.paths <- read.table(file.path(funcdir, "df_paths.txt"), header=T)
  inds     <- df.paths$inindex[df.paths$name == "unfam_vids"]
  
  fpaths   <- sprintf("%s/mc/func_run%02i_fd_rel.1D", funcdir, inds)
  fds_rel  <- ldply(fpaths, function(fpath) {
    x <- read.table(fpath)
    x <- as.matrix(x)
    x <- scale(x, scale=F, center=T)
    x
  })
  fds_rel <- as.matrix(fds_rel)
  
  fpaths   <- sprintf("%s/mc/func_run%02i_fd_abs.1D", funcdir, inds)
  fds_abs  <- ldply(fpaths, function(fpath) {
    x <- read.table(fpath)
    x <- as.matrix(x)
    x <- scale(x, scale=F, center=T)
    x
  })
  fds_abs <- as.matrix(fds_abs)
  
  cbind(rel=fds_rel, abs=fds_abs)
}

# only use top5
load.compcor <- function(subj) {
  funcdir <- file.path("/data1/famface01/analysis/preprocessed", subj, "func")
  fpath <- file.path(funcdir, "unfam_vids", "compcor_top10", "compcor_top10.1D")
  comps <- read.table(fpath)
  comps <- as.matrix(comps)[,1:5]
  
  comps
}

load.covars <- function(subj) {
  mcs <- load.mc(subj)
  fds <- load.fds(subj)
  #comps <- load.compcor(subj)
  cbind(mcs, fds)
}



#--- RUN ---#

cat("load all hyperalignments\n")
lst.hyps <- llply(roi.names, load.hyperalign, .progress="text")
names(lst.hyps) <- roi.names

lst.hyps.dat <- list()
for (si in 1:length(subs)) {
  sub <- subs[si]
  cat(sub, "\n")
  
  # Load the functional data
  cat("...load the functional data\n")
  func_files <- sort(Sys.glob(file.path(pre.base, sub, "func/unfam_vids/std_filtered_func*.nii.gz")))
  sdat <- ldply(func_files, function(func_file) {
    dat <- read.big.nifti(func_file)
    dat <- dat[,grp.mask]
    dat <- scale(dat, scale=F, center=T)
    dat
  }, .progress="text")
  sdat <- as.matrix(sdat)
  
  cat("...load covariates, remove but keep mean\n")
  #covars <- load.covars(sub)
  covars <- load.mc(sub)
  rdat <- lm(sdat ~ covars)$residuals
  mean.sdat <- colMeans(sdat)
  sdat2 <- sweep(rdat, 2, mean.sdat, FUN="+")
  rm(rdat)
  
  # Apply the alignment for each roi
  cat("...apply for each roi\n")
  hyp.dat <- llply(roi.names, function(rname) {
    hyper <- lst.hyps[[rname]]
    pret  <- hyper$transforms[[which.subs[si]]]
    x     <- sdat2[,hyper$mask]
    x.new <- pret$s * x %*% pret$R
    x.new <- sweep(x.new, 2, FUN="+", pret$tt)
    return(x.new)
  }, .parallel=T)
  names(hyp.dat) <- roi.names
  
  lst.hyps.dat[[sub]] <- hyp.dat
}

# Save
ofile <- "/data1/famface01/analysis/misc/320_roi_task_activity/"

# Load in the whole data and save it to it
cat("...save into rda\n")
indir  <- "/data1/famface01/analysis/encoding/SubjRois_Unfam/data"
infile <- sprintf("%s/roi_n_more_%s.rda", indir, sub)
ret    <- NULL
load(infile)
ret$basics$fmri$hyp.dat <- hyp.dat
outfile <- infile
save(ret, file=outfile)




