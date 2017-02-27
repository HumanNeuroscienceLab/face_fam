# This script will generate the ROI time-series data in standard space


# Setup --------------------------------------------------------------------

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(niftir)
library(plyr)


# Load ROIs ---------------------------------------------------------------

# ROI Information
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
mask <- mrois!=0
rois <- mrois[mask]
roi.rnames <- factor(rois, levels=1:length(rnames), labels=rnames)


# Load Data ---------------------------------------------------------------

# Now that we have the rois, we need to load the ROI data
load.fmri <- function(subj, rois, mask) {
  preproc.base <- "/data1/famface01/analysis/preprocessed"
  
  # Masks
  rois      <- rois[mask]
  
  # Functional
  func.fnames <- Sys.glob(sprintf("%s/%s/func/unfam_vids/std_filtered_func_run*.nii.gz", preproc.base, subj))
  func.fnames <- sort(func.fnames)
  
  brain.dat <- ldply(func.fnames, function(funcfile) {
    func <- read.big.nifti(funcfile)
    func <- func[,mask]
    func <- scale(func, center=T, scale=F)
    func
  }, .progress="text")
  
  rois.dat <- llply(1:nrois, function(i) {
    brain.dat[,rois == i]
  }, .progress="text")
  rm(brain.dat)
  names(rois.dat) <- roi.names
  
  # also save which time-point is which run
  runs <- llply(1:length(func.fnames), function(i) {
    hdr <- read.nifti.header(func.fnames[1])
    rep(i, hdr$dim[4])
  }, .progress="text")
  runs <- unlist(runs)
  
  list(mask=mask, rois=rois, dat=rois.dat, runs=runs)
}
std.lst.rdats <- lapply(subjects, function(subj) {
  cat(subj, "\n")
  load.fmri(subj, mrois, mask)
})
names(std.lst.rdats) <- subjects


# Get ROI TS ---------------------------------------------------------------

## collapse across subjects
std.rdats.all <- llply(rnames, function(rname) {
  ldat <- lapply(std.lst.rdats, function(s.rdats) {
    s.rdats$dat[[rname]]
  })
  do.call(rbind, ldat)
}, .progress="text")
names(std.rdats.all) <- rnames


# Get Averages -------------------------------------------------------------

## get average voxelwise data
std.rdats <- sapply(std.rdats.all, rowMeans)

# Try with std.rdats based on thresholding of data
## cortex
locfiles1 <- sprintf("/data1/famface01/analysis/roi/Functional_v3/subjects/%s_zstats_thresh.nii.gz", subjects)
locs1 <- sapply(locfiles1, read.mask, NULL)
## amygdala
locfiles2 <- sprintf("/data1/famface01/analysis/roi/Functional_v3/amygdala/subjects/%s_zstats_thresh.nii.gz", subjects)
locs2 <- sapply(locfiles2, read.mask, NULL)
## combine
locs <- locs1
locs[locs2!=0] <- locs2[locs2!=0]
## average
std.rdats2 <- sapply(1:length(rnames), function(ri) {
  tmp <- lapply(1:ncol(locs), function(si) {
    subj <- subjects[si]
    loc.roi <- locs[rois == ri,si]
    rowMeans(std.rdats.all[[ri]][ssubs==subj,loc.roi!=0])
  })
  unlist(tmp)
})
colnames(std.rdats2) <- rnames
## threshold locs
locs.rois <- locs[rois!=0,]


# Save -------------------------------------------------------------

ssubs <- rep(subjects, each=nrow(std.rdats)/length(subjects))

## save for future
save(rnames, roi.rnames, mrois, mask, rois, std.rdats.all, std.rdats, 
     std.rdats2, locs.rois, ssubs, 
     file="/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda")

## load
#load("/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda", verbose=T)




