# Extracts the time-series from each subjects ROIs

if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

suppressMessages(library(niftir))
library(plyr)
idir <- "/data1/famface01/analysis/roi/asap/subjects"

subjects <- sprintf("sub%02i", 1:6)

for (subj in subjects) {
  cat(subj, "\n")
  roifile <- sprintf("%s/%s_asap_ventral_peaks.nii.gz", idir, subj)
  
  funcdir  <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/fam_pics", subj)
  funcfiles <- sort(list.files(funcdir, pattern="^filtered_func_run[0-9]{2}[.]nii[.]gz"))
  
  tsdir  <- file.path(funcdir, "rois_asap_ventral_peaks")
  dir.create(tsdir, showWarnings=F)
  
  for (funcfile in funcfiles) {
    tsfile <- sub("filtered_func", "asap_ventral_peaks", funcfile)
    tsfile <- sub("[.]nii[.]gz", ".1D", tsfile)
    
    cmd <- sprintf("3dROIstats -mask %s -nobriklab -1Dformat %s > %s", roifile, file.path(funcdir, funcfile), 
                   file.path(tsdir, tsfile))
    cat(cmd, "\n")
    system(sprintf(". ~/.bash_profile; %s", cmd))
  }
  
  cat("\n")
}

# loop through to combine the different run ROIs
# make sure to center each run as well
for (subj in subjects) {
  cat(subj, "\n")
  
  funcdir  <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/fam_pics", subj)
  tsdir    <- file.path(funcdir, "rois_asap_ventral_peaks")
  tsfiles  <- sort(list.files(tsdir, pattern="^asap_ventral_peaks_run[0-9]{2}[.]1D"))
  
  ts.mat <- ldply(tsfiles, function(tsfile) {
    ts <- read.table(file.path(tsdir, tsfile))
    scale(ts, center=T, scale=F)
  })
  ts.mat <- as.matrix(ts.mat)
  
  # save
  ofile <- file.path(tsdir, "asap_ventral_peaks_all_centered.1D")
  write.table(ts.mat, file=ofile, row.names=F, col.names=F, quote=F)
  
  cat("\n")
}

## not centered
for (subj in subjects) {
  cat(subj, "\n")
  
  funcdir  <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/fam_pics", subj)
  tsdir    <- file.path(funcdir, "rois_asap_ventral_peaks")
  tsfiles  <- sort(list.files(tsdir, pattern="^asap_ventral_peaks_run[0-9]{2}[.]1D"))
  
  ts.mat <- ldply(tsfiles, function(tsfile) {
    ts <- read.table(file.path(tsdir, tsfile))
    ts
  })
  ts.mat <- as.matrix(ts.mat)
  
  # save
  ofile <- file.path(tsdir, "asap_ventral_peaks_all.1D")
  write.table(ts.mat, file=ofile, row.names=F, col.names=F, quote=F)
  
  # save runs...
  runs <- ldply(tsfiles, function(tsfile) {
    ts <- read.table(file.path(tsdir, tsfile))
    run <- as.numeric(sub(".1D", "", sub("asap_ventral_peaks_run", "", tsfile)))
    as.matrix(rep(run, nrow(ts)))
  })
  ofile <- file.path(tsdir, "runs.1D")
  write.table(runs, file=ofile, row.names=F, col.names=F, quote=F)
  
  cat("\n")
}

## easier to use with AFNI
for (subj in subjects) {
  cat(subj, "\n")
  
  funcdir  <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/fam_pics", subj)
  tsdir    <- file.path(funcdir, "rois_asap_ventral_peaks")
  tsfiles  <- sort(list.files(tsdir, pattern="^asap_ventral_peaks_run[0-9]{2}[.]1D"))
  
  # uncentered
  ts.mat <- ldply(tsfiles, function(tsfile) {
    ts <- read.table(file.path(tsdir, tsfile))
    
    ofile <- file.path(tsdir, sub(".1D", ".nii.gz", tsfile))
    tmp <- t(ts)
    dim(tmp) <- c(dim(tmp)[1], 1, 1, dim(tmp)[2])
    tmp <- as.nifti(tmp)
    write.nifti(tmp@.Data, tmp@header, outfile=ofile, overwrite=T)
    
    ts
  })
  ts.mat <- as.matrix(ts.mat)
  ofile <- file.path(tsdir, "asap_ventral_all.nii.gz")
  tmp <- as.nifti(ts.mat)
  write.nifti(tmp@.Data, tmp@header, outfile=ofile, overwrite=T)
  
  # centered
  ts.mat <- ldply(tsfiles, function(tsfile) {
    ts <- read.table(file.path(tsdir, tsfile))
    scale(ts, center=T, scale=F)
  })
  ts.mat <- as.matrix(ts.mat)
  ofile <- file.path(tsdir, "asap_ventral_all_centered.nii.gz")
  tmp <- as.nifti(ts.mat)
  write.nifti(tmp@.Data, tmp@header, outfile=ofile, overwrite=T)
  
  cat("\n")
}
