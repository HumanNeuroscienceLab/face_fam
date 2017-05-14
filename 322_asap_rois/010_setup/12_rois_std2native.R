

idir <- "/data1/famface01/analysis/roi/asap"
odir <- file.path(idir, "subjects")
if (!file.exists(odir)) dir.create(odir)
std.roi <- sprintf("%s/asap_ventral_peaks.nii.gz", idir)

subjects <- sprintf("sub%02i", 1:6)

for (subj in subjects) {
  regdir  <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/reg", subj)
  mask    <- sprintf("/data1/famface01/analysis/preprocessed/%s/func/mask.nii.gz", subj)
  ofile   <- sprintf("%s/%s_asap_ventral_peaks.nii.gz", odir, subj)
  cmd <- sprintf("gen_applywarp.sh -f -i %s -m %s -r %s -w 'standard-to-exfunc' -o %s -t nn", std.roi, mask, regdir, ofile)
  cat(cmd, "\n")
  system(sprintf(". ~/.bash_profile; %s", cmd))
  cat("\n")
}

