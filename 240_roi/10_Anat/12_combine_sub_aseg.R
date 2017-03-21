
library(niftir)
library(plyr)

subjects <- sprintf("sub%02i", 1:6)

indirs <- "/data1/famface01/analysis/preprocessed/%s/anat/atlases"
indirs <- sprintf(indirs, subjects)
stdfile   <- sprintf("%s/data/standard/MNI152_T1_2mm_brain.nii.gz", Sys.getenv("FSLDIR"))
asegfiles <- sprintf("%s/aseg_to_std.nii.gz", indirs)

outdir <- "/data1/famface01/analysis/roi/10_Anat"
dir.create(outdir)
file.copy(stdfile, file.path(outdir, "standard_2mm.nii.gz"))

hdr <- read.nifti.header(stdfile)
asegs <- t(laply(asegfiles, read.mask, NULL, .progress="text"))

table(asegs[,1])

aseg.inds <- c(3, 8, 10:13, 16, 17:18, 26, 28, 42, 47, 49:54, 58, 60)
tab <- table(asegs[,2])
as.numeric(names(tab))[!(as.numeric(names(tab)) %in% aseg.inds)]

asegs.grey <- apply(asegs, 2, function(x) (x %in% aseg.inds))
asegs.grey <- (rowSums(asegs.grey) > 0) * 1
write.nifti(asegs.grey, hdr, outfile=file.path(outdir, "fs_grey_mask_allsubs.nii.gz"))

# Dilate
infile <- file.path(outdir, "fs_grey_mask_allsubs.nii.gz")
outfile <- file.path(outdir, "fs_grey_mask_allsubs_dil1-1.nii.gz")
cmd <- sprintf("3dmask_tool -input %s -dilate_inputs 1 -1 -prefix %s", infile, outfile)
cat(cmd, "\n")
system(cmd)






# THIS DOESNT WORK BECAUSE THE APARC+ASEG IS MESSEG UP
aparcfiles <- sprintf("%s/aparc+aseg_2_to_std.nii.gz", indirs)
aparcs <- t(laply(aparcfiles, read.mask, NULL, .progress="text"))
 
aparc0.inds <- c(1000:1035)
aparc.inds <- c(aseg.inds, aparc0.inds, aparc0.inds+1000, aparc0.inds+2000, aparc0.inds+3000)
aparcs.grey <- apply(aparcs, 2, function(x) (x %in% aparc.inds))
aparcs.grey <- (rowSums(aparcs.grey) > 0) * 1
write.nifti(aparcs.grey, hdr, outfile=file.path(outdir, "fs_grey2_mask_allsubs.nii.gz"))

