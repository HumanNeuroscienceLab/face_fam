#!/usr/bin/env bash


# Setup/Load --------------------------------------------------------------

# copy over the original peak coordinates and read them in
file.copy("/mnt/nfs/share/rois/face_gt_house+scene.txt", "/data1/famface01/analysis/roi/asap/")
orig.peaks <- read.table("/data1/famface01/analysis/roi/asap/face_gt_house+scene.txt")
colnames(orig.peaks) <- c("x", "y", "z", "num")



# Select Needed Ones ------------------------------------------------------

# rh: V1=21, OFA=3, pFFA=7, mFFA=1, aFFA=69, vATL=32, amygdala=5
# lh: V1=24, OFA=8, pFFA=(flip), mFFA=2, aFFA=62, vATL=26, amygdala=9

select.names <- c("V1", "OFA", "pFFA-1", "mFFA-1", "aFFA-1", "vATL", "amygdala")
select.nums.rh <- c(21, 3, 7, 1, 69, 32, 5)
select.nums.lh <- c(24, 8, 7, 2, 62, 26, 9)

rh.peaks <- orig.peaks[select.nums.rh,]
lh.peaks <- orig.peaks[select.nums.lh,]

rh.peaks$hemi <- "rh"
rh.peaks$name <- select.names
lh.peaks$hemi <- "lh"
lh.peaks$name <- select.names



# Add/Replace -----------------------------------------------------------

## replace right vatl=32 with 34 -6 -34 (this is now a mirror of the left vatl)
rh.peaks[rh.peaks$name=="vATL",1:3] <- c(34,-6,-34)

## pFFA=(flip) # -42 -66 -16
lh.peaks[lh.peaks$name=="pFFA-1",1] <- -1 * lh.peaks[lh.peaks$name=="pFFA-1",1]

## intermediate ones to create on both hemis...below are right hemi
# 42 -56 -12 # pFFA2
# 42 -38 -18 # mFFA2
# 38 -14 -28 # aFFA2
rh.peaks <- rbind(rh.peaks, c(42, -56, -12, 97, "rh", "pFFA-2"))
rh.peaks <- rbind(rh.peaks, c(42, -38, -18, 98, "rh", "mFFA-2"))
rh.peaks <- rbind(rh.peaks, c(38, -14, -28, 99, "rh", "aFFA-2"))
## lh
lh.peaks <- rbind(lh.peaks, c(-42, -56, -12, 97, "lh", "pFFA-2"))
lh.peaks <- rbind(lh.peaks, c(-42, -38, -18, 98, "lh", "mFFA-2"))
lh.peaks <- rbind(lh.peaks, c(-38, -14, -28, 99, "lh", "aFFA-2"))



# Reorder ----------------------------------------------------------------

# We want the numbers to be 100s for right and 200s for left
# Also we want to reorder them from posterior (to anterior)
new.order <- c("V1", "OFA", "pFFA-1", "pFFA-2", "mFFA-1", "mFFA-2", 
               "aFFA-1", "aFFA-2", "vATL", "amygdala")

rh.oinds <- sapply(new.order, function(x) which(rh.peaks$name==x))
rh.peaks <- rh.peaks[rh.oinds,]
rh.peaks$num <- 1:nrow(rh.peaks) + 100

lh.oinds <- sapply(new.order, function(x) which(lh.peaks$name==x))
lh.peaks <- lh.peaks[lh.oinds,]
lh.peaks$num <- 1:nrow(lh.peaks) + 200

# combine
all.peaks <- rbind(rh.peaks, lh.peaks)
rownames(all.peaks) <- 1:nrow(all.peaks)



# Save --------------------------------------------------------------------

odir <- "/data1/famface01/command/misc/face_representations/240_roi"
write.csv(all.peaks, file=file.path(odir, "z_asap_allpeaks.csv"))
write.table(all.peaks[,1:4], file=file.path(odir, "z_asap_peak_crds.txt"), 
            row.names=F, col.names=F, quote=F)

odir2 <- "/data1/famface01/analysis/roi/asap"
system(sprintf(". ~/.bashrc; 3dUndump -overwrite -srad 4 -master $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -xyz -orient LPI -prefix %s/asap_ventral_peaks.nii.gz %s/z_asap_peak_crds.txt", odir2, odir))


