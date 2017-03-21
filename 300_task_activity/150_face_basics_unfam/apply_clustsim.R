#!/usr/bin/env Rscript

# sdir="/mnt/nfs/psych/faceMemoryMRI/analysis/subjects/tb9226/NoQuestions/task/preproc_spmg1.reml"
# ./apply_clustsim.R ${sdir}/ClustSim.NN3 ${sdir}/stats/zstat_bio.nii.gz ztmp_output.nii.gz

# Arguments

args      <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) cat("usage: apply_clustsim.R ClustSim-prefix vox-thr(pval) clust-thr(pval) input-file output-file\n")

cprefix   <- args[1]
vox_thr   <- as.numeric(as.character(args[2]))
clust_thr <- as.numeric(as.character(args[3]))
ifile     <- args[4]
ofile     <- args[5]


# Setup

file_1d   <- sprintf("%s.1D", cprefix)
file_niml <- sprintf("%s.niml", cprefix)

# get voxel thresh index
pthr      <- system(sprintf("grep pthr %s", file_niml), intern=T)
pthr      <- gsub("\"", "", sub("\ *pthr=", "", pthr))
pthr      <- as.numeric(as.character(strsplit(pthr, ",")[[1]]))
pind      <- which(pthr == vox_thr)
if (length(pind) == 0) stop("voxel thresh not found:", vox_thr)

# get cluster thresh index
athr      <- system(sprintf("grep athr %s", file_niml), intern=T)
athr      <- gsub("\"", "", sub("\ *athr=", "", athr))
athr      <- as.numeric(as.character(strsplit(athr, ",")[[1]]))
aind      <- which(athr == clust_thr)
if (length(aind) == 0) stop("voxel thresh not found:", clust_thr)

# get cluster size
# pthr (vox_thr) x athr (vox_thr)
x         <- as.matrix(read.table(file_1d))
clust_size<- as.integer(round(x[pind,aind]))
cat("cluster size:", clust_size, "\n")


# Run

# apply voxel size threshold
vox_zthr <- qt(vox_thr/2, Inf, lower.tail=F)
cmd <- sprintf("3dcalc -overwrite -a %s -expr 'a*astep(a,%f)' -prefix %s", ifile, vox_zthr, ofile)
cat(cmd, "\n")
system(cmd)

# apply cluster threshold
cmd <- sprintf("3dclust -overwrite -prefix %s -dxyz=1 1 %i %s", ofile, clust_size, ofile)
cat(cmd, "\n")
system(cmd)

#3dcalc -a /mnt/nfs/psych/faceMemoryMRI/analysis/subjects/tb9226/NoQuestions/task/preproc_spmg1.reml/stats/zstat_bio.nii.gz -expr 'step(a-1.96)' -prefix ztmp_output.nii.gz
#3dclust -overwrite -prefix ztmp_output2.nii.gz -1thresh 1.96 -dxyz=1 1 10 /mnt/nfs/psych/faceMemoryMRI/analysis/subjects/tb9226/NoQuestions/task/preproc_spmg1.reml/stats/zstat_bio.nii.gz
