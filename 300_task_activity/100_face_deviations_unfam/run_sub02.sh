#!/usr/bin/env bash

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=8
subj="sub02"

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/task_activity/unfam_face_errors/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/face_deviations_unfam"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/whole_face.reml"
rm -r ${outdir} 2> /dev/null

echo ${outdir}

${tog}/task_analysis.rb -i ${fundir}/unfam_vids/filtered_func_run??.nii.gz \
  -m ${fundir}/mask.nii.gz \
  -b ${fundir}/mean_func.nii.gz \
  --output ${outdir} \
  --global \
  --tr 1 \
  --polort 0 \
  --motion ${tdir}/motion.1D \
  --stim-am2 all ${tdir}/stimam2_all.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}

#  --stim-am2 whole_face ${tdir}/stimam_whole_dists.txt 'SPMG1(2)' \
#  --stim-am2 scale_1 ${tdir}/stimam_scale1.txt 'SPMG1(2)' \
#  --stim-am2 scale_2 ${tdir}/stimam_scale2.txt 'SPMG1(2)' \
#  --stim-am2 pose_x ${tdir}/stimam_pose1.txt 'SPMG1(2)' \
#  --stim-am2 pose_y ${tdir}/stimam_pose2.txt 'SPMG1(2)' \
