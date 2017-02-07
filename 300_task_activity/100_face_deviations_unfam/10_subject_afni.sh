#!/usr/bin/env bash

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 10_subject_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=4
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/face_deviations_unfam"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/whole_face_alt.reml"
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
  --stim faces ${tdir}/stim_faces.txt 'SPMG1(2)' \
  --stim-am1 whole_face ${tdir}/stimam_whole_dists.txt 'SPMG1(2)' \
  --stim-am1 mirror_transform ${tdir}/stimam_mirror_dists.txt 'SPMG1(2)' \
  --stim-am1 frame_diff ${tdir}/stimam_frame_dists.txt 'SPMG1(2)' \
  --stim-am1 scale ${tdir}/stimam_scale.txt 'SPMG1(2)' \
  --stim-am1 pose_x ${tdir}/stimam_pose1.txt 'SPMG1(2)' \
  --stim-am1 pose_y ${tdir}/stimam_pose2.txt 'SPMG1(2)' \
  --stim-am1 mouth_comp1 ${tdir}/stimam_mouth1.txt 'SPMG1(2)' \
  --stim-am1 mouth_comp2 ${tdir}/stimam_mouth2.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
