#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 10a_main_subj_afni.sh {} ::: ${subjects[@]}

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
outdir="${outbase}/main_face_measures.reml"
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
  --stim-am1 scale ${tdir}/stimam_scale.txt 'SPMG1(2)' \
  --stim-am1 pose_x ${tdir}/stimam_pose_scores.txt 'SPMG1(2)' \
  --stim-am1 pose_y ${tdir}/stimam_pose_scores.1.txt 'SPMG1(2)' \
  --stim-am1 pose_r1 ${tdir}/stimam_pose_scores.2.txt 'SPMG1(2)' \
  --stim-am1 pose_r2 ${tdir}/stimam_pose_scores.3.txt 'SPMG1(2)' \
  --stim-am1 pose_mse ${tdir}/stimam_pose_mse.txt 'SPMG1(2)' \
  --stim-am1 mouth_s1 ${tdir}/stimam_mouth_scores.txt 'SPMG1(2)' \
  --stim-am1 mouth_s2 ${tdir}/stimam_mouth_scores.1.txt 'SPMG1(2)' \
  --stim-am1 asym ${tdir}/stimam_asym.txt 'SPMG1(2)' \
  --stim-am1 mean_face ${tdir}/stimam_mean_face.txt 'SPMG1(2)' \
  --stim-am1 pca_texture ${tdir}/stimam_pca_texture.txt 'SPMG1(2)' \
  --stim-am1 rel_scale ${tdir}/stimam_rel_scale.txt 'SPMG1(2)' \
  --stim-am1 rel_pose_x ${tdir}/stimam_rel_pose_scores.txt 'SPMG1(2)' \
  --stim-am1 rel_pose_y ${tdir}/stimam_rel_pose_scores.1.txt 'SPMG1(2)' \
  --stim-am1 rel_pose_r1 ${tdir}/stimam_rel_pose_scores.2.txt 'SPMG1(2)' \
  --stim-am1 rel_pose_r2 ${tdir}/stimam_rel_pose_scores.3.txt 'SPMG1(2)' \
  --stim-am1 rel_pose_mse ${tdir}/stimam_rel_pose_mse.txt 'SPMG1(2)' \
  --stim-am1 rel_mouth_s1 ${tdir}/stimam_rel_mouth_scores.txt 'SPMG1(2)' \
  --stim-am1 rel_mouth_s2 ${tdir}/stimam_rel_mouth_scores.1.txt 'SPMG1(2)' \
  --stim-am1 rel_mean_fds ${tdir}/stimam_rel_mean_fds.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
