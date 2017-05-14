#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 342a_avg1_subj.sh {} ::: ${subjects[@]}

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
outdir="${outbase}/nnet2_pca_avg+knn+covars.reml"
rm -r ${outdir} 2> /dev/null

echo ${outdir}

${tog}/task_analysis.rb -i ${fundir}/unfam_vids/filtered_func_run??.nii.gz \
  -m ${fundir}/mask.nii.gz \
  -b ${fundir}/mean_func.nii.gz \
  --output ${outdir} \
  --oresiduals \
  --global \
  --tr 1 \
  --polort 0 \
  --motion ${tdir}/motion.1D \
  --stim faces ${tdir}/stim_faces.txt 'SPMG1(2)' \
  --stim-am1 avg_dist ${tdir}/stimam_nnet2_pca_mean_diff.txt 'SPMG1(2)' \
  --stim-am1 knn_dist ${tdir}/stimam_nnet2_pca_knn_diff.txt 'SPMG1(2)' \
  --stim-am1 pose ${tdir}/stimam_nnet2_pca_pose.txt 'SPMG1(2)' \
  --stim-am1 frame_diff ${tdir}/stimam_nnet2_pca_frame_mean_diff.txt 'SPMG1(2)' \
  --stim-am1 frame_pose ${tdir}/stimam_nnet2_pca_frame_pose.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
