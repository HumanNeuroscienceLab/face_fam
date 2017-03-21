#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 52a_pred_subj_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=4
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/150_face_basics_unfam/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/face_basics_unfam"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/pred_demo_trait_feats.reml"
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
  --stim-am1 pred_demo_age ${tdir}/stimam_age_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_demo_gender ${tdir}/stimam_gender_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_demo_glasses ${tdir}/stimam_glasses_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_demo_makeup ${tdir}/stimam_makeup_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_t01_unemotional ${tdir}/stimam_trait01_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_t02_competent ${tdir}/stimam_trait02_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_t03_trustworthy ${tdir}/stimam_trait03_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_t04_memorable ${tdir}/stimam_trait04_pred.txt 'SPMG1(2)' \
  --stim-am1 pred_t05_attractive ${tdir}/stimam_trait05_pred.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
