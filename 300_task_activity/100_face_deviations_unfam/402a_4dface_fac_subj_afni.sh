#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 402a_4dface_fac_subj_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=5
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/face_deviations_unfam"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/4dface_factor_analysis.reml"
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
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --stim-am1 static ${tdir}/stimam_4dface_fac_static.txt 'SPMG1(2)' \
  --stim-am1 texture ${tdir}/stimam_4dface_fac_texture.txt 'SPMG1(2)' \
  --stim-am1 dynamic ${tdir}/stimam_4dface_fac_dynamic.txt 'SPMG1(2)' \
  --stim-am1 expression ${tdir}/stimam_4dface_fac_expression.txt 'SPMG1(2)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
