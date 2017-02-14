
#!/usr/bin/env bash

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 32_demos_gs_subj_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=4
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/150_face_basics_unfam/timings/${subj}"
tdir2="${base}/command/misc/face_representations/300_task_activity/100_face_deviations_unfam/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/face_basics_unfam"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/agegender_givenshape.reml"
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
  --stim-am1 whole_face ${tdir2}/stimam_whole_dists.txt 'SPMG1(2)' \
  --stim-am1 shape_gender ${tdir}/stimam_shape_gender.txt 'SPMG1(2)' \
  --stim-am1 resid_gender ${tdir}/stimam_resid_gender.txt 'SPMG1(2)' \
  --stim-am1 shape_age ${tdir}/stimam_shape_age.txt 'SPMG1(2)' \
  --stim-am1 resid_age ${tdir}/stimam_resid_age.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --glt shape 'SYM: +0.5*shape_gender +0.5*shape_age' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
