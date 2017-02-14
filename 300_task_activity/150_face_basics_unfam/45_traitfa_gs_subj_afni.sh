
#!/usr/bin/env bash

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 45_traitfa_gs_subj_afni.sh {} ::: ${subjects[@]}

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
outdir="${outbase}/traitsfa_givenshape.reml"
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
  --stim-am1 shape_trait1 ${tdir}/stimam_shape.trait1.txt 'SPMG1(2)' \
  --stim-am1 shape_trait2 ${tdir}/stimam_shape.trait2.txt 'SPMG1(2)' \
  --stim-am1 shape_trait3 ${tdir}/stimam_shape.trait3.txt 'SPMG1(2)' \
  --stim-am1 shape_trait4 ${tdir}/stimam_shape.trait4.txt 'SPMG1(2)' \
  --stim-am1 shape_trait5 ${tdir}/stimam_shape.trait5.txt 'SPMG1(2)' \
  --stim-am1 shape_trait6 ${tdir}/stimam_shape.trait6.txt 'SPMG1(2)' \
  --stim-am1 resid_trait1 ${tdir}/stimam_resid.trait1.txt 'SPMG1(2)' \
  --stim-am1 resid_trait2 ${tdir}/stimam_resid.trait2.txt 'SPMG1(2)' \
  --stim-am1 resid_trait3 ${tdir}/stimam_resid.trait3.txt 'SPMG1(2)' \
  --stim-am1 resid_trait4 ${tdir}/stimam_resid.trait4.txt 'SPMG1(2)' \
  --stim-am1 resid_trait5 ${tdir}/stimam_resid.trait5.txt 'SPMG1(2)' \
  --stim-am1 resid_trait6 ${tdir}/stimam_resid.trait6.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
