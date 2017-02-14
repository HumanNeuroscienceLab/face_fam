
#!/usr/bin/env bash

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 12_trait_subj_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=2
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/150_face_basics_unfam/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/face_basics_unfam"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/traits.reml"
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
  --stim-am1 comp01 ${tdir}/stimam_trait_comp01.txt 'SPMG1(2)' \
  --stim-am1 comp02 ${tdir}/stimam_trait_comp02.txt 'SPMG1(2)' \
  --stim-am1 comp03 ${tdir}/stimam_trait_comp03.txt 'SPMG1(2)' \
  --stim-am1 comp04 ${tdir}/stimam_trait_comp04.txt 'SPMG1(2)' \
  --stim-am1 comp05 ${tdir}/stimam_trait_comp05.txt 'SPMG1(2)' \
  --stim-am1 comp06 ${tdir}/stimam_trait_comp06.txt 'SPMG1(2)' \
  --stim-am1 comp07 ${tdir}/stimam_trait_comp07.txt 'SPMG1(2)' \
  --stim-am1 comp08 ${tdir}/stimam_trait_comp08.txt 'SPMG1(2)' \
  --stim-am1 comp09 ${tdir}/stimam_trait_comp09.txt 'SPMG1(2)' \
  --stim-am1 comp10 ${tdir}/stimam_trait_comp10.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
