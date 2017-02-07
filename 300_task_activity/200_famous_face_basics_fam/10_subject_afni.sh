#!/usr/bin/env bash

# njobs=5
# subjects=( sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 10_subject_afni.sh {} ::: ${subjects[@]}


tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=4
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/200_famous_face_basics_fam/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/famous_faces"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/basics.reml"
rm -r ${outdir} 2> /dev/null

echo ${outdir}

${tog}/task_analysis.rb -i ${fundir}/fam_vids/filtered_func_run??.nii.gz \
  -m ${fundir}/mask.nii.gz \
  -b ${fundir}/mean_func.nii.gz \
  --output ${outdir} \
  --global \
  --tr 1 \
  --polort 0 \
  --motion ${tdir}/motion.1D \
  --stim famous ${tdir}/stim_famous.txt 'SPMG1(2)' \
  --stim nonfamous ${tdir}/stim_nonfamous.txt 'SPMG1(2)' \
  --stim quests ${tdir}/stim_questions.txt 'SPMG1(4)' \
  --glt faces 'SYM: +famous +nonfamous' \
  --glt fam_gt_nonfam 'SYM: +famous -nonfamous' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
