#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub02 sub03 sub04 sub05 sub06 ) # see 10a_sub01_fix.sh
# parallel --no-notice -j $njobs --eta bash 10a_basic_subj_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=4
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/300_fampics/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/fampics"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/fampics_main_effect.reml"
rm -r ${outdir} 2> /dev/null

echo ${outdir}

${tog}/task_analysis.rb -i ${fundir}/fam_pics/filtered_func_run??.nii.gz \
  -m ${fundir}/mask.nii.gz \
  -b ${fundir}/mean_func.nii.gz \
  --local \
  --output ${outdir} \
  --tr 1 \
  --polort 0 \
  --motion ${tdir}/motion.1D \
  --stim faces ${tdir}/stim_faces.txt 'SPMG1(1)' \
  --stim incorrect ${tdir}/stim_incorrect.txt 'SPMG1(1)' \
  --stim noresp ${tdir}/stim_noresp.txt 'SPMG1(1)' \
  --stim-am1 rt_amp ${tdir}/stimam_rt.txt 'SPMG1(1)' \
  --stim-am1 rt_dur ${tdir}/stimdur_rt.txt 'dmBLOCK' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
