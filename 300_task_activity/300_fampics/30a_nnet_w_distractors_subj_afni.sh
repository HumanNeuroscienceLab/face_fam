#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 22a_nnet_deviation_subj_afni.sh {} ::: ${subjects[@]}

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
outdir="${outbase}/fampics_nnet_deviate_w_distractor.reml"
rm -r ${outdir} 2> /dev/null

echo ${outdir}

if [[ "$subj" -eq "sub01"  ]]; then
  ${tog}/task_analysis.rb -i ${fundir}/fam_pics/filtered_func_run??.nii.gz \
    -m ${fundir}/mask.nii.gz \
    -b ${fundir}/mean_func.nii.gz \
    --local \
    --output ${outdir} \
    --tr 1 \
    --polort 0 \
    --motion ${tdir}/motion.1D \
    --stim faces2 ${tdir}/stim_faces_distractor.txt 'SPMG1(1)' \
    --stim-am1 avg_diff_all2 ${tdir}/stimam_nnet_avg_deviation_distractor.txt 'SPMG1(1)' \
    --stim faces ${tdir}/stim_faces.txt 'SPMG1(1)' \
    --stim incorrect ${tdir}/stim_incorrect.txt 'SPMG1(1)' \
    --stim-am1 rt_amp ${tdir}/stimam_rt.txt 'SPMG1(1)' \
    --stim-am1 avg_diff_all ${tdir}/stimam_nnet_avg_deviation_all.txt 'SPMG1(1)' \
    --stim-am1 avg_diff_person ${tdir}/stimam_nnet_avg_deviation_person_resids.txt 'SPMG1(1)' \
    --glt faces_tot 'SYM: +0.5*faces +0.5*faces2' \
    --glt avg_diff_tot 'SYM: +0.5*avg_diff_all +0.5*avg_diff_all2' \
    --regdir ${fundir}/reg \
    --tostandard \
    --threads ${nthreads}
  # --stim noresp ${tdir}/stim_noresp.txt 'SPMG1(1)' \
else
  ${tog}/task_analysis.rb -i ${fundir}/fam_pics/filtered_func_run??.nii.gz \
    -m ${fundir}/mask.nii.gz \
    -b ${fundir}/mean_func.nii.gz \
    --local \
    --output ${outdir} \
    --tr 1 \
    --polort 0 \
    --motion ${tdir}/motion.1D \
    --stim faces2 ${tdir}/stim_faces_distractor.txt 'SPMG1(1)' \
    --stim-am1 avg_diff_all2 ${tdir}/stimam_nnet_avg_deviation_distractor.txt 'SPMG1(1)' \
    --stim faces ${tdir}/stim_faces.txt 'SPMG1(1)' \
    --stim incorrect ${tdir}/stim_incorrect.txt 'SPMG1(1)' \
    --stim noresp ${tdir}/stim_noresp.txt 'SPMG1(1)' \
    --stim-am1 rt_amp ${tdir}/stimam_rt.txt 'SPMG1(1)' \
    --stim-am1 avg_diff_all ${tdir}/stimam_nnet_avg_deviation_all.txt 'SPMG1(1)' \
    --stim-am1 avg_diff_person ${tdir}/stimam_nnet_avg_deviation_person_resids.txt 'SPMG1(1)' \
    --glt faces_tot 'SYM: +0.5*faces +0.5*faces2' \
    --glt avg_diff_tot 'SYM: +0.5*avg_diff_all +0.5*avg_diff_all2' \
    --regdir ${fundir}/reg \
    --tostandard \
    --threads ${nthreads}
fi
