#!/usr/bin/env bash

# See 00_setup_timing_all.ipynb to make the timing files

# njobs=6
# subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
# parallel --no-notice -j $njobs --eta bash 10a_basic_subj_afni.sh {} ::: ${subjects[@]}

tog=/home/zshehzad/guntherxr/bin

# paths
base='/data1/famface01'

nthreads=2
subj=$1
echo $subj

fundir="${base}/analysis/preprocessed/${subj}/func"
tdir="${base}/command/misc/face_representations/300_task_activity/400_imagine/timings/${subj}"
outbase="${base}/analysis/task_activity/${subj}/imagine"
mkdir ${outbase} 2> /dev/null
outdir="${outbase}/imagine_main+identity.reml"
rm -r ${outdir} 2> /dev/null

echo ${outdir}

${tog}/task_analysis.rb -i ${fundir}/fam_pics/filtered_func_run??.nii.gz \
  -m ${fundir}/mask.nii.gz \
  -b ${fundir}/mean_func.nii.gz \
  --output ${outdir} \
  --tr 1 \
  --polort 0 \
  --motion ${tdir}/motion.1D \
  --stim Angelina_Jolie ${tdir}/stim_celeb_Angelina_Jolie.txt 'SPMG2(12)' \
  --stim Brad_Pitt ${tdir}/stim_celeb_Brad_Pitt.txt 'SPMG2(12)' \
  --stim Jennifer_Aniston ${tdir}/stim_celeb_Jennifer_Aniston.txt 'SPMG2(12)' \
  --stim Johnny_Depp ${tdir}/stim_celeb_Johnny_Depp.txt 'SPMG2(12)' \
  --stim Julia_Roberts ${tdir}/stim_celeb_Julia_Roberts.txt 'SPMG2(12)' \
  --stim Justin_Timberlake ${tdir}/stim_celeb_Justin_Timberlake.txt 'SPMG2(12)' \
  --stim Oprah_Winfrey ${tdir}/stim_celeb_Oprah_Winfrey.txt 'SPMG2(12)' \
  --stim Will_Smith ${tdir}/stim_celeb_Will_Smith.txt 'SPMG2(12)' \
  --glt imagine 'SYM: +0.125*Angelina_Jolie +0.125*Brad_Pitt +0.125*Jennifer_Aniston +0.125*Johnny_Depp +0.125*Julia_Roberts +0.125*Justin_Timberlake +0.125*Oprah_Winfrey +0.125*Will_Smith' \
  --glt gender 'SYM: +0.25*Angelina_Jolie -0.25*Brad_Pitt +0.25*Jennifer_Aniston -0.25*Johnny_Depp +0.25*Julia_Roberts -0.25*Justin_Timberlake +0.25*Oprah_Winfrey -0.25*Will_Smith' \
  --regdir ${fundir}/reg \
  --tostandard \
  --threads ${nthreads}
