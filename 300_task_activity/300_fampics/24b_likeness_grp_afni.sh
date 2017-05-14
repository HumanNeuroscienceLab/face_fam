#!/usr/bin/env bash

###
# SUPPORT FUNCIONS
###

run() {
  echo "$@"
  eval "$@"
  return $?
}

die() {
  echo "$@"
  exit 2
}


###
# ARGS
###

njobs=30
subjects=( sub01 sub02 sub03 sub04 sub05 sub06 )
outdir="/data1/famface01/analysis/task_activity/group/fampics/fampics_likeness.mema"
taskdir_template="/data1/famface01/analysis/task_activity/__SUBJECT__/fampics/fampics_likeness.reml"


###
# SETUP
###

export OMP_NUM_THREADS=$njobs

scriptdir=$(pwd)

# Make the output directory
echo
echo "creating output directory"
[ -e ${outdir} ] && run "mv ${outdir} ${outdir}+"
mkdir -p ${outdir} 2> /dev/null
echo "cd ${outdir}"
cd ${outdir}

# Make subjects directory
mkdir subjects 2> /dev/null

echo
echo "copying or linking input subject data"
for subject in ${subjects[@]}; do
  taskdir=$( echo "${taskdir_template}" | sed -e s/__SUBJECT__/${subject}/g )
  [ ! -e ${taskdir} ] && die "Task directory (${taskdir}) does not exist"
  
  rdir="${taskdir}/reg_standard"
  run "ln -sf ${rdir}/stats subjects/${subject}_stats"
  run "ln -sf ${rdir}/mask.nii.gz subjects/${subject}_mask.nii.gz"
  run "tail -1 ${taskdir}/blur_est.1D > subjects/${subject}_blur_errts.1D"
done

echo
echo "creating group mask + grey-matter mask"
run "3dMean -mask_inter -prefix mask_func.nii.gz subjects/*_mask.nii.gz"
run "3dcalc -a /data1/famface01/analysis/roi/10_Anat/fs_grey_mask_allsubs_dil1-1.nii.gz -expr 'a' -prefix mask_grey.nii.gz"
#run "3dcalc -a $FSLDIR/data/standard/tissuepriors/2mm/gray_10perc.nii.gz -expr 'a' -prefix mask_grey.nii.gz"
run "3dcalc -a mask_func.nii.gz -b mask_grey.nii.gz -expr 'step(a)*step(b)' -prefix mask.nii.gz"

echo
echo "gathering contrasts to run (just the different coef names)"
contrasts=( $( cd subjects/${subjects[0]}_stats; ls coef_*.nii.gz | sed s/coef_//g | sed s/.nii.gz//g ) )
echo "=> ${contrasts[@]}"
contrasts=( likeness )
echo "selecting ${contrasts[@]}"

###
# ANALYSIS
###

echo
echo "looping through contrasts to run 3dMEMA"
# note: use cio option so I can use nifti files
for cname in ${contrasts[@]}; do
  echo
  echo "========"
  echo "${cname}"
  
  cmd="3dMEMA \
          -mask mask.nii.gz
          -prefix ${cname} \
          -jobs ${njobs} \
          -missing_data 0 \
          -HKtest         \
          -model_outliers \
          -cio \
          -residual_Z \
          -verb 1 \
          -set  ${cname}"
  for subject in ${subjects[@]}; do
    cmd+=" ${subject} subjects/${subject}_stats/coef_${cname}.nii.gz subjects/${subject}_stats/tstat_${cname}.nii.gz"
  done
  run $cmd
  
  echo "refit"
  run "3drefit -view tlrc -space MNI ${cname}_ICC+orig"
  run "3drefit -view tlrc -space MNI ${cname}_resZ+orig"
  run "3drefit -view tlrc -space MNI ${cname}+orig"
  
  echo "to nifti"
  run "3dAFNItoNIFTI -overwrite -prefix ${cname}_ICC.nii.gz ${cname}_ICC+tlrc"
  run "3dAFNItoNIFTI -overwrite -prefix ${cname}_resZ.nii.gz ${cname}_resZ+tlrc"
  run "3dAFNItoNIFTI -overwrite -prefix ${cname}.nii.gz ${cname}+tlrc"
  
  echo "extract tstats"
  run "3dcalc -a ${cname}.nii.gz'[1]' -expr a -prefix tstats_${cname}.nii.gz"
  
  echo "tstat => zstat"
  run "3dcalc -a tstats_${cname}.nii.gz -expr 'fitt_t2z(a,${#subjects[@]})' -prefix zstats_${cname}.nii.gz"
  
  echo "mask"
  run "3dcalc -a tstats_${cname}.nii.gz -expr 'step(abs(a))' -prefix ${cname}_mask.nii.gz"
  
  echo "========"
  echo
done

echo
echo "standard underlays"
reses="0.5 1 2"
for res in ${reses}; do
  run "ln -sf ${FSLDIR}/data/standard/MNI152_T1_${res}mm.nii.gz standard_${res}mm.nii.gz"
  run "ln -sf ${FSLDIR}/data/standard/MNI152_T1_${res}mm_brain.nii.gz standard_brain_${res}mm.nii.gz"
done
run "rm standard_brain_0.5mm.nii.gz"


###
# CLUSTER CORRECT (various methods)
###

echo
echo "easythresh"
run "mkdir easythresh 2> /dev/null"
run "cd easythresh"
for cname in ${contrasts[@]}; do
  echo "contrast: ${cname}"
  ## positive
  run "easythresh ../zstats_${cname}.nii.gz ../mask.nii.gz 1.96 0.05 ../standard_2mm.nii.gz zstat_${cname}"
  run "easythresh ../zstats_${cname}.nii.gz ../mask.nii.gz 1.645 0.1 ../standard_2mm.nii.gz liberal_zstat_${cname}"
  ## negative
  run "fslmaths ../zstats_${cname}.nii.gz -mul -1 ztmp_flip_zstats_${cname}.nii.gz"
  run "easythresh ztmp_flip_zstats_${cname}.nii.gz ../mask.nii.gz 1.96 0.05 ../standard_2mm.nii.gz flipped_zstat_${cname}"
  run "easythresh ztmp_flip_zstats_${cname}.nii.gz ../mask.nii.gz 1.645 0.1 ../standard_2mm.nii.gz flipped_liberal_zstat_${cname}"
  ## combine
  run "3dcalc -overwrite -a thresh_zstat_${cname}.nii.gz -b thresh_flipped_zstat_${cname}.nii.gz -expr 'a-b' -prefix combined_thresh_zstat_${cname}.nii.gz"  
  run "3dcalc -overwrite -a thresh_liberal_zstat_${cname}.nii.gz -b thresh_flipped_zstat_${cname}.nii.gz -expr 'a-b' -prefix combined_thresh_liberal_zstat_${cname}.nii.gz"
  ## clean
  run "rm -f *flipped*"
  run "rm -f ztmp_flip_zstats_${cname}.nii.gz"
  echo
done
run "cd -"

echo
echo "now do FDR"
mkdir fdr 2> /dev/null
for cname in ${contrasts[@]}; do
  run "fslmaths zstats_${cname}.nii.gz -abs -ztop pvals_${cname}.nii.gz"
  run "fdr -i pvals_${cname}.nii.gz -m ${cname}_mask.nii.gz -q 0.05 --othresh=fdr_mask_${cname}.nii.gz"
  run "3dcalc -overwrite -a zstats_${cname}.nii.gz -b fdr_mask_${cname}.nii.gz -expr 'a*step(1-b)' -prefix fdr/fdr_zstats_${cname}.nii.gz"
  run "rm fdr_mask_${cname}.nii.gz pvals_${cname}.nii.gz"
done

echo
echo "run cluster threshold simulations"
# gather all the blurs
run "cat subjects/*_blur_errts.1D > tmp_blur_errts.1D"
# compute average blur and append
blurs=( `3dTstat -mean -prefix - tmp_blur_errts.1D\'` )
echo "average errts blurs: ${blurs[@]}"
echo "${blurs[@]}" > blur_est.1D
# clustsim
fxyz=( `tail -1 blur_est.1D` )
run "3dClustSim -both -NN 123 -mask mask.nii.gz \
           -fwhmxyz ${fxyz[@]:0:3} -prefix ClustSim"
echo "apply cluster results to each output file"
for cname in ${contrasts[@]}; do
  run "3drefit -atrstring AFNI_CLUSTSIM_MASK file:ClustSim.mask                \
          -atrstring AFNI_CLUSTSIM_NN1  file:ClustSim.NN1.niml            \
          -atrstring AFNI_CLUSTSIM_NN2  file:ClustSim.NN2.niml            \
          -atrstring AFNI_CLUSTSIM_NN3  file:ClustSim.NN3.niml            \
          ${outdir}/${cname}+tlrc"
done

echo
echo "do the same but on nifti (+ have a liberal threshold)"
run "mkdir afnithresh"
for cname in ${contrasts[@]}; do
  echo "contrast: ${cname}"
  #run "rm -f thresh_zstats_${cname}.nii.gz thresh_liberal_zstats_${cname}.nii.gz"
  run "rm -f afnithresh/thresh_zstats_${cname}.nii.gz afnithresh/thresh_liberal_zstats_${cname}.nii.gz"
  run "${scriptdir}/./apply_clustsim.R ClustSim.NN3_2sided 0.05 0.05 zstats_${cname}.nii.gz afnithresh/thresh_zstats_${cname}.nii.gz"
  run "${scriptdir}/./apply_clustsim.R ClustSim.NN3_2sided 0.1 0.1 zstats_${cname}.nii.gz afnithresh/thresh_liberal_zstats_${cname}.nii.gz"
  echo
done

echo
echo "combine easythresh together"
run "${scriptdir}/./combine_zstats.sh ${outdir}"
run "${scriptdir}/./combine_zstats_afni.sh ${outdir}"

echo "END"

