#!/usr/bin/env bash

function parcel_reg2std {
  subject=$1
  
  anatdir="/data1/famface01/analysis/preprocessed/${subject}/anat"
  regdir="${anatdir}/reg"
  atlasdir="${anatdir}/atlases"
  mkdir ${atlasdir2} 2> /dev/null
  
  gen_applywarp.sh -i ${atlasdir}/aparc+aseg_2.nii.gz -r ${regdir} -w 'highres-to-standard' -o ${atlasdir}/aparc+aseg_2_to_std.nii.gz -t nn
  
  gen_applywarp.sh -i ${atlasdir}/aseg.nii.gz -r ${regdir} -w 'highres-to-standard' -o ${atlasdir}/aseg_to_std.nii.gz -t nn
}

subjects='sub01 sub02 sub03 sub04 sub05 sub06'

for subject in ${subjects}; do
  parcel_reg2std ${subject}
done

