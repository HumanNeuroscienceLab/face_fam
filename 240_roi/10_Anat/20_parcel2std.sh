#!/usr/bin/env bash

# This will register a bunch of parcels to 2mm standard
# and also it will create a version that are dilated

function run() {
  echo "$@"
  eval "$@"
  echo $?
}


function parcel_reg2std {
  subject=$1
  
  anatdir="/data1/famface01/analysis/preprocessed/${subject}/anat"
  regdir="${anatdir}/reg"
  atlasdir="${anatdir}/atlases/aparc"
  atlasdir2="${anatdir}/atlases/aparc_std"
  mkdir ${atlasdir2} 2> /dev/null
  
  fns=$( cd ${atlasdir}; ls *nii.gz )
  
  for fn in ${fns}; do
    run gen_applywarp.sh -i ${atlasdir}/${fn} -r ${regdir} -w 'highres-to-standard' -o ${atlasdir2}/${fn} -t nn
  done  
}

function parcel_reg2std_DKT {
  subject=$1
  
  anatdir="/data1/famface01/analysis/preprocessed/${subject}/anat"
  regdir="${anatdir}/reg"
  atlasdir="${anatdir}/atlases/aparc_DKTatlas40"
  atlasdir2="${anatdir}/atlases/aparc_DKTatlas40_std"
  mkdir ${atlasdir2} 2> /dev/null
  
  fns=$( cd ${atlasdir}; ls *nii.gz )
  
  for fn in ${fns}; do
    run gen_applywarp.sh -i ${atlasdir}/${fn} -r ${regdir} -w 'highres-to-standard' -o ${atlasdir2}/${fn} -t nn
  done  
}

function parcel_reg2std_subcort {
  subject=$1
  
  anatdir="/data1/famface01/analysis/preprocessed/${subject}/anat"
  regdir="${anatdir}/reg"
  atlasdir="${anatdir}/segment/aseg"
  atlasdir2="${anatdir}/segment/aseg_std"
  mkdir ${atlasdir2} 2> /dev/null
  
  fns=$( cd ${atlasdir}; ls *nii.gz )
  
  for fn in ${fns}; do
    run gen_applywarp.sh -i ${atlasdir}/${fn} -r ${regdir} -w 'highres-to-standard' -o ${atlasdir2}/${fn} -t nn
  done    
}


subjects='sub01 sub02 sub03 sub04 sub05 sub06'

for subject in ${subjects}; do
  parcel_reg2std ${subject}
  parcel_reg2std_DKT ${subject}
  parcel_reg2std_subcort ${subject}
done

