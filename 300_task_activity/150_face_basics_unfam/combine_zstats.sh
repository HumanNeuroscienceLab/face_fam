#!/usr/bin/env bash

# Setup
grpdir="$1"
thresh=1.96

if [[ $# -ne 1 ]]; then
  echo "usage: $0 grpdir"
  echo "assume thresh of ${thresh}"
  exit 1
fi
if [[ ! -e "$grpdir" ]]; then
  echo "grpdir: $grpdir does not exist"
  exit 2
fi

function run {
  echo "$@"
  eval "$@"
  return $?
}


# Process
cd ${grpdir}/easythresh

## paths
contrasts=( $( ls combined_thresh_zstat_*.nii.gz | sed s/combined_thresh_zstat_//g | sed s/.nii.gz//g ) )
fnames=( $( ls combined_thresh_zstat_*.nii.gz ) )

## remove the threshold for better viewing in afni
tmp_fnames=()
for (( i = 0; i < ${#fnames[@]}; i++ )); do
  tmp_fnames[i]="ZTMP_COMBINE_${fnames[i]}"
  run "3dcalc -a ${fnames[i]} -expr 'step(a)*(a-${thresh}+0.2) - step(-1*a)*((-1*a)-${thresh}+0.2)' -prefix ${tmp_fnames[i]}"
done

## combine
run 3dbucket -overwrite -prefix ../bucket_combined_thresh_zstat -fbuc ${tmp_fnames[@]}

## add labels
cmds=()
for (( i = 0; i < ${#contrasts[@]}; i++ )); do
  cmds[i]="-sublabel $i ${contrasts[i]}"
done
run 3drefit ${cmds[@]} ../bucket_combined_thresh_zstat+tlrc.

## remove
run rm ${tmp_fnames[@]}
