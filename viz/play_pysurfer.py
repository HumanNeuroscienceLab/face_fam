import os
import numpy as np
from surfer import Brain, project_volume_data

brain = Brain("mni152", "rh", "inflated")
brain.add_annotation("aparc")
brain.save_montage('z_mni152.png', order=['lat','ven'])

brain = Brain("fsaverage", "rh", "inflated")
brain.add_annotation("aparc")
brain.save_montage('z_fsaverage.png', order=['lat','ven'])

brain = Brain("cvs_avg35_inMNI152", "rh", "inflated")
brain.add_annotation("aparc")
brain.save_montage('z_cvs.png', order=['lat','ven'])

brain = Brain("mni152", "rh", "midthickness")
brain.add_annotation("aparc")
brain.save_montage('z_mni152.png', order=['lat','ven'])


base = "/Volumes/hnl17"
asap_base = "%s/mnt/nfs/share/rois/asap/ASAP_maps" % base

prob_file1 = "%s/facescene_pmap_N124_stat3.nii.gz" % asap_base
prob_file2 = "%s/facescene_pmap_N124_stat4.nii.gz" % asap_base


prob1 = project_volume_data(prob_file1, "rh", subject_id="mni152", smooth_fwhm=2)
prob2 = project_volume_data(prob_file2, "rh", subject_id="mni152", smooth_fwhm=2)

brain = Brain("mni152", "rh", "inflated")
brain.add_annotation("aparc")

brain.add_overlay(prob1, min=0.2, max=0.6)