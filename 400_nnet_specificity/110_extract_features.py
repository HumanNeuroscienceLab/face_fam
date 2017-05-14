#!/usr/bin/env python

import numpy as np
import lutorpy as lua
import scipy.io as sio # for saving
from glob import glob
from tqdm import tqdm
import os

# rpy2
#import rpy2
#import rpy2.robjects as robjects
#r = robjects.r
#import rpy2.robjects.numpy2ri
#rpy2.robjects.numpy2ri.activate()
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

require('torch')
require('nn')
require('dpnn')
require('image')
require('paths')

torch.setdefaulttensortype('torch.FloatTensor')

#require('cunn')
#require('cudnn')
#require('cutorch')
#torch.setdefaulttensortype('torch.CudaTensor')

modelPath = "/home/zshehzad/Downloads/openface/models/openface/nn4.small2.v1.t7"
imgDim    = 96
imgDir    = "/home/zshehzad/Downloads/tmp_sf_nn/single_frames"
imgFiles  = sorted(glob("%s/*.png" % imgDir))
outDir    = "/data1/famface01/analysis/misc/openface/layer_features"
if not os.path.exists(outDir): os.mkdir(outDir)

# Load neural network
net = torch.load(modelPath)
net._evaluate()
print(net)
nlayers = net._size()

# Loop through each image and save the features separately
for imgPath in tqdm(imgFiles):  
  ## load the image
  img = torch.Tensor(1, 3, imgDim, imgDim)
  img_orig = image.load(imgPath, 3)
  img[0] = image.scale(img_orig, imgDim, imgDim)
  _ = net._forward(img)
  
  # Loop through each layer and get features into list
  dlayers = {}
  for i in range(nlayers):
    ## extract features
    layer = net.modules[i].output
    layer = np.squeeze(layer.asNumpyArray())
    dlayers['layer%02i' % (i+1)] = layer
  ## save
  vidfile = os.path.basename(imgPath).replace(".png", ".mat")
  sio.savemat('%s/%s' % (outDir, vidfile), dlayers)
  
