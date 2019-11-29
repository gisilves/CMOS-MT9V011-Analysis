#!/usr/bin/env python
# coding: utf-8
"""CMOS MT9V011 Cluster Detection"""

#Import libraries
import os
import sys
import fnmatch
import pickle
import pandas as pd
from math import sqrt
from tqdm import tqdm
from skimage.feature import blob_log

import numpy as np
from numpy import interp

#Files directory
#DATADIR = "/nfs/NASPG/BTData/Jul2012_CMOS_data/MT9V011_Firenze_2012_07_13/TEST/"
DATADIR = sys.argv[1]

PEDS_VALUE = np.zeros(307200, dtype="int")
PEDS_EVTS = np.zeros(307200, dtype="int")
PEDS = np.zeros(307200, dtype="int")

CLUSTERS = pd.DataFrame({
    "Max": [],
    "Sum": [],
    "Size": [],
    "row": [],
    "column": []
})

#Calculate or (load) Pedestals
if not os.path.isfile(DATADIR + "pedestals.npy"):
    THRESH = 100
    for _filename in tqdm(sorted(fnmatch.filter(os.listdir(DATADIR),
                                                '*.txt'))):
        print("Processing " + _filename + " for pedestals computation")
        im = np.loadtxt(DATADIR + _filename, dtype="int")
        if im.shape == PEDS.shape:
            mask = im < THRESH  # mask of pixels under threshold (good for pedestals)
            im[im >= THRESH] = 0  # to prevent adding events with signal
            PEDS_VALUE += im
            PEDS_EVTS += mask.astype(int)
        else:
            continue
    PEDS_EVTS[PEDS_EVTS == 0] = 1  #Mask pixels with 0 good events
    PEDS = PEDS_VALUE / PEDS_EVTS
    PEDS.dump(DATADIR + "pedestals.npy")
else:
    PEDS = np.load(DATADIR + "pedestals.npy", allow_pickle="True")

#Find clusters in every frame
for _filename in tqdm(sorted(fnmatch.filter(os.listdir(DATADIR), '*.txt'))):
    print("Processing " + _filename)
    frame = np.loadtxt(DATADIR + _filename, dtype="int")

    if frame.shape == PEDS.shape:
        frame = frame - PEDS
    else:
        continue

    frame = frame.reshape(480, 640)
    _temp = interp(frame, [0, 1024], [0, 1])

    #Find "blobs"
    blobs_log = blob_log(_temp, max_sigma=50, num_sigma=50, threshold=.05)
    blobs_log[:, 2] = blobs_log[:, 2] * 2 * sqrt(2)

    for _x, _y, _r in blobs_log:
        _min_y = (_y - _r).astype(int)
        _max_y = (_y + _r).astype(int)
        _min_x = (_x - _r).astype(int)
        _max_x = (_x + _r).astype(int)
        _subim = frame[_min_x:_max_x, _min_y:_max_y]
        CLUSTERS = CLUSTERS.append(
            {
                'Max': np.max(_subim),
                'Sum': np.sum(_subim),
                'Size': _subim.size,
                'row': _x,
                'column': _y
            },
            ignore_index=True)
    CLUSTERS.to_feather(DATADIR + "clusters.feather")