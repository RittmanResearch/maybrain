# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 21:59:59 2015
This module deals with the braineac gene expression data.
@author: tim
"""

import maybrain.brainObjs as mbo

import nibabel as nb
import numpy as np


import csv
from os import path,rename
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from glob import glob

import rpy2.robjects as ro

class braineac:
    def __init__(self,brain, AALrois="ROI_MNI_V4.txt"):
        
        # get dictionary of nodes ready
        f = open(AALrois, "r")
        lines = f.readlines()
        
        # dictionary of regions defined by braineac and key words in AAL
        self.brD = {"Cerebelum":"CRBL",
                  "Vermis":"CRBL",
                  "Frontal":"FCTX",
                  "Hippocampus":"HIPP",
                  "Occipital":"OCTX",
                  "Putamen":"PUTM",
                  "Temporal":"TCTX",
                  "Thalamus":"THAL"}

        self.ROIs= {v:[] for v in brD.values()}
        
        for l in lines:
            bits = l.rstrip('\n').split()
            for k,v in brD.iteritems():
                if k in bits[1]:
                    self.ROIs[v].append(bits[2])
        f.close()

    def makeTemplate(self, AAL="ROI_MNI_V4.nii.gz", outFile="AAL_braineac.nii",
                     logName="braineacRegionLog.txt"):
        '''
        Create an AAL based template of the braineac regions
        '''
        # load data from aal atlas
        nbaal = nb.load(AAL)
        aal = nbaal.get_data()
        Header = nbaal.get_header() # get header
        
        # Create new template image with new parcel values
        out = np.zeros(aal.shape, dtype="int")  # create output matrix
        
        # create dictionary of which the values are parcels in ALL atlas lie within the keys of braineac regions
        # create log file with names and regions
        log = open(logName, "w")
        self.regValDict = {}
        for n,k in enumerate(self.ROIs.keys()):
            for x in ROIs[k]:
                out[np.where(aal==float(x))] = n+1
            log.writelines(' ' .join([str(v) for v in [k,n+1]])+'\n')
            self.regValDict[n+1] = k
        
        # Save new brain template
        outNii = nb.Nifti1Image(out, nbaal.get_affine(), header=Header)
        nb.save(outNii, outFile)
        
        log.close() # close log file





data = open('t3723687_rda/3723687.rda').read()


# collate metrics within parcels
mFile = "Control/brain_d2_degree_wt.txt"

m = open(mFile, "r")
mDict = 