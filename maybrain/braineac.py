# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 21:59:59 2015
This module deals with the braineac gene expression data available from http://www.braineac.org/

It also requires an AAL atlas and AAL regions of interest file (contact me on tr332@medschl.cam.ac.uk if you need either of these)
@author: tim
"""

import maybrain.brainObj as mbo

import nibabel as nb
import numpy as np
from os import path

class braineac:
    """
    A class the defines a braineac object. This reads braineac data and defines
    regions based on their definitions in the AAL atlas. It can then extract data
    from a maybrain produced set of metrics and aggregate node-wise measures over
    the braineac regions.
    
    Don't forget to set the path to your braineac diretory (braineacDir) which
    should contain the downloaded braineac data.
    """
    def __init__(self,AALrois="ROI_MNI_V4.txt", braineacDir = "/home/tim/Documents/braineac"):
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

        self.ROIs= {v:[] for v in self.brD.values()}
        self.braineacDir = braineacDir
        
        for l in lines:
            bits = l.rstrip('\n').split()
            for k,v in self.brD.iteritems():
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
        self.aff = nbaal.get_affine()
        self.Header = nbaal.get_header() # get header
        
        # Create new template image with new parcel values
        self.out = np.zeros(aal.shape, dtype="int")  # create output matrix
        
        # create dictionary of which the values are parcels in AAL atlas lie within the keys of braineac regions
        # create log file with names and regions
        log = open(logName, "w")
        self.regValDict = {}
        for n,k in enumerate(self.ROIs.keys()):
            for x in self.ROIs[k]:
                self.out[np.where(aal==float(x))] = n+1
            log.writelines(' ' .join([str(v) for v in [k,n+1]])+'\n')
            self.regValDict[n+1] = k
        
        # Save new brain template
        outNii = nb.Nifti1Image(self.out, self.aff, header=self.Header)
        nb.save(outNii, outFile)
        
        log.close() # close log file

    def comparison(self, assMat, delim=" ", spatialFile="parcel_500.txt", parcellation="parcels/parcel_500.nii",
                   outFile="braineacNodes.txt", midLine=44.5):
        """
        now let's match up the imaging data to the braineac data        
        to create a brain object and import an association matrix.
        """
        self.a = mbo.brainObj()
        self.a.importAdjFile(assMat, delimiter=delim)
        
        # get the spatial information and parcellation scheme
        self.a.importSpatialInfo(spatialFile)
        self.a.importISO(parcellation)
        self.a.parcels(self.a.G.nodes())
        
        # set up output file
        outFile = open(outFile, "w")
        
        for node in self.a.G.nodes():
            self.a.G.node['braineac'] = None
        
        self.nodeDict = {n:{"L":[], "R":[]} for n in range(1,8)}
        for n in self.nodeDict.keys():
            tmp = self.a.parcelList[np.where(self.out==n)]
            nodeList = [int(v) for v in np.unique(tmp.data.flatten()) if not v==0]
            
            # check if 50% of parcel is in the braineac region
            for node in nodeList:
                rat = float(len(self.a.parcelList[np.where(self.a.parcelList==node) and np.where(self.a.parcelList==n)])) / float(len(self.a.parcelList[np.where(self.a.parcelList==node)]))
                if rat < 0.5:
                    print str(node)+" less than 50% in "+str(n)
                    
                elif self.a.G.node[node-1]["xyz"][0] < midLine:
                    self.nodeDict[n]["L"].append(node-1)
                    self.a.G.node[node-1]['braineac'] = self.regValDict[n]+"L"
                else:
                    self.nodeDict[n]["R"].append(node-1)
                    self.a.G.node[node-1]['braineac'] = self.regValDict[n]+"R"
                    
            # write the output in to a file
            outFile.writelines(self.regValDict[n]+"L" + ' ' + ' '.join([str(v) for v in self.nodeDict[n]["L"]])+'\n')
            outFile.writelines(self.regValDict[n]+"R" + ' ' + ' '.join([str(v) for v in self.nodeDict[n]["R"]])+'\n')
            
        outFile.close()

    def createParcellation(self, outFile="AAL_braineac_parcels.nii"):
        """
        Create a parcellation based on groups of imaging data assigned a single value based on the braineac cluster
        they fall in to. The 'comparison' function must have been run first.
        """
        # create parcellation image
        isoPC = np.zeros(self.a.iso.shape, dtype="int")
        for pc in self.nodeDict.keys():
            for n,s in enumerate(["L", "R"]):
                fillVal = float(str(n+1)+str(pc))
                for ar in self.nodeDict[pc][s]:
                    isoPC[np.where(self.a.iso==ar+1)] = fillVal
                    
                    
        pcNii = nb.Nifti1Image(isoPC, self.aff, header=self.Header)
        nb.save(pcNii, outFile)

    def braineacData(self, probeList, naVal="NA"):
        """ 
        Get braineac data for the specified list of probes.
        """
        # create a brain object
        self.b = mbo.brainObj()
        
        # the list of brineac nodes
        nodeList = [v for v in self.regValDict.keys()]
        nodeList.sort()
        
        # iterate through the braineac regions
        for n,k in enumerate(nodeList):
            # open the data file for the appropriate regions
            inFile = "expr_"+self.regValDict[k]+'.txt'
            f = open(path.join(self.braineacDir,inFile), "r")
            lines = f.readlines()
            
            # have a look through each line of the data file to look for the appropriate probe data            
            for line in lines:
                bits = line.split()
                if bits[0] in probeList:
                    pMean = np.mean([float(v) for v in bits[1:] if not v==naVal])
                    self.b.G.add_node(n, attr_dict={bits[0]:pMean})
                
                
            
            
            

         
