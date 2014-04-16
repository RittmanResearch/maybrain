# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 19:05:26 2014

@author: tim
"""

import mayBrainTools as mbt


import csv
from os import path,rename
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

class allenBrain:
    def __init__(self, allenSubj, assocMat, delim=",",
                 spatialFile="parcel_500_xyz.txt", nodesToExclude=None):
        
        self.subj = allenSubj
        self.fName = "SampleAnnot.csv"
        self.maFile = "MicroarrayExpression.csv"
        
        # set up brain for expression data
        self.a = mbt.brainObj()
        
        # import probe data
        f = open(path.join(self.subj, self.fName), "rb")        
        reader = csv.DictReader(f, delimiter=",", quotechar='"')
        
        self.headers = ['probe']
        for l in reader:
            n = int(l["structure_id"])
            self.a.G.add_node(n)
            self.a.G.node[n] = l
            self.headers.append(l["structure_id"])
        f.close()
        
        # convert location data for Allen brain from MNI space
        for n in self.a.G.nodes():
            x = 45 - (float(self.a.G.node[n]['mni_x'])/2)
            y = 63 + (float(self.a.G.node[n]['mni_y'])/2)
            z = 36 + (float(self.a.G.node[n]['mni_z'])/2)
            self.a.G.node[n]['xyz'] = (x,y,z)
        
        # set up brain with graph properties
        self.c = mbt.brainObj()
        self.c.readAdjFile(assocMat, delimiter=delim,
                           excludedNodes=nodesToExclude)
        self.c.readSpatialInfo("parcel_500_xyz.txt")
        
        # set up dictionary to link nodes from probe data and graph
        nodeDictMRIs = {}
        
        for node in self.c.G.nodes():
            dOther = (None, 999.)
            dOwn = (None, 999.)
        
            for n in self.a.G.nodes():
                d = np.linalg.norm(np.array(self.c.G.node[node]['xyz'] - np.array(self.a.G.node[n]['xyz'])))
                if d < dOther[1]:
                    dOther = (n,d)
            
            for n in [v for v in self.c.G.nodes() if not v==node]:
                d = np.linalg.norm(np.array(self.c.G.node[node]['xyz'] - np.array(self.c.G.node[n]['xyz'])))
                if d < dOwn[1]:
                    dOwn = (n,d)
            nodeDictMRIs[node] = {"allen":dOther, "MRIs":dOwn}
        
        # set up dictionary to link nodes from probe data and graph
        nodeDictAllen = {}
        for node in self.a.G.nodes():
            dOther = (None, 999.)
            dOwn = (None, 999.)
            
            for n in self.c.G.nodes():
                d = np.linalg.norm(np.array(self.a.G.node[node]['xyz'] - np.array(self.c.G.node[n]['xyz'])))
                if d < dOther[1]:
                    dOther = (n,d)
            
            for n in [v for v in self.a.G.nodes() if not v==node]:
                d = np.linalg.norm(np.array(self.a.G.node[node]['xyz'] - np.array(self.a.G.node[n]['xyz'])))
                if d < dOwn[1]:
                    dOwn = (n,d)
            nodeDictAllen[node] = {"allen":dOwn, "MRIs":dOther}
        
        nodePairs = []
        for node in nodeDictMRIs.keys():
            # check if no other MRI nodes closer
            if nodeDictMRIs[node]['allen'][1] < nodeDictMRIs[node]['MRIs'][1]:
                # check if the closest MRI node
                n = nodeDictMRIs[node]['allen'][0]
                if nodeDictAllen[n]['MRIs'][0] == node:
                    self.c.G.node[node]['pair'] = n
                    self.a.G.node[n]['pair'] = node
                    nodePairs.append((node,n))
                else:
                    self.c.G.remove_node(node)
            else:
                self.c.G.remove_node(node)
        
        for node in self.a.G.nodes():
            if not 'pair' in self.a.G.node[node].keys():
                self.a.G.remove_node(node)
                    
    def doPlot(self):
        self.a.importSkull("/usr/share/data/fsl-mni152-templates/MNI152_T1_2mm_brain.nii.gz")
        p = mbt.plotObj()
        p.plotSkull(self.a, contourVals = [3000,9000])
        p.plotBrainNodes(self.c, nodes =  self.c.G.nodes(),
                                     col=(1,0,0), sizeList=5)
        p.plotBrainNodes(self.a, nodes = self.a.G.nodes(),
                                     col=(0,0,1), sizeList=5)
        self.saveFig()
                                     
    def probeData(self, propDict, graphMetric="gm", nodeList=None, plot=False,
                  probeList=[], sigVal=1.0, T=False):
        self.gm=graphMetric
        self.sigVal=sigVal
        
        probes = open(path.join(self.subj,"Probes.csv"), "rb")
        probeReader = csv.DictReader(probes, delimiter=",",
                                     quotechar='"')
        self.probeDict = {l['probe_id']:[l['gene_symbol'], l['gene_name']] for l in probeReader}
        if probeList:
            probeNumbers = []
            for p in probeList:
                probeNumbers.extend([v for v in self.probeDict.keys() if any([p in self.probeDict[v][1], p in self.probeDict[v][0]])])
            print " ".join(["Probe numbers:", ' '.join(probeNumbers)])
        
        else:
            probeNumbers = None
        
        self.outFile = path.join(self.subj, self.gm+'.txt')
        if path.exists(self.outFile):
            rename(self.outFile, self.outFile+'.old')
        f = open(self.outFile,"w")
        f.writelines(','.join(['probe_id', 'gene_name', 'r','p\n']))
        f.close()

        self.propDict = propDict
                                
        if nodeList:
            for node in self.c.G.nodes():
                if not node in nodeList:
                    self.a.G.remove_node(self.c.G.node[node]['pair'])
                    self.c.G.remove_node(node)
        
        # import probe data
        f = open(path.join(self.subj, self.maFile), "rb")
        reader = csv.DictReader(f, delimiter=",",
                                fieldnames=self.headers,
                                quotechar='"')
        
        for l in reader:
            if T:
                self.probeSubT(l, plot, probeNumbers)
            else:
                self.probeSub(l, plot, probeNumbers)
        f.close()
    
#    def brainProbeVis(self):
#        
    
    def probeSub(self, l, plot, probeNumbers=None):
        probe = l['probe']
        if probeNumbers:
            if probe in probeNumbers:
                # assign probe values to sample numbers
                for node in self.a.G.nodes():
                    if l[str(node)]:
                        self.a.G.node[node][probe] = l[str(node)]
                    else:
                        self.a.G.node[node][probe] = None
                
                aa = np.zeros((len(self.a.G.nodes()),3))
                
                for n,node in enumerate(self.a.G.nodes()):
                    if self.propDict[self.a.G.node[node]['pair']]:
                        aa[n,0] = self.a.G.node[node][probe]
                        aa[n,1] = self.propDict[self.a.G.node[node]['pair']]
                        aa[n,2] = self.a.G.node[node]['pair']
                    else:
                        aa[n,:] = [np.nan, np.nan, np.nan]
                
                aa = aa[~np.isnan(aa)]
                aa.shape = (len(aa)/3,3)
                    
                r,p = stats.pearsonr(aa[:,0], aa[:,1])
                
                if p < self.sigVal:
                    print probe
                    # plot graph
                    out = open(self.outFile, "a")
                    out.writelines(','.join([str(v) for v in [probe, '"'+self.probeDict[probe][1]+'"', r, p]])+'\n')
                    out.close()
                    if plot:
                        plt.scatter(aa[:,1], aa[:,0])
                        plt.savefig(self.outFile.replace('.txt',probe+'.png'), dpi=300)
                        plt.close()
                    
                    # save data
                    datFile = open(self.outFile.replace('.txt', probe+self.gm+'.txt'), "wb")
                    datFile.writelines(' '.join([probe, self.gm, "node", "subj"])+'\n')
                    datFile.writelines('\n'.join([' '.join([str(aa[n,0]), str(aa[n,1]), str(aa[n,2]), self.subj]) for n in range(len(aa[:,1]))]))
                    
        else:
            pass
    
    def probeSubT(self, l, plot, probeNumbers=None):
        '''
        l is a line from the probe file.
        The purpose of this function is to write thresholded data to a datafile
        eg for use in ANOVA
        '''
        probe = l['probe']
        datFile = None
        if probeNumbers:
            if probe in probeNumbers:
                # assign probe values to sample numbers
                for node in self.a.G.nodes():
                    if l[str(node)]:
                        self.a.G.node[node][probe] = l[str(node)]
                    else:
                        self.a.G.node[node][probe] = None
                                      
                    outDict = {probe:probe, 'subj':self.subj}
                    for p in self.propDict.keys():
                        outDict[p] = self.propDict[p]
                    
                    if not datFile:
                        headers = [probe, "subj"]
                        gmSubjs = self.propDict[self.probeDict.keys()[0]].keys()
                        gmSubjs.sort()
                        headers.extend(gmSubjs)

                        datFile = open(self.outFile.replace('.txt', probe+self.gm+'.txt'), "wb")
                        writer = csv.DictWriter(datFile, fieldnames=headers, delimiter=" ")
                        writer.writeheader()
                        
                    writer.writerow(outDict)
        
    def norm(self,x):
        xMin = np.min(x)
        xMax = np.max(x)
        
        x = np.array((x)) - xMin
        x = (x/xMax) * 20
        return([v for v in x])
