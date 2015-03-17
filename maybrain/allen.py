# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 19:05:26 2014
Module for linking with the Allen brain atlas
@author: tim
"""

import maybrain.brainObjs as mbo

import csv
from os import path,rename
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from glob import glob

class allenBrain:
    def __init__(self, allenSubj, assocMat, delim=",",
                 spatialFile="parcel_500.txt", nodesToExclude=None):
       
        self.subj = allenSubj
        self.fName = "SampleAnnot.csv"
        self.maFile = "MicroarrayExpression.csv"
       
        # set up brain for expression data
        self.a = mbo.brainObj()
       
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
        self.c = mbo.brainObj()
        self.c.importAdjFile(assocMat, delimiter=delim,
                             exclnodes=nodesToExclude)
        self.c.importSpatialInfo(spatialFile)
   
    def comparison(self):
        # set up dictionary to link nodes from probe data and graph
        nodeDictMRIs = {}
       
        for node in self.c.G.nodes():
            dOther = (None, 999.) # dummy length of 999
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
        p = mbo.plotObj()
        p.plotSkull(self.a, contourVals = [3000,9000])
        p.plotBrainCoords(self.c, nodes =  self.c.G.nodes(),
                                     col=(1,0,0), sizeList=5)
        p.plotBrainCoord(self.a, nodes = self.a.G.nodes(),
                                     col=(0,0,1), sizeList=5)
        self.saveFig()
                                     
    def probeData(self, propDict, graphMetric="gm", nodeList=None, plot=False,
                  probeList=[], probeNumbers=[], sigVal=1.0, T=False):
        '''
       
        '''
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

class multiSubj:
    def __init__(self, assocMat, nodesToExclude=[], delim=" ",
                 subjList=None, spatialFile="parcel_500.txt"):
        if subjList:
            self.subjList = subjList
        else:
            self.subjList = [v for v in glob("17823*") if path.isdir(v)]
           
        self.fName = "SampleAnnot.csv"
        self.maFile = "MicroarrayExpression.csv"
        self.probeFile = "Probes.csv"
       
        # set up brain for expression data
        self.a = mbo.brainObj()
       
        self.headers={}
        for subj in self.subjList:
            # import probe data
            f = open(path.join(subj, self.fName), "rb")        
            reader = csv.DictReader(f, delimiter=",", quotechar='"')
           
            self.headers[subj] = ['probe']
            for l in reader:
                n = int(l["structure_id"])
                if not self.a.G.has_node(n):
                    self.a.G.add_node(n)
                    self.a.G.node[n] = l
                self.headers[subj].append(l["structure_id"])
            f.close()
           
            # convert location data for Allen brain from MNI space
            for n in self.a.G.nodes():
                x = 45 - (float(self.a.G.node[n]['mni_x'])/2)
                y = 63 + (float(self.a.G.node[n]['mni_y'])/2)
                z = 36 + (float(self.a.G.node[n]['mni_z'])/2)
                self.a.G.node[n]['xyz'] = (x,y,z)
       
        # set up brain with graph properties
        self.c = mbo.brainObj()
        self.c.importAdjFile(assocMat, delimiter=delim, exclnodes=nodesToExclude)
        self.c.importSpatialInfo(spatialFile)
   
    def comparison(self):
        # set up dictionary to link nodes from probe data and graph
        nodeDictMRIs = {}
       
        for node in self.c.G.nodes():
            dOther = (None, 999.) # dummy length of 999
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

    def probeData(self, probeNumbers=[]):
        for subj in self.subjList:
            # import probe data
            f = open(path.join(subj, self.maFile), "rb")
            reader = csv.DictReader(f, delimiter=",",
                                    fieldnames=self.headers[subj],
                                    quotechar='"')
            for l in reader:
                probe = l['probe']
                if probe in probeNumbers:
                    # assign probe values to sample numbers
                    for node in self.a.G.nodes():
                        if not probe in self.a.G.node[node].keys():
                            self.a.G.node[node][probe] = {}
                           
                        if str(node) in l.keys():
                            self.a.G.node[node][probe][subj] = l[str(node)]
            f.close()
            del(reader)
                                   
        for n in self.a.G.nodes():
            for probe in probeNumbers:
                self.a.G.node[n][probe] = np.mean([float(v) for v in self.a.G.node[n][probe].values()])
               
    def writeXMatrix(self, outFile="Xmatrix.csv", probeNumbers=None):
        # get all probes if otherwise unspecified
        if not probeNumbers:
            f = open(path.join(self.subjList[0], self.probeFile))
            reader = csv.DictReader(f, delimiter=",", quotechar='"')
            probeNumbers = [l['probe_id'] for l in reader]
            del(reader)
            f.close()  
   
        # set up out file
        out = open(outFile, "wb")
        headers = ["Gene"]
        headers.extend([str(v) for v in self.c.G.nodes()])
        writer = csv.DictWriter(out, fieldnames = headers, delimiter=" ")
        writer.writeheader()
       
        # set up matrix
        x = len(self.subjList)
        y = len(probeNumbers)
        z = len(self.c.G.nodes())
        print(x,y,z)
        probeMat = np.memmap("tempMat.txt",
                             dtype="float64",
                             mode="w+",
                             shape=(x,y,z))
                             
        # get the corresponding node names in the MRI graph
        cNodes = {str(self.a.G.node[v]['pair']):v for v in self.a.G.nodes()}
       
        # set up gene list
        geneFlag = True
        pFile = open(path.join(self.subjList[0], self.probeFile))
        pReader = csv.DictReader(pFile, delimiter=",", quotechar='"')
        pDict = {l['probe_id']:l['gene_symbol'] for l in pReader}
        geneList = pDict.values()
        set(geneList)
        geneList = {gene:[] for gene in geneList}
       
        # assign values to matrix
        for x,subj in enumerate(self.subjList):
            # import probe data
            f = open(path.join(subj, self.maFile), "rb")
            reader = csv.DictReader(f, delimiter=",",
                                    fieldnames=self.headers[subj],
                                    quotechar='"')
            y = 0
            for l in reader:
                probe = l['probe']
                if probe in probeNumbers:
                    # assign probe values to sample numbers
                    for z,cNode in enumerate(self.c.G.nodes()):
                        aNode = cNodes[str(cNode)]
                        if str(aNode) in l.keys():
                            probeMat[x,y,z] = float(l[str(aNode)])
                    if geneFlag:
                        geneList[pDict[probe]].append(y)
                    y+=1
            f.close()
            del(reader)

            # normalise expression levels for each probe within subject
            for y in range(probeMat.shape[1]):
                # create a masked array removing the 0. values
                subjMat = np.ma.array(probeMat[x,y,:],
                                      mask=probeMat[x,y,:]==0.,
                                      dtype="float64")

                if len(subjMat[subjMat.mask]>1):
                    subjMat = (subjMat - np.mean(subjMat)) / np.std(subjMat)
                else:
                    subjMat = subjMat - np.mean(subjMat)
                probeMat[x,y,:] = subjMat
            geneFlag=False
       
        # collapse across subjects and probes by gene
        geneNames = geneList.keys()
        geneNames.sort() # sort in to alphabetical order
       
        for g,gene in enumerate(geneNames):
            if (geneList[gene]):
                x = probeMat.shape[2] # number of nodes
                y = probeMat.shape[0]*len(geneList[gene]) # subjects * probes
                               
                geneMat = np.zeros(shape=(x,y), dtype="float64")
               
                for n,p in enumerate(geneList[gene]):
                    for s,subj in enumerate(self.subjList):
                        geneMat[:,n*len(self.subjList)+s] = probeMat[s,p,:]
                   
                geneMat = np.ma.array(geneMat, mask=geneMat==0.)
               
                self.geneMat = geneMat
                meanGene = np.mean(geneMat, axis=1)
               
                outDict = dict(zip([str(v) for v in self.c.G.nodes()], ["{:10.20f}".format(v) for v in meanGene]))
                outDict["Gene"] = gene
                writer.writerow(outDict)
       
        self.geneList = geneList
        self.probeMat = probeMat
        out.close()
   
    def writeYMatrixGroup(self, metricDict, subj="Control", outFile="YmatrixGroup.csv"):
        '''
        Collates metrics in to a matrix for use in partial least squares analysis.
        Note, the metricDict contains the metric name as a key and filename as
        the value. Takes group level measures
        '''
        out = open(outFile, "wb")
        headers = ["Metric"]
        headers.extend([str(v) for v in self.c.G.nodes()])
        writer = csv.DictWriter(out, fieldnames=headers)
        writer.writeheader()
       
        # iterate through the metrics
        for m in metricDict.keys():
            f = open(path.join(subj, metricDict[m]), "r")
            reader = csv.DictReader(f, delimiter=" ")
            l = reader.next()
           
            # remove non-numeric keys
            for v in l.keys():
                try:
                    int(v)
                except:
                    del(l[v])

            mDict = {v:l[v] for v in l.keys() if int(v) in self.c.G.nodes()}
            f.close()
           
            # normalise within metric
            meanMetric = np.mean([float(v) for v in mDict.values()])
            sdMetric = np.std([float(v) for v in mDict.values()])

            mDict = {str(v):(float(mDict[str(v)])-meanMetric)/sdMetric for v in mDict.keys()}
           
            mDict["Metric"] = m
           
            writer.writerow(mDict)
           
    def writeYMatrixIndividuals(self, metricDict, subjList, outFile="YmatrixInd.csv"):
        '''
        Collates metrics in to a matrix for use in partial least squares analysis.
        Note, the metricDict contains the metric name as a key and filename as
        the value. Takes metrics for individual subjects defined in the subject list.
        '''
        out = open(outFile, "wb")
        headers = ["Metric", "Subject"]
        headers.extend([str(v) for v in self.c.G.nodes()])
        writer = csv.DictWriter(out, fieldnames=headers)
        writer.writeheader()
       
        # iterate through the metrics
        for m in metricDict.keys():
            for subj in subjList:            
                f = open(path.join(subj, metricDict[m]), "r")
                reader = csv.DictReader(f, delimiter=" ")
                l = reader.next()
               
                # remove non-numeric keys
                for v in l.keys():
                    try:
                        int(v)
                    except:
                        del(l[v])
   
                mDict = {v:l[v] for v in l.keys() if int(v) in self.c.G.nodes()}
                f.close()
               
                # normalise within metric
                meanMetric = np.mean([float(v) for v in mDict.values()])
                sdMetric = np.std([float(v) for v in mDict.values()])
   
                mDict = {str(v):(float(mDict[str(v)])-meanMetric)/sdMetric for v in mDict.keys()}
               
                mDict["Metric"] = m
                mDict["Subject"] = subj
               
                writer.writerow(mDict)
