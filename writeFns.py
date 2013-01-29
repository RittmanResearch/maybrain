# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 10:23:18 2013

@author: martyn

Functions to write output of a brain object. Now object oriented.

"""

from os import path, rename
from csv import writer, DictWriter
from networkx.algorithms import centrality
from networkx import degree_histogram
from numpy import array, sqrt, mean, sum, shape

class writeFns():
    ''' a set of functions to write various things from a maybrain brainObj to file '''

    def fileCheck(self, fname, boolVal):
        ''' perform a check to see if old file is there and move if necessary '''
        if boolVal:
            return path.exists(fname)
        
        if path.exists(fname):
            print 'moving ' + fname + ' to ' + fname + '.old'
            rename(fname, fname + '.old')
            return True
        return False
        
    def degreewrite(self, brain, outfilebase="brain", append=True):
        # node degrees
        outfile = outfilebase+'_degrees_nodes'    
        boolVal = self.fileCheck(outfile, append)
            
        # open file to write to
        if append and boolVal:
            # case where file is appended
            f = open(outfile,"ab")
            nodewriter = DictWriter(f,fieldnames = brain.G.nodes())
            
        else:
            # case where a new file is created
            f= open(outfile,"wb")
            nodewriter = DictWriter(f,fieldnames = brain.G.nodes())
            # write header based on node names
            headers = dict((n,n) for n in brain.G.nodes())
            nodewriter.writerow(headers)
            
        # get degrees of each node
        degs = brain.G.degree()
        degstowrite = dict((n,None) for n in brain.G.nodes())
        for node in degs.keys():
            degstowrite[node] = degs[node]
        # write out
        nodewriter.writerow(degstowrite)
        f.close()
        
        ## ===========================================================================    
        
        ## hub degrees
        outfile = outfilebase+'_degrees_hubs'
        hubidfile = outfilebase+'_degrees_hubs_ids'

        # move old file if necessary
        boolOF = self.fileCheck(outfile, append)
        self.fileCheck(hubidfile, append)
    
        # open files to write to
        # bit dodgy, should check both files before "ab" opening...
        if append and boolOF:
            f = open(outfile,"ab")
            g = open(hubidfile,"ab")
            
        else:
            f= open(outfile,"wb")
            g = open(hubidfile,"wb")
        
        # hubs within largest connected graph component
        deghubs = [hub for hub in brain.hubs if hub in brain.G] 
    
        # write hub identifies to file    
        idwriter = DictWriter(f,fieldnames = brain.hubs)
        hubwriter = DictWriter(g,fieldnames = brain.hubs)
        
        headers = dict((n,n) for n in brain.hubs)
        hubwriter.writerow(headers)
        
        # empty dictionary to populate with degree data
        degstowrite = dict((n,None) for n in brain.hubs) 
    
        # get degrees of hubs
        try:
            degs = brain.G.degree(deghubs)
            for node in degs.keys():
                degstowrite[node] = degs[node]
        except:
            # why not check for this specifically beforehand than using try...except?
            print "no hubs in largest connected component"
        
        # write to file
        idwriter.writerow(degstowrite)
        f.close()
        g.close()
        
        ## =======================================================================    
        
        ## write degree histogram
        outfileHist = outfilebase+'_degreehistogram'
        histBool = self.fileCheck(outfileHist, append)    
    
        # open file    
        if append and histBool:
            f = open(outfileHist,"ab")    
        else:
            f= open(outfileHist,"wb")
            
        # get histogram
        degreeHist = degree_histogram(brain.G)
        
        # write histogram to file
        histwriter = writer(f)
        histwriter.writerow(degreeHist)
        f.close()
        
     
    
                
    
    def betweennesscentralitywrite(self, brain,outfilebase = "brain", append=True):
        """ Calculates node and hub betweenness centralities. For hub centralities there are two files, one with the values in and another
        with the hub identities in corresponding rows.
        """
        
        ## betweenness centrality
        # node centrality
        outfile = outfilebase+'_betweenness_centralities_nodes'
        boolVal = self.fileCheck(outfile, append)
        
        if append and boolVal:
            f= open(outfile,"ab")
            writeObj = DictWriter(f,fieldnames = brain.G.nodes())
            
        else:
            f = open(outfile,"wb")
            writeObj = DictWriter(f,fieldnames = brain.G.nodes())
            headers = dict((n,n) for n in brain.G.nodes())
            writeObj.writerow(headers)
            
        centralities = centrality.betweenness_centrality(brain.G)  # calculate centralities for largest connected component
        nodecentralitiestowrite = dict((n,None) for n in brain.G.nodes())   # create a blank dictionary of all nodes in the graph
        for node in centralities:
            nodecentralitiestowrite[node] = centralities[node]    # populate the blank dictionary with centrality values
        writeObj.writerow(nodecentralitiestowrite)                    # write out centrality values
        f.close()
        
        ## ==================================================================        
        
        ## hub centrality
        outfile = outfilebase+'_betweenness_centralities_hubs'
        hubidfile = outfilebase+'_betweenness_centralities_hubs_ids'
        
        OFbool = self.fileCheck(outfile, append)
        self.fileCheck(hubidfile, append)
        
        if append and OFbool:
            f = open(outfile,"ab")
            g = open(hubidfile,"ab")
            
        else:
            f= open(outfile,"wb")
            g = open(hubidfile,"wb")
            
        centhubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component
    
        # write hub identifies to file    
        writeObj = DictWriter(f,fieldnames = brain.hubs)
        hubwriter = DictWriter(g,fieldnames = brain.hubs)
        
        headers = dict((n,n) for n in brain.hubs)         # dictionary of all hubs in network to write
        hubwriter.writerow(headers)
        
        hubcentralitieistowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with centralities data
    
        for hub in centhubs:
            hubcentralitieistowrite[hub] = nodecentralitiestowrite[hub]
            
        writeObj.writerow(hubcentralitieistowrite)
        f.close()
        g.close()
        
        
    def closenesscentralitywrite(self, brain,outfilebase = "brain", append=True):
        """ Calculates node and hub closeness centralities. For hub centralities there are two files, one with the values in and another
        with the hub identities in corresponding rows.
        """
        ## closeness centrality
        # node centrality
        outfile = outfilebase+'_closeness_centralities_nodes'
        boolVal = self.fileCheck(outfile, append)        
        
        # open file and write headers if necessary
        if append and boolVal:
            f= open(outfile,"ab")
            headers = None            
        else:
            f = open(outfile,"wb")
            headers = dict((n,n) for n in brain.G.nodes())
                        
        writeObj = DictWriter(f,fieldnames = brain.G.nodes())            
        if headers:
            writeObj.writerow(headers)

        # make calculations            
        centralities = centrality.closeness_centrality(brain.G)  # calculate centralities for largest connected component
        nodecentralitiestowrite = dict((n,None) for n in brain.G.nodes())   # create a blank dictionary of all nodes in the graph
        for node in centralities:
            nodecentralitiestowrite[node] = centralities[node]    # populate the blank dictionary with centrality values
        writeObj.writerow(nodecentralitiestowrite)                    # write out centrality values
        f.close()
        
        # hub centrality
        outfile = outfilebase+'_closeness_centralities_hubs'
        hubidfile = outfilebase+'_closeness_centralities_hubs_ids'
        
        OFbool = self.fileCheck(outfile, append)
        self.fileCheck(hubidfile, append)        
        
        if append and OFbool:
            f = open(outfile,"ab")
            g = open(hubidfile,"ab")
            
        else:
            f= open(outfile,"wb")
            g = open(hubidfile,"wb")
            
        centhubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component
    
        # write hub identifies to file    
        writeObj = DictWriter(f,fieldnames = brain.hubs)
        hubwriter = DictWriter(g,fieldnames = brain.hubs)
        
        headers = dict((n,n) for n in brain.hubs)         # dictionary of all hubs in network to write
        hubwriter.writerow(headers)
        
        hubcentralitieistowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with centralities data
    
        for hub in centhubs:
            hubcentralitieistowrite[hub] = nodecentralitiestowrite[hub]
            
        writeObj.writerow(hubcentralitieistowrite)
        f.close()
        g.close()
    
    def efficiencywrite(self, brain, localEffs, globalEffs, outfilebase = "brain", append=True):
        """
        Writes to a file the output of global and local efficiencies for the 
        largest connected component of a network.
        """
        outfile = outfilebase+'_efficiency'
        boolVal = self.fileCheck(outfile, append)
        
        if append and boolVal:
            f = open(outfile,"ab")        
        else:
            f= open(outfile,"wb")
        
        # set data
        effs = {"Globalefficiency":globalEffs ,"Localefficiency":localEffs}
        
        # write results to file
        writeObj = DictWriter(f,fieldnames = effs.keys())
        if brain.iter == None:
            f.writelines("Localefficiency,Globalefficiency\n")
    #        writer.writeheader()
        writeObj.writerow(effs)
        f.close()
    
    def modularstructure(self, brain,outfilebase = "brain", redefine_clusters=True, append=True):
        """
        Writes to a file the output of cluster membership and modularity.
        """
        # redefine graph clusters
        if redefine_clusters == True:
            brain.clusters()
        
        # write out modularity
        outfile = outfilebase+'_modularity'
        boolVal = self.fileCheck(outfile, append)
        
        if append and boolVal:
            f = open(outfile,"ab")        
        else:
            f= open(outfile,"wb")
    
        # write to file
        writeObj = writer(f,delimiter='\t')
        writeObj.writerow([brain.modularity])
        f.close()
    
    def writeEdgeLengths(self, brain,outfilebase = "brain", append=True):
        ''' change log: 29/1/2013, removed unnecessary looping when calculating avg edge lengths'''
        outfile = outfilebase+'_edgeLengths'
        boolVal = self.fileCheck(outfile, append)
        # remove this try...except - it's just laziness! Check the attribute you want to test for first rather than relying on an error
           
        if append and boolVal:
            f = open(outfile,"ab")            
        else:
            f= open(outfile,"wb")

        # not sure exactly what I should be checking here - double check what happens
        if brain.lengthEdgesRemoved:                
            writeObj = writer(f,delimiter='\t')
            writeObj.writerow(brain.lengthEdgesRemoved)
            f.close()

        # get edges
        n1 = []
        n2 = []
        for edge in brain.G.edges():            
            n1.append(brain.G.node[edge[0]]['xyz'])
            n2.append(brain.G.node[edge[1]]['xyz'])
        n1 = array(n1)
        n2 = array(n2)
        print("in writeEdgeLengths")
        print(shape(n1), shape(n2))
        
        # calculate mean edge length
        l = sqrt(sum((n1-n2)**2, axis=1))
        meanEdgeLengths = mean(l)        
        
        # do the same for hubs
        hn1 = []
        hn2 = []
        for edge in brain.G.edges(brain.hubs):
            hn1.append(brain.G.node[edge[0]]['xyz'])
            hn2.append(brain.G.node[edge[1]]['xyz'])
        hn1 = array(hn1)
        hn2 = array(hn2)
        
        l = sqrt(sum((hn1-hn2)**2, axis=1))
        meanHubLengths = mean(l)            
            
        # open output file
        outfile = outfilebase+'_meanEdgeLengths'
        self.fileCheck(outfile)        
        
        if append and path.exists(outfile):
            f = open(outfile,"ab")
        else:
            f= open(outfile,"wb")
            f.writelines("MeanEdgeLengths\tMeanHubEdgeLengths\n")
        
        # write output to file
        writeObj = writer(f,delimiter='\t')
        writeObj.writerow([meanEdgeLengths, meanHubLengths])
        f.close()
        
    def writeEdgeNumber(self, brain, outfilebase = "brain", append=True):
        ''' Write number of edges to file'''
        outfile = outfilebase+'_EdgeNumber'
        self.fileCheck(outfile)
        
        if append and path.exists(outfile):
            f = open(outfile,"ab")            
        else:
            f= open(outfile,"wb")
            f.writelines("EdgeNumber\n")
            
        edgeNum = len(brain.G.edges())
        writeObj = writer(f,delimiter='\t')
        writeObj.writerow([edgeNum])
        f.close()
    
