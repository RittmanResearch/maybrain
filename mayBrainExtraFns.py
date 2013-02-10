# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 22:25:22 2012

Functions that used to be at the end of networkutils_bin_camb

"""
import os
import numpy as np
import networkx as nx
import csv
from networkx.algorithms import cluster
from networkx.algorithms import centrality

class extraFns():
    
    def globalefficiency(self, G):
        """
        A set of definitions to calculate global efficiency and local efficiency. Definitions are taken from Latora and Marchiori
        2001, Physical Review Letters. 
        """
    
        N = len(G.nodes())      # count nodes
        ssl = 0            # sum of inverse of the shortest path lengths
    
        # nicked from nx.single_source_shortest_path_length to sum shortest lengths
        for node in G.nodes():
            seen={}                  # level (number of hops) when seen in BFS
            level=0                  # the current level
            nextlevel={node:1}  # dict of nodes to check at next level
            while nextlevel:
                thislevel=nextlevel  # advance to next level
                nextlevel={}         # and start a new list (fringe)
                for v in thislevel:
                    if v not in seen: 
                        seen[v]=level # set the level of vertex v
                        nextlevel.update(G[v]) # add neighbors of v
                level=level+1
            if sum(seen.values())>0:
                invpls = [1/float(v) for v in seen.values() if v!=0 ]
                invpl = np.sum(invpls)
                ssl += invpl         # sum inverse shortest path lengths
        
        if N>1:
            Geff = (1/(float(N)*(float(N-1))))*float(ssl)
            return Geff
        else:
            "Number of nodes <1, can't calculate global efficiency"
            return None
        
    
    def localefficiency(self, G):
        nodecalcs = []
        
        for node in G.nodes():
            Ginodes = nx.neighbors(G,node)
            Giedges = G.edges(Ginodes)
            Gi = nx.Graph()
            Gi.add_nodes_from(Ginodes)
            Gi.add_edges_from(Giedges)
            
            Gi_globeff = self.globalefficiency(Gi)
            
            nodecalcs.append(Gi_globeff)
        
        nodecalcs = [float(v) for v in nodecalcs if v]
        
        try:
            loceff = np.sum(nodecalcs) / float(len(G.nodes()))
            return(loceff)
        
        except:
            print "Can not calculate local efficiency"
            print "Local efficiency list for each node:"+','.join(loceff)
            return None
        
    
    def smallworldparameters(self, brain,outfilebase = "brain", append=True):
        """A calculation for average cluster coefficient and average shortest path length (as defined in Humphries
        2008 http://www.plosone.org/article/info:doi/10.1371/journal.pone.0002051).
        """
        outfile = outfilebase+'_smallworldparameters'
        
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            os.rename(outfile,outfile+'.old')
        
        if append and os.path.exists(outfile):
            f = open(outfile,"ab")
        
        else:
            f= open(outfile,"wb")
            f.writelines('ClusterCoeff\tAvShortPathLength\n')
        
        writer = csv.writer(f,delimiter='\t')       
        cc_pl =  (None,None)
        
        # check biggest connected component is define
        if not brain.bigconnG:
            brain.largestConnComp()
    
        # calculate parameters
        try:
            cc_pl = (cluster.average_clustering(brain.bigconnG),nx.average_shortest_path_length(brain.bigconnG))
        
        except:
            if brain.bigconnG.nodes() == []:
                print "No nodes left in network, can not calculate path length, clustering coeff or smallworldness"
            
            else:
                print "No edges, can not calculate shortest path length"
                print "Writing clustering coefficient only"
                cc_pl = (cluster.average_clustering(brain.bigconnG),None)
                
        writer.writerow(cc_pl)
        f.close()
        
        return cc_pl
    
## Everything below this line will soon be deprecated. It's currently kept for backwards
## compatibility.     
    

def degreewrite(brain,outfilebase="brain", append=True):
    # node degrees
    outfile = outfilebase+'_degrees_nodes'

    try:
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            os.rename(outfile,outfile+'.old')
    except:
        pass
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        nodewriter = csv.DictWriter(f,fieldnames = brain.G.nodes())
        
    else:
        f= open(outfile,"wb")
        nodewriter = csv.DictWriter(f,fieldnames = brain.G.nodes())
        headers = dict((n,n) for n in brain.G.nodes())
        nodewriter.writerow(headers)
        
    degs = brain.G.degree()
    degstowrite = dict((n,None) for n in brain.G.nodes())
    for node in degs.keys():
        degstowrite[node] = degs[node]
    nodewriter.writerow(degstowrite)
    f.close()
    
    # hub degrees
    outfile = outfilebase+'_degrees_hubs'
    hubidfile = outfilebase+'_degrees_hubs_ids'

    try:
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            try:
                os.rename(outfile,outfile+'.old')
                os.rename(hubidfile,hubidfile+'.old')
                print "Moving existing degrees file to "+outfile+'.old'
                print "Moving existing hub ID file to "+hubidfile+'.old'
                
            except IOError, error:
                (errorno, errordetails) = error
                print "Error moving files "+errordetails
                print "Carrying on anyway"
    except:
        pass
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        g = open(hubidfile,"ab")
        
    else:
        f= open(outfile,"wb")
        g = open(hubidfile,"wb")
    
    deghubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component

    # write hub identifies to file    
    idwriter = csv.DictWriter(f,fieldnames = brain.hubs)
    hubwriter = csv.DictWriter(g,fieldnames = brain.hubs)
    
    headers = dict((n,n) for n in brain.hubs)
    hubwriter.writerow(headers)
    
    degstowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with degree data

    try:
        degs = brain.G.degree(deghubs)
        for node in degs.keys():
            degstowrite[node] = degs[node]
    except:
        print "no hubs in largest connected component"
    idwriter.writerow(degstowrite)
    f.close()
    g.close()
    
    # write degree histogram
    outfileHist = outfilebase+'_degreehistogram'
    if not append and os.path.exists(outfileHist):
        print "Moving existing file to "+outfileHist+'.old'
        os.rename(outfileHist,outfileHist+'.old')
    
    if append and os.path.exists(outfileHist):
        f = open(outfileHist,"ab")
    
    else:
        f= open(outfileHist,"wb")
        
    degreeHist = nx.degree_histogram(brain.G)
    
    histwriter = csv.writer(f)
    histwriter.writerow(degreeHist)
    f.close()
        

            

def betweennesscentralitywrite(brain,outfilebase = "brain", append=True):
    """ Calculates node and hub betwenness centralities. For hub centralities there are two files, one with the values in and another
    with the hub identities in corresponding rows.
    """
    
    ## betweenness centrality
    # node centrality
    outfile = outfilebase+'_betweenness_centralities_nodes'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f= open(outfile,"ab")
        writer = csv.DictWriter(f,fieldnames = brain.G.nodes())
        
    else:
        f = open(outfile,"wb")
        writer = csv.DictWriter(f,fieldnames = brain.G.nodes())
        headers = dict((n,n) for n in brain.G.nodes())
        writer.writerow(headers)
        
    centralities = centrality.betweenness_centrality(brain.G)  # calculate centralities for largest connected component
    nodecentralitiestowrite = dict((n,None) for n in brain.G.nodes())   # create a blank dictionary of all nodes in the graph
    for node in centralities:
        nodecentralitiestowrite[node] = centralities[node]    # populate the blank dictionary with centrality values
    writer.writerow(nodecentralitiestowrite)                    # write out centrality values
    f.close()
    
    # hub centrality
    outfile = outfilebase+'_betweenness_centralities_hubs'
    hubidfile = outfilebase+'_betweenness_centralities_hubs_ids'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        try:
            os.rename(outfile,outfile+'.old')
            os.rename(hubidfile,hubidfile+'.old')
            print "Moving existing degrees file to "+outfile+'.old'
            print "Moving existing hub ID file to "+hubidfile+'.old'
            
        except IOError, error:
            (errorno, errordetails) = error
            print "Error moving files "+errordetails
            print "Carrying on anyway"
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        g = open(hubidfile,"ab")
        
    else:
        f= open(outfile,"wb")
        g = open(hubidfile,"wb")
        
    centhubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component

    # write hub identifies to file    
    writer = csv.DictWriter(f,fieldnames = brain.hubs)
    hubwriter = csv.DictWriter(g,fieldnames = brain.hubs)
    
    headers = dict((n,n) for n in brain.hubs)         # dictionary of all hubs in network to write
    hubwriter.writerow(headers)
    
    hubcentralitieistowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with centralities data

    for hub in centhubs:
        hubcentralitieistowrite[hub] = nodecentralitiestowrite[hub]
        
    writer.writerow(hubcentralitieistowrite)
    f.close()
    g.close()
    
def closenesscentralitywrite(brain,outfilebase = "brain", append=True):
    """ Calculates node and hub betwenness centralities. For hub centralities there are two files, one with the values in and another
    with the hub identities in corresponding rows.
    """
    ## closeness centrality
    # node centrality
    outfile = outfilebase+'_closeness_centralities_nodes'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f= open(outfile,"ab")
        writer = csv.DictWriter(f,fieldnames = brain.G.nodes())
        
    else:
        f = open(outfile,"wb")
        writer = csv.DictWriter(f,fieldnames = brain.G.nodes())
        headers = dict((n,n) for n in brain.G.nodes())
        writer.writerow(headers)
        
        
    centralities = centrality.closeness_centrality(brain.G)  # calculate centralities for largest connected component
    nodecentralitiestowrite = dict((n,None) for n in brain.G.nodes())   # create a blank dictionary of all nodes in the graph
    for node in centralities:
        nodecentralitiestowrite[node] = centralities[node]    # populate the blank dictionary with centrality values
    writer.writerow(nodecentralitiestowrite)                    # write out centrality values
    f.close()
    
    # hub centrality
    outfile = outfilebase+'_closeness_centralities_hubs'
    hubidfile = outfilebase+'_closeness_centralities_hubs_ids'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        try:
            os.rename(outfile,outfile+'.old')
            os.rename(hubidfile,hubidfile+'.old')
            print "Moving existing degrees file to "+outfile+'.old'
            print "Moving existing hub ID file to "+hubidfile+'.old'
            
        except IOError, error:
            (errorno, errordetails) = error
            print "Error moving files "+errordetails
            print "Carrying on anyway"
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        g = open(hubidfile,"ab")
        
    else:
        f= open(outfile,"wb")
        g = open(hubidfile,"wb")
        
    centhubs = [hub for hub in brain.hubs if hub in brain.G] # hubs within largest connected graph component

    # write hub identifies to file    
    writer = csv.DictWriter(f,fieldnames = brain.hubs)
    hubwriter = csv.DictWriter(g,fieldnames = brain.hubs)
    
    headers = dict((n,n) for n in brain.hubs)         # dictionary of all hubs in network to write
    hubwriter.writerow(headers)
    
    hubcentralitieistowrite = dict((n,None) for n in brain.hubs) # empty dictionary to populate with centralities data

    for hub in centhubs:
        hubcentralitieistowrite[hub] = nodecentralitiestowrite[hub]
        
    writer.writerow(hubcentralitieistowrite)
    f.close()
    g.close()

def efficiencywrite(brain,outfilebase = "brain", append=True):
    """
    Writes to a file the output of global and local efficiencies for the largest connected component of a network.
    """
    from maybrain import analysis
    outfile = outfilebase+'_efficiency'
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")
    
    # set headers
    effs = {"Globalefficiency":None,"Localefficiency":None}
            
    # local efficiency
    effs["Localefficiency"] = analysis.localefficiency(brain.G)
        
    # global efficiency
    effs["Globalefficiency"] = analysis.globalefficiency(brain.G)
    
    # write results to file
    writer = csv.DictWriter(f,fieldnames = effs.keys())
    if brain.iter == None:
        f.writelines("Localefficiency,Globalefficiency\n")
#        writer.writeheader()
    writer.writerow(effs)
    f.close()

def modularstructure(brain,outfilebase = "brain", redefine_clusters=True, append=True):
    """
    Writes to a file the output of cluster membership and modularity.
    """
    # redefine graph clusters
    if redefine_clusters == True:
        brain.clusters()
    
    # write out modularity
    outfile = outfilebase+'_modularity'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
    
    else:
        f= open(outfile,"wb")

    writer = csv.writer(f,delimiter='\t')
    writer.writerow([brain.modularity])
    f.close()

def writeEdgeLengths(brain,outfilebase = "brain", append=True):
    outfile = outfilebase+'_edgeLengths'
    try:
        if not append and os.path.exists(outfile):
            print "Moving existing file to "+outfile+'.old'
            os.rename(outfile,outfile+'.old')
        
        if append and os.path.exists(outfile):
            f = open(outfile,"ab")
        
        else:
            f= open(outfile,"wb")
            
        writer = csv.writer(f,delimiter='\t')
        writer.writerow([brain.lengthEdgesRemoved])
        f.close()
        
    except AttributeError:
        pass
        
    # mean edge lengths
    edgeLengths = []
    for edge in brain.G.edges():
        edgeLengths.append(np.linalg.norm(np.array(brain.G.node[edge[0]]['xyz']) - np.array(brain.G.node[edge[1]]['xyz'])))
        
    meanEdgeLengths = np.mean(np.absolute(edgeLengths))
    
    hubEdgeLengths = []
    for edge in brain.G.edges(brain.hubs):
        hubEdgeLengths.append(np.linalg.norm(np.array(brain.G.node[edge[0]]['xyz']) - np.array(brain.G.node[edge[1]]['xyz'])))

    meanHubLengths = np.mean(np.absolute(hubEdgeLengths))
    
    outfile = outfilebase+'_meanEdgeLengths'
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
            
    else:
        f= open(outfile,"wb")
        f.writelines("MeanEdgeLengths\tMeanHubEdgeLengths\n")
        
    writer = csv.writer(f,delimiter='\t')
    writer.writerow([meanEdgeLengths, meanHubLengths])
    f.close()
    
def writeEdgeNumber(brain, outfilebase = "brain", append=True):
    outfile = outfilebase+'_EdgeNumber'
    
    
    if not append and os.path.exists(outfile):
        print "Moving existing file to "+outfile+'.old'
        os.rename(outfile,outfile+'.old')
    
    if append and os.path.exists(outfile):
        f = open(outfile,"ab")
        
    else:
        f= open(outfile,"wb")
        f.writelines("EdgeNumber\n")
        
    edgeNum = len(brain.G.edges())
    writer = csv.writer(f,delimiter='\t')
    writer.writerow([edgeNum])
    f.close()
    
