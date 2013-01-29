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



def globalefficiency(G):
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
    

def localefficiency(G):
    nodecalcs = []
    
    for node in G.nodes():
        Ginodes = nx.neighbors(G,node)
        Giedges = G.edges(Ginodes)
        Gi = nx.Graph()
        Gi.add_nodes_from(Ginodes)
        Gi.add_edges_from(Giedges)
        
        Gi_globeff = globalefficiency(Gi)
        
        nodecalcs.append(Gi_globeff)
    
    nodecalcs = [float(v) for v in nodecalcs if v]
    
    try:
        loceff = np.sum(nodecalcs) / float(len(G.nodes()))
        return(loceff)
    
    except:
        print "Can not calculate local efficiency"
        print "Local efficiency list for each node:"+','.join(loceff)
        return None
    

def smallworldparameters(brain,outfilebase = "brain", append=True):
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
    effs["Localefficiency"] = localefficiency(brain.G)
        
    # global efficiency
    effs["Globalefficiency"] = globalefficiency(brain.G)
    
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
        writer.writerow(brain.lengthEdgesRemoved)
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
    
    
def writeout(brain,outfilebase="brain"):
    if not brain.iter:
        outfile = '_'.join(outfilebase,'0.txt')
    else:
        outfile = '_'.join(outfilebase, str(brain.iter)+'.txt')
    
    
#def mayavivis(brain,outfile,nbunch = None):
#    """
#    draws 3D version of brain
#    """
#    if nbunch:
#        # set up arrays of data for 3D visualisation
#        visarray = np.array([[0.0]*len(nbunch)]*4)
#        for n in range(len(nbunch)):
#            node = nbunch[n]
#            for x in range(3):
#                visarray[x:x+1,n:n+1] = float(brain.G.node[node]['xyz'][x])/200
#            visarray[3:4,n:n+1] = brain.G.node[node]['cluster']
#        
#        hubs = [v for v in brain.hubs if v in nbunch]
#        
#        visarrayhubs = np.array([[0.0]*len(hubs)]*4)
#        for n in range(len(hubs)):
#            hub = hubs[n]
#            for x in range(3):
#                visarrayhubs[x:x+1,n:n+1] = float(brain.G.node[hub]['xyz'][x])/200
#            visarrayhubs[3:4,n:n+1] = brain.G.node[hub]['cluster']
#            brain.clustervisarrayhubs = visarrayhubs
#            
#        degennodes = [v for v in nbunch if brain.G.node[v]['degenerating']]
#        visarraydegens = np.array([[0.0]*len(degennodes)]*4)
#        for n in range(len(degennodes)):
#            degennode = degennodes[n]
#            for x in range(3):
#                visarraydegennodes[x:x+1,n:n+1] = float(brain.G.node[degennode]['xyz'][x])/200
#            visarraydegennodes[3:4,n:n+1] = brain.G.node[degennode]['cluster']
#            brain.clustervisarraydegennodes = visarraydegennodes
#        edges = [v for v in brain.G.edges(nbunch) if v[0] in nbunch and v[1] in nbunch]
#
#    else:    
##    set up hub 3D locations
#        visarray = np.array([[0.0]*len(brain.G.nodes())]*4)
#        for n in range(len(brain.G.nodes())):
#            node = brain.G.nodes()[n]
#            for x in range(3):
#                visarray[x:x+1,n:n+1] = float(brain.G.node[node]['xyz'][x])/200
#            visarray[3:4,n:n+1] = brain.G.node[node]['cluster']
#
#        visarrayhubs = np.array([[0.0]*len(brain.hubs)]*4)
#        for n in range(len(brain.hubs)):
#            hub = brain.hubs[n]
#            for x in range(3):
#                visarrayhubs[x:x+1,n:n+1] = float(brain.G.node[hub]['xyz'][x])/200
#            visarrayhubs[3:4,n:n+1] = brain.G.node[hub]['cluster']
#            
#        
#        try:
#            degennodes = [v for v in brain.G.nodes() if brain.G.node[v]['degenerating']]
#            visarraydegennodes = np.array([[0.0]*len(degennodes)]*4)
#            for n in range(len(degennodes)):
#                degennode = degennodes[n]
#                for x in range(3):
#                    visarraydegennodes[x:x+1,n:n+1] = float(brain.G.node[degennode]['xyz'][x])/200
#                visarraydegennodes[3:4,n:n+1] = brain.G.node[degennode]['cluster']
#        except:
#            pass
#        
#        
#        edges = brain.G.edges()
#        
#    
#    mlab.figure(1,bgcolor=(1,1,1),size=(1200,900))
#    
#    pts = mlab.points3d(visarray[0],visarray[1],visarray[2],visarray[3],scale_factor=0.01,scale_mode='none',colormap='Greens',opacity=0)
#    pts.mlab_source.dataset.lines = np.array(edges)
#    
#    ptshubs = mlab.points3d(visarrayhubs[0],visarrayhubs[1],visarrayhubs[2],scale_factor=0.01,color=(1,0,0))
#    try:
#        ptsdegens = mlab.points3d(visarraydegennodes[0],visarraydegennodes[1],visarraydegennodes[2],scale_factor=0.02,color=(1,1,0))
#    except:
#        pass
#    
#    tube = mlab.pipeline.tube(pts,tube_radius=0.0002)
#    mlab.pipeline.surface(tube,colormap='Accent',opacity=1)
#    brain.picturecount = 1
#    while os.path.exists(outfile+'_3D_'+str(brain.picturecount)+'.png'):
#        brain.picturecount += 1
#    mlab.savefig(outfile+'_3D_'+str(brain.picturecount)+'.png')
#    scene.scene.save_x3d(outfile+".x3d")
#    mlab.close(all=True)
#    mlab.show()
#    
#def flatpictures(brain,outfile='Brain',nbunch=None):
#    """
#    draws 2D version of brain, or subset of nodes
#    """
#    if nbunch:
#        hubs = [v for v in brain.hubs if v in nbunch]
#        try:
#            degennodes = [v for v in nbunch if brain.G.node[v]['degenerating']]
#        except:
#            pass
#        nodes = [v for v in nbunch if v not in hubs]
#        
#        crossclusteredges = [v for v in brain.G.edges(nbunch) if v[0] in nbunch and v[1] in nbunch]
#        noncrossingedges = []
#            
#    else:
#        hubs = brain.hubs
#        nodes = [v for v in brain.G.nodes() if v not in brain.hubs]
#        try:
#            degennodes = [v for v in brain.G.nodes() if brain.G.node[v]['degenerating']]
#        except:
#            pass
#        
#        # define edges crossing or not crossing clusters
#        crossclusteredges = [v for v in brain.G.edges() if brain.G.node[v[0]]['cluster'] != brain.G.node[v[1]]['cluster']]
#        noncrossingedges = [v for v in brain.G.edges() if brain.G.node[v[0]]['cluster'] == brain.G.node[v[1]]['cluster']]
#    
#    # draw whole brain/cluster
#    cm = plt.get_cmap('Autumn')
#        
#    plt.figure(1)
#    plt.clf()
#    draw_networkx_edges(brain.G,brain.pos,edgelist=crossclusteredges,edge_color='Black',width=0.2)
#    draw_networkx_nodes(brain.G,brain.pos,nodelist=nodes,node_size=0,)
#    try:
#        draw_networkx_nodes(brain.G,brain.pos,nodelist=degennodes,node_size=20,node_color="blue")
#    except:
#        pass
#    
#    draw_networkx_nodes(brain.G,brain.pos,nodelist=hubs,node_size=10)
#    colors = [brain.G.node[v[0]]['cluster'] for v in noncrossingedges]
#    draw_networkx_edges(brain.G,brain.pos,edgelist=noncrossingedges,edge_color=colors,edge_cmap=cm,width=0.2)
#    plt.axis('off')
#    
#    brain.picturecount = 1
#    while os.path.exists(outfile+'_2D_'+str(brain.picturecount)+'.png'):
#        brain.picturecount += 1
#    
#    plt.savefig(outfile+'_2D_'+str(brain.picturecount)+'.png',dpi=75,figsize=(2,1.5))
#
#    plt.close()

#def glassBrain(brain, bgImage='../average4_thr.mnc',outfile="glassbrain", nbunch=None):
#    try:
#        engine = mayavi.engine
#    except NameError:
#        from mayavi.api import Engine
#        engine = Engine()
#        engine.start()
#    if len(engine.scenes) == 0:
#        engine.new_scene()
#    
#    scene = engine.scenes[0]
#    scene.scene.background = (1.0, 1.0, 1.0)
#    scene.scene.show_axes = True
#           
#    mlab.set_engine(engine)
#    
#    # import background image
#    bg = nibabel.load(bgImage)
#    a = bg.get_data()
#    src = mlab.contour3d(a, colormap='bone', figure=scene)
#    
#    src.actor.property.representation = 'points'    
#
#    
#    # draw network
#    posCorrectionDict = {0:91, 1:109, 2:91}
#    if nbunch:
#        # set up arrays of data for 3D visualisation
#        visarray = np.array([[0.0]*len(nbunch)]*4)
#        for n in range(len(nbunch)):
#            node = nbunch[n]
#            for x in range(3):
#                visarray[x:x+1,n:n+1] = float(brain.G.node[node]['xyz'][x])##*2-posCorrectionDict[x]
#            visarray[3:4,n:n+1] = brain.G.node[node]['cluster']
#        
#        hubs = [v for v in brain.hubs if v in nbunch]
#        
#        visarrayhubs = np.array([[0.0]*len(hubs)]*4)
#        for n in range(len(hubs)):
#            hub = hubs[n]
#            for x in range(3):
#                visarrayhubs[x:x+1,n:n+1] = float(brain.G.node[hub]['xyz'][x])#*2-posCorrectionDict[x]
#            visarrayhubs[3:4,n:n+1] = brain.G.node[hub]['cluster']
#            brain.clustervisarrayhubs = visarrayhubs
#            
#        degennodes = [v for v in nbunch if brain.G.node[v]['degenerating']]
#        visarraydegens = np.array([[0.0]*len(degennodes)]*4)
#        for n in range(len(degennodes)):
#            degennode = degennodes[n]
#            for x in range(3):
#                visarraydegennodes[x:x+1,n:n+1] = float(brain.G.node[degennode]['xyz'][x])#*2-posCorrectionDict[x]
#            visarraydegennodes[3:4,n:n+1] = brain.G.node[degennode]['cluster']
#            brain.clustervisarraydegennodes = visarraydegennodes
#        edges = [v for v in brain.G.edges(nbunch) if v[0] in nbunch and v[1] in nbunch]
#
#    else:    
##    set up hub 3D locations
#        visarray = np.array([[0.0]*len(brain.G.nodes())]*4)
#        for n in range(len(brain.G.nodes())):
#            node = brain.G.nodes()[n]
#            for x in range(3):
#                visarray[x:x+1,n:n+1] = float(brain.G.node[node]['xyz'][x])#*2-posCorrectionDict[x]
#            visarray[3:4,n:n+1] = brain.G.node[node]['cluster']
#
#        visarrayhubs = np.array([[0.0]*len(brain.hubs)]*4)
#        for n in range(len(brain.hubs)):
#            hub = brain.hubs[n]
#            for x in range(3):
#                visarrayhubs[x:x+1,n:n+1] = float(brain.G.node[hub]['xyz'][x])#*2-posCorrectionDict[x]
#            visarrayhubs[3:4,n:n+1] = brain.G.node[hub]['cluster']
#            
#        
#        try:
#            degennodes = [v for v in brain.G.nodes() if brain.G.node[v]['degenerating']]
#            visarraydegennodes = np.array([[0.0]*len(degennodes)]*4)
#            for n in range(len(degennodes)):
#                degennode = degennodes[n]
#                for x in range(3):
#                    visarraydegennodes[x:x+1,n:n+1] = float(brain.G.node[degennode]['xyz'][x])#*2-posCorrectionDict[x]
#                visarraydegennodes[3:4,n:n+1] = brain.G.node[degennode]['cluster']
#        except:
#            pass
#        
#        
#        edges = brain.G.edges()
#        
#    
##    mlab.figure(1,bgcolor=(1,1,1),size=(1200,900))
#    
#    pts = mlab.points3d(visarray[0],visarray[1],visarray[2],visarray[3],scale_factor=1,scale_mode='none',colormap='Greens',opacity=0,figure=scene)
#    pts.mlab_source.dataset.lines = np.array(edges)
#    
#    ptshubs = mlab.points3d(visarrayhubs[0],visarrayhubs[1],visarrayhubs[2],scale_factor=1,color=(1,0,0),figure=scene)
#    try:
#        ptsdegens = mlab.points3d(visarraydegennodes[0],visarraydegennodes[1],visarraydegennodes[2],scale_factor=1,color=(1,1,0),figure=scene)
#    except:
#        pass
#    
#    tube = mlab.pipeline.tube(pts,tube_radius=0.03,figure=scene)
#    mlab.pipeline.surface(tube,colormap='Accent',opacity=1)    
##    mlab.show()
#    # set a few environmental variables
#    scene.scene.camera.position = [145, 145, 130]
#    scene.scene.camera.focal_point = [-100, -100, -100]
#    scene.scene.camera.view_angle = 30.0
#    scene.scene.camera.view_up = [0.0, 0.0, 1.0]
#    scene.scene.camera.clipping_range = [100, 600]
#    scene.scene.camera.compute_view_plane_normal()
#    scene.scene.set_size((1600,1200))
#    scene.scene.render()
#    ## scene.scene.x_minus_view()
#    scene.scene.save(outfile+'.png')
#    scene.scene.save_x3d(outfile+".x3d")
#    
    
          
    
