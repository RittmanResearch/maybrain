class highlightObj():
    ''' object to hold information to highlight a subsection of a brainObj '''
    
    def __init__(self, nodes = [], edges = []):
        ''' Points refer to the indices of node coordinates in the brain object to which it
        is related. Edges is a set of pairs of coordinates of the same brain object '''
        
#        self._mode = 'pe' # p, e or pe for points only, edges only or points and edges
        self.nodeIndices = nodes # indices of points used from a brain object
        self.edgeIndices = edges # list of ordered pairs of edges
        self.colour = (0.5, 0.5, 0.5)
        self.opacity = 1.0
        self.edgeOpacity = None
       
        
    def getEdgeCoordsToPlot(self, brain):
        ''' turn list of edges into lists of coordinates - note that the second set of coordinates are the 
        vector from the first to the second points of the edge '''
        
        # initialise outputs
        x1 = []
        x2 = []
        y1 = []
        y2 = []
        z1 = []
        z2 = []
        s = []
                
        for e in self.edgeIndices:
            # get coordinates of edge indices
            p1 = brain.G.node[e[0]]['xyz']
            p2 = brain.G.node[e[1]]['xyz']    
            
            x1.append(p1[0])
            x2.append(p2[0]-p1[0]) 
            
            y1.append(p1[1])
            y2.append(p2[1]-p1[1])
            
            z1.append(p1[2])
            z2.append(p2[2]-p1[2])
            
            s.append(brain.G.edge[e[0]][e[1]]['weight'])
            
        return x1, y1, z1, x2, y2, z2, s
        
    
    def getCoords(self, brain):
        ''' turn indices of points into coordinates '''

        x = []
        y = []
        z = []        
        
        for p in self.nodeIndices:
            pc = brain.G.node[p]['xyz']
            # pc = brain.G.node[p[1]]['xyz']
            x.append(pc[0])
            y.append(pc[1])
            z.append(pc[2])
            
        return x, y, z
        
    
    def hlPrint(self, brain):
        ''' print the contents of the highlight object '''
        
        # introduction
        print('# ================ #')
        print('\n'+'highlight information:')
        # colour
        print('colour', self.colour)
        # opacity
        print('opacity', self.opacity) 
        # edge opacity
        print('edge opacity', self.edgeOpacity)
        
        # nodes
        print('\n nodes')
        for n in self.nodeIndices:
            print(str(n), brain.G.node[n]['xyz'])
            
        # edges
        print('\n edges')
        for e in self.edgeIndices:
            print(str(e), brain.G.node[e[0]]['xyz'], brain.G.node[e[1]]['xyz'])
            
        print('\n # ================ # \n')
