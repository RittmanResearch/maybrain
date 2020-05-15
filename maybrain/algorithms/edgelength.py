import numpy as np

def edgelength(G, node_wise=False, edge_wise=False, summary="mean"):
    """
    This function calculates the physical distance between pairs of nodes. The default behaviour is to
    return a dictionary of edges.

    nodeWise:    if True, then returns a dictionary of the sum of distance of all edges for each node
    edgewise:    if True returns a dictionary of all edge values
    summary:     values are either "mean" or "median", returns a single value for all edges
    """
    el =  dict(zip(G.edges(),
                   [ np.absolute(np.linalg.norm(np.array(G.node[edge[0]]['xyz']) - np.array(G.node[edge[1]]['xyz']))) for edge in G.edges()]))
    
    if node_wise:
        eln = {}
        for n in G.nodes():
            eln[n] = np.sum([el[e] for e in el.keys() if n in e])
        return(eln)
    else:
        if edge_wise:
            return(el)
        else:
            if summary=="mean":
                return(np.mean([v for v in el.values()]))
            elif summary=="median":
                return(np.median([v for v in el.values()]))

