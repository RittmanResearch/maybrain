import numpy as np
import random


def modularity(brain, hierarchy=False, diag_val=0., nodes_to_exclude=None):
    """
    Modularity function borrowed (after asking nicely!) from
    https://sites.google.com/site/bctnet/measures/list and converted from
    matlab to python code.

    The main modification is to allow NA values in the association matrix.
    The code is being integrated in to maybrain: http://code.google.com/p/maybrain/

    The function only returns a hierarchical dictionary of matrices and
    modularities if hierarchy is True. Otherwise, labels are added to
    individual nodes and the modularity is assigned as 'Q', eg brain.Q
    """

    w = brain.adjMat.copy()
    n0 = len(w)  # number of nodes

    w = np.ma.array(w, mask=False)  # convert to masked array
    w.mask = w.data
    w.mask = False
    w[np.isnan(brain.adjMat)] = 0.

    h = 0  # hierarchy index
    ci = {h: np.ma.array(np.zeros(n0), mask=False,
                         dtype=int)}  # create dictionary of hierarchy assignments and blank arrays
    if nodes_to_exclude:
        ci[h].mask = ci[h].data
        ci[h].mask = False
        for i in [int(v) for v in nodes_to_exclude]:
            ci[h].mask[i] = True
            w.mask[i, :] = True
            w.mask[:, i] = True

    # change diagonals to d only in non-masked rows/columns and assign
    # initial values
    count = 0
    for i in range(n0):
        if np.ma.is_masked(ci[h][i]):
            pass
        else:
            ci[h][i] = int(count)
            count += 1
            w[i, i] = diag_val
            w.mask[i, i] = False

    q = {h: -1}

    # get rid of nan's
    w = w[np.invert(w.mask)]
    w.shape = np.repeat(np.sqrt(len(w)), 2)
    n = len(w)

    s = np.sum(w)  # weight of edges

    while 1:
        k = np.sum(w, axis=1)  # node degree
        km = k.copy()  # module degree
        knm = w.copy()  # node-to-module degree

        m = np.array([v for v in range(n)])  # initial module assignments

        nm = np.ones(n)  # number of nodes in modules

        flag = True  # flag for within network hierarchy search

        while flag:
            flag = False
            nlist = [v for v in range(n)]
            random.shuffle(nlist)
            while nlist:
                i = nlist.pop()
                dq = (knm[i, :] - knm[i, m[i]] + w[i, i]) - k[i] * (km - km[m[i]] + k[i]) / s  # algorithm condition
                #            dQ=(Knm(i,:)-Knm(i,M(i))+W(i,i)) - K(i).*(Km-Km(M(i))+K(i))/s;

                dq[m[i]] = 0

                max_dq = np.max(dq)  # find maximal increase in modularity

                if max_dq > 0:  # if maximal increase is positive
                    j = np.argmax(dq)

                    knm[:, j] = knm[:, j] + w[:, i]  # change node-to-module degrees
                    knm[:, m[i]] = knm[:, m[i]] - w[:, i]

                    km[j] = km[j] + k[i]
                    km[m[i]] = km[m[i]] - k[i]  # change module degrees

                    nm[j] += 1  # change number of nodes in modules
                    nm[m[i]] -= 1

                    m[i] = j  # reassign module

                    flag = True

        x, m1 = np.unique(m, return_inverse=True)

        h += 1
        ci[h] = np.ma.array(np.zeros(n0), dtype=int)

        for i in range(n):
            ci[h][ci[h - 1] == i] = int(m[i])
        ci[h].mask = ci[0].mask.copy()

        n = len(x)  # new number of modules

        w1 = np.zeros((n, n))  # new weighted matrix

        for i in range(n):
            for j in range(i, n):  # pool weights of nodes in same module w=sum(sum(W(M1==i,M1==j)));
                a = np.zeros(w.shape)
                ind_row = np.array([z for z, v in enumerate(m1) if v == i])
                ind_col = np.array([z for z, v in enumerate(m1) if v == j])

                for x in ind_row:
                    for y in ind_col:
                        a[x, y] = w[x, y]

                w = np.sum(a)
                #                print w

                w1[i, j] = w
                w1[j, i] = w

        w = w1.copy()
        del w1

        q[h] = np.sum(np.diagonal(w)) / s - np.sum(np.sum(w / s, axis=0) ** 2)  # compute modularity
        if q[h] <= q[h - 1]:  # if modularity does not increase
            break

    for node in brain.G.nodes():
        brain.G.node[node]['module'] = ci[h - 1][node]

    brain.Q = q[h - 1]

    # return hierarchy only if desired
    if hierarchy:
        return ci, q


def within_module_degree(brain, ci, weight=None):
    """
    To calculate mean within module degree
    """
    #    moduleList = [v for v in set(ci.values())] # get module list
    withinDegDict = {}  # output dictionary of nodes and mean within module degree
    #    modDict = {m:[v for v in ci.keys() if ci[v]==m] for m in moduleList} # sort nodes in to modules

    for n in brain.G.nodes():
        m = ci[n]  # select module

        eList = brain.G.edges([n])
        eList = [e for e in eList if all([ci[e[0]] == m, ci[e[1]] == m])]  # find edges exclusively within the module

        if weight:
            wts = np.sum([float(G.edge[e[0]][e[1]]['weight']) for e in eList])  # get weights/degree
            wts = wts / float(len(eList))
        else:
            wts = float(len(eList))

        withinDegDict[n] = wts

    #        if len(modDict[m]) > 1:
    #            withinDegDict[n] = wts/(len(modDict[m])-1)  # mean of weight/degree, ie average degree within module
    #        else:
    #            withinDegDict[n] = 0.

    return withinDegDict
