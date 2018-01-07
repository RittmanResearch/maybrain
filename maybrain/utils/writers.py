from os import path, rename, remove
from csv import writer, DictWriter
from networkx.algorithms import centrality
from networkx import degree_histogram
from numpy import array, sqrt, mean, sum, shape, zeros


def output_adj_matrix(brain, filename):
    """
    Outputs the adjacency matrix to file

    brain: Instance of Brain class
    filename: the filename to which the adjacency matrix will be written
    """

    try:
        with open(filename, "w") as f:
            for line in brain.adjMat:
                for num in line[:-1]:
                    f.write(str(num) + '\t')
                f.write(str(line[-1]) + '\n')
    except IOError as error:
        error.strerror = 'Problem with opening file "' + filename + '": ' + error.strerror
        raise error


def output_edges(brain, filename, properties=[]):
    """
    Outputs the edges of a brain to file

    brain: Instance of Brain class
    filename: the filename to which the edges will be written
    properties: the list of properties you want to save from each edge
    """
    try:
        with open(filename, "w") as f:
            # write column headers
            line = 'n1' + '\t' + 'n2'
            for p in properties:
                line = line + '\t' + p
            line = line + '\n'
            f.write(line)

            for e in brain.G.edges(data=True):
                # add coordinates
                line = str(e[0]) + '\t' + str(e[1])
                # add other properties
                for p in properties:
                    try:
                        line = line + '\t' + str(e[2][p])
                    except:
                        continue  # Property doesn't exist, just continue
                line = line + '\n'
                # write out
                f.write(line)
    except IOError as error:
        error.strerror = 'Problem with opening file "' + filename + '": ' + error.strerror
        raise error
