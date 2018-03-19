# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 19:05:26 2014
Module for linking with the Allen brain atlas
@author: tim
"""

import csv
from os import path, rename, remove
from glob import glob
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

from maybrain import brain
# from maybrain.plotting import mayavi_wrapper as plot

class AllenBrain:
    """
    An object the combines a network generated from imaging data with gene
    expression data taken from the Allen brain atlas. Genetic data can be
    downloaded from http://human.brain-map.org/static/download.

    An example will be available on the wiki soon.
    """
    def __init__(self, allen_subj,
                 assoc_matrix,
                 allen_dir="allen_data",
                 delim=None,
                 imaging_spatial_file="atlas471_xyz_flip_xy.txt",
                 nodesToExclude=[],
                 symmetrise=False,
                 mirror=False,
                 signif_value=1.0,
                 convertMNI=False):
        """
        This object contains two 'brain' network objects, one for the imaging
        data and one for the Allen data. Embedded functions make a comparison
        between the imaging and Allen data by pairing nodes between the two
        network objects.
        """

        self.allen_dir = allen_dir # directory where data is located
        self.allen_subj = allen_subj # sujbect number, eg 178236545

        # files in the subject directory to assign gene properties
        self.annotation_file = "SampleAnnot.csv"
        self.microarray_expression_file = "MicroarrayExpression.csv"
        self.probe_file = "Probes.csv"

        self.signif_value = signif_value # set significance value for gene correlations

        # if symmetrise is true then regions are identified by the structure
        # name, if symmetrise is false, then regions are identified by the
        # structural acronym these refer to headers in the SampleAnnot.csv file
        if symmetrise:
            self.s_lab = "structure_acronym"
        else:
            self.s_lab = "structure_name"

        self.mirror = mirror
        if symmetrise and self.mirror:
            print("Please select either a symmetrised or mirrored graph, \
                  ignoring mirror=True")
            self.mirror = False

        # set up brain for expression data
        self.expr_brain = brain.Brain()

        node_counter = 0

        # import probe data
        file = open(path.join(self.allen_dir, self.allen_subj, self.annotation_file), "r")
        reader = csv.DictReader(file, delimiter=",", quotechar='"')

        self.headers = ['probe']
        self.s_id_dict = {}

        for line in reader:
            s_id = line[self.s_lab]
            node_count = node_counter # GIVE NODES UNIQUE INCREMENTAL sa_id (across all subjects)
            print(node_count)
            if not self.expr_brain.G.has_node(node_count):
                self.expr_brain.G.add_node(node_count)
                brain.nx.set_node_attributes(self.expr_brain.G, {node_count: line})

                # store the structure_acronym/structure_name for the node
                brain.nx.set_node_attributes(self.expr_brain.G, {node_count: {'s_id': s_id}})
#                node_counter += 1

                #STORE structure_acronym or structure_name depending on symmetrise
                self.headers.append(s_id)

                if not s_id in list(self.s_id_dict.keys()):
                    self.s_id_dict[s_id] = [node_count]
                else:
                    self.s_id_dict[s_id].append(node_count)
            node_counter += 1
        file.close()

        if convertMNI:
            # convert location data for Allen brain from MNI space
            for node_count in self.expr_brain.G.nodes():
                x = 45 - (float(self.expr_brain.G.node[node_count]['mni_x'])/2)
                y = 63 + (float(self.expr_brain.G.node[node_count]['mni_y'])/2)
                z = 36 + (float(self.expr_brain.G.node[node_count]['mni_z'])/2)
                self.expr_brain.G.node[node_count]['xyz'] = (x, y, z)
        else:
            for node_count in self.expr_brain.G.nodes():
                x = float(self.expr_brain.G.node[node_count]['mni_x'])
                y = float(self.expr_brain.G.node[node_count]['mni_y'])
                z = float(self.expr_brain.G.node[node_count]['mni_z'])
                self.expr_brain.G.node[node_count]['xyz'] = (x, y, z)

        # copy hemisphere if required
        if self.mirror and len(self.expr_brain.G.nodes()) < 600:
            self.expr_brain.copy_hemisphere()

        # set up brain with graph properties
        self.imaging_brain = brain.Brain()
        self.imaging_brain.import_adj_file(assoc_matrix, delimiter=delim,
                                           nodes_to_exclude=nodesToExclude)
        self.imaging_brain.import_spatial_info(imaging_spatial_file)

    def comparison(self):
        """
        set up dictionary to link nodes from probe data and graph
        keys are the mri nodes and values are disctionaries containing two keys:
        key 1= allen, value= (node_count=sa_id of closest allen node, dist=distance to closest allen node)
        key 2= mri, value= (node_count=sa_id of closest other mri node, dist=distance to closest mri node)
        """
        imaging_node_dict = {}

        # iterate through the imaging nodes and find the closest Allen brain nodes
        for node in self.imaging_brain.G.nodes():
            d_other = (None, 999.) # dummy length of 999
            d_own = (None, 999.)

            # iterate through the Allen nodes
            for node_count in self.expr_brain.G.nodes():
                # find the distance between the imaging and Allen nodes
                dist = np.linalg.norm(np.array(self.imaging_brain.G.node[node]['xyz']
                                               - np.array(self.expr_brain.G.node[node_count]['xyz'])))
                if dist < d_other[1]:
                    d_other = (node_count, dist)

            # now iterate through the imaging nodes to find the closest imaging node
            for node_count in [v for v in self.imaging_brain.G.nodes() if not v == node]:
                dist = np.linalg.norm(np.array(self.imaging_brain.G.node[node]['xyz']
                                               - np.array(self.imaging_brain.G.node[node_count]['xyz'])))
                if dist < d_own[1]:
                    d_own = (node_count, dist)
            imaging_node_dict[node] = {"allen":d_other, "MRIs":d_own}

        # set up dictionary to link nodes from probe data and graph
        allen_node_dict = {}
        for node in self.expr_brain.G.nodes():
            d_other = (None, 999.)
            d_own = (None, 999.)

            # iterate through the imaging nodes to find the closest nodes
            for node_count in self.imaging_brain.G.nodes():
                dist = np.linalg.norm(np.array(self.expr_brain.G.node[node]['xyz']
                                               - np.array(self.imaging_brain.G.node[node_count]['xyz'])))
                if dist < d_other[1]:
                    d_other = (node_count, dist)

            # iterate through the allen nodes to find closest nodes
            for node_count in [v for v in self.expr_brain.G.nodes() if not v == node]:
                dist = np.linalg.norm(np.array(self.expr_brain.G.node[node]['xyz']
                                               - np.array(self.expr_brain.G.node[node_count]['xyz'])))
                if dist < d_own[1]:
                    d_own = (node_count, dist)
            allen_node_dict[node] = {"allen":d_own, "MRIs":d_other}

        node_pairs = [] # create a list of node pairs

        # for each MRI node..
        for node in list(imaging_node_dict.keys()):
            # ...find closest allen node 'node_count'
            node_count = imaging_node_dict[node]['allen'][0]
            self.imaging_brain.G.node[node]['pair'] = node_count
            self.expr_brain.G.node[node_count]['pair'] = node
            node_pairs.append((node, node_count))

        # iterate through the list of nodes to assign pairs and remove unpaired nodes
        node_list = [v for v in self.expr_brain.G.nodes()]
        for node in node_list:
            if 'pair' not in list(self.expr_brain.G.node[node].keys()):
                self.expr_brain.G.remove_node(node)

    def probedata(self, prop_dict, graph_metric_name="gm", node_list=None, scatter_plot=False,
                  probe_list=[], probe_numbers=[], thresh=False):
        '''
        Explanation required.
        '''
        probes = open(path.join(self.allen_dir, self.allen_subj, "Probes.csv"), "r")
        probe_reader = csv.DictReader(probes, delimiter=",",
                                      quotechar='"')
        probe_dict = {line['probe_id']:[line['gene_symbol'], line['gene_name']] for line in probe_reader}
        if probe_list:
            probe_numbers = []
            for probe in probe_list:
                probe_numbers.extend([v for v in list(probe_dict.keys()) if any([probe in probe_dict[v][1], probe in probe_dict[v][0]])])
            print((" ".join(["Probe numbers:", ' '.join(probe_numbers)])))

        else:
            probe_numbers = None

        out_file = path.join(self.allen_subj, graph_metric_name+'.txt')
        print(("Saving data in:"+out_file))
        if path.exists(out_file):
            rename(out_file, out_file+'.old')
        file = open(out_file, "w")
        file.writelines(', '.join(['probe_id', 'gene_name', 'r', 'probe\node_count']))
        file.close()

        if node_list:
            for node in self.imaging_brain.G.nodes():
                if not node in node_list:
                    self.expr_brain.G.remove_node(self.imaging_brain.G.node[node]['pair'])
                    self.imaging_brain.G.remove_node(node)

        # import probe data
        file = open(path.join(self.allen_dir, self.allen_subj, self.microarray_expression_file), "r")
        reader = csv.DictReader(file, delimiter=",",
                                fieldnames=self.headers,
                                quotechar='"')

        for line in reader:
            if thresh:
                self.probesubt(line, prop_dict, probe_dict, graph_metric_name, out_file, probe_numbers)
            else:
                self.probesub(line, prop_dict, probe_dict, graph_metric_name, out_file,
                              scatter_plot, probe_numbers)
        file.close()

    def xmatrix(self, out_file="Xmatrix.csv", probe_numbers=None, temp_mat_name="tempMat.txt",
                st_dev=False, st_dev_file="NodesSd.txt"):
        """
        Needs writing
        """
        # get all probes if otherwise unspecified
        if not probe_numbers:
            file = open(path.join(self.allen_dir, self.allen_subj, self.probe_file))
            probe_numbers = [line['probe_id'] for line in csv.DictReader(file, delimiter=",", quotechar='"')]
            file.close()

        # set up out file
        out = open(self.allen_subj+out_file, "wb")
        headers = ["Gene"]
        headers.extend([str(v) for v in self.imaging_brain.G.nodes()])
        writer = csv.DictWriter(out, fieldnames=headers, delimiter=" ")
        writer.writeheader()

        # set up matrix
        y = len(probe_numbers)
        z = len(self.imaging_brain.G.nodes())
        st = np.max([len(list(self.s_id_dict.values()))]) # max numbers of nodes for any region

        probe_mat = np.memmap(temp_mat_name,
                              dtype="float64",
                              mode="w+",
                              shape=(y, z, st))

        # set up gene list
        gene_flag = True
        probe_file = open(path.join(self.allen_dir, self.allen_subj, self.probe_file))
        probe_reader = csv.DictReader(probe_file, delimiter=",", quotechar='"')
        probe_dict = {line['probe_id']:line['gene_symbol'] for line in probe_reader}
        gene_list = list(probe_dict.values())
        set(gene_list)
        gene_list = {gene:[] for gene in gene_list}

        s_id_list = [v for v in self.headers if v != 'probe']
        s_id_list = set(s_id_list)

        # get the corresponding node names in the MRI graph
        # c_nodes is a dict whose keys are all the MRI nodes and values are the matched alen nodes
        #### PV modified line below which constructed c_nodes by looping through allen nodes
        # but with PV's lax matching criteria several mri nodes can be matched to same allen node
        # the mri pair of these allen nodes gets overwritten in self.expr_brain.G.nodes and so
        # not all mri nodes will appear as pairs of allen nodes in this dict...
        # need to look up pairs in self.imaging_brain.G.nodes instead, where
        # each mri node is matched to an allen region
        # c_nodes = {str(self.expr_brain.G.node[v]['pair']):v for v in self.expr_brain.G.nodes()}

        c_nodes = {str(v):self.imaging_brain.G.node[v]['pair'] for v in self.imaging_brain.G.nodes()}

        # assign values to matrix
        print((str(self.allen_subj)))
        print('\node_count')
        # Generate custom fieldnames list for DictReader which doesn'thresh rely on structure_id
        # *************************************
        fieldnames_pv = ['probe']
        my_node_dict = {}
        temp_headers = self.headers
        temp_headers.remove('probe')

        for p, q in enumerate(self.headers):
            my_node_dict[p] = q
            fieldnames_pv.append(str(p))   #####
        # *************************************

        # import probe data
        file = open(path.join(self.allen_dir, self.allen_subj, self.microarray_expression_file), "r")
        reader = csv.DictReader(file, delimiter=",",
                                fieldnames=fieldnames_pv,
                                quotechar='"')
        y = 0
        for line in reader:
            probe = line['probe']
            if probe in probe_numbers:
                # assign probe values to sample numbers
                for z, c_node in enumerate(self.imaging_brain.G.nodes()):
                    a_node = c_nodes[str(c_node)]
                    # *************************************
                    # Find structure_acronym corresponding to the matched allen node
                    acronym = self.expr_brain.G.node[a_node][self.s_lab]
                    # Initialise list for summing expression values for allen nodes with
                    # same structure_acronym
                    total_expression = []
                    # Loop over allen nodes to find all those with the correct acronym for
                    # the current MRI node
                    # and add their expression values to the list for averaging
                    s = 0
                    for sa_id, struct_ac in list(my_node_dict.items()):
                        if struct_ac == acronym:
                            # print(sa_id, struct_ac, line[str(sa_id)])
                            # print('\node_count')
                            # if line[str(sa_id)]:
                            total_expression.append(float(line[str(sa_id)]))
                            probe_mat[y, z, s] = float(line[str(sa_id)])
                            s += 1

                if gene_flag:
                    # records the position of the probe in a dictionary with genes as a key
                    gene_list[probe_dict[probe]].append(y)
                y += 1
        file.close()
        del reader

        # get values to normalise expression levels for each probe within subject
        for y in range(probe_mat.shape[0]):
            # create a masked array removing the 0. values
            subj_mat = np.ma.array(probe_mat[y, :, st:st],
                                   mask=probe_mat[y, :, st:st] == 0.,
                                   dtype="float64")

            subj_mat = (subj_mat - np.mean(np.ma.array(subj_mat, mask=subj_mat == 0.))) / np.std(np.ma.array(subj_mat, mask=subj_mat == 0.))
            probe_mat[y, :, st:st] = subj_mat

        gene_flag = False

        # collapse across subjects and probes by gene
        gene_names = list(gene_list.keys())
        gene_names.sort() # sort in to alphabetical order

        # collapse across nodes within regions (averaging across all subjects)
        sh = probe_mat.shape
        probe_mat_temp = np.memmap("probe_mat_temp.txt", mode="w+", dtype="float64", shape=sh[:2])

        # write out the standard deviation for each probe if specified
        if st_dev:
            st_dev_out = open(st_dev_file, "wb")
            st_dev_out.writelines("Probe Node st_dev\node_count")
        for y in range(sh[0]):
            for z in range(sh[1]):
                # mask out unused values, ie where there are less than the maximum number
                # of homolous nodes in a structural region
                probe_mat_temp[y, z] = np.mean(np.ma.array(probe_mat[y, z],
                                                           mask=probe_mat[y, z] == 0.))
                if st_dev:
                    std = np.std(np.ma.array(probe_mat[y, z], mask=probe_mat[y, z] == 0.))
                    st_dev_out.writelines(' '.join([str(int(probe_numbers[y])),
                                                    str(self.imaging_brain.G.nodes()[z]),
                                                    "{:2.5f}".format(std)])+'\node_count')
        if st_dev:
            st_dev_out.close()

        # reassign the probe matrix and delete temporary memory-mapped file
        probe_mat = probe_mat_temp
        del probe_mat_temp
        remove(temp_mat_name)

        for gene in gene_names:
            if gene_list[gene]:
                x = probe_mat.shape[1] # number of nodes
                y = len(gene_list[gene]) # number of probes

                gene_mat = np.zeros(shape=(x, y), dtype="float64")

                for node_count, probe in enumerate(gene_list[gene]):
                    gene_mat[:, node_count] = probe_mat[probe, :]

                # collapse values across probes for each gene
                mean_gene = np.mean(np.ma.array(gene_mat, mask=np.isnan(gene_mat)), axis=1)

                out_dict = dict(list(zip([str(v) for v in self.imaging_brain.G.nodes()],
                                         ["{:10.20f}".format(v) for v in mean_gene])))
                out_dict["Gene"] = gene
                writer.writerow(out_dict)

        out.close()
        remove("probe_mat_temp.txt") # delete memory map file

    def probesub(self, line, prop_dict, probe_dict, scatter_plot,
                 graph_metric_name, out_file, probe_numbers=None):
        """
        Needs completing
        """
        if probe_numbers:
            if line['probe'] in probe_numbers:
                # assign probe values to sample numbers
                for node in self.expr_brain.G.nodes():
                    if line[self.expr_brain.G.node[node][self.s_lab]]:
                        self.expr_brain.G.node[node][line['probe']] = line[self.expr_brain.G.node[node][self.s_lab]]
                    else:
                        self.expr_brain.G.node[node][line['probe']] = None

                aa_mat = np.zeros((len(self.imaging_brain.G.nodes()), 3))

                for node_count, c_node in enumerate(self.imaging_brain.G.nodes()):
                    node = self.imaging_brain.G.node[c_node]['pair']
                    if prop_dict[self.expr_brain.G.node[node]['pair']]:
                        aa_mat[node_count, 0] = self.expr_brain.G.node[node][line['probe']]
                        aa_mat[node_count, 1] = prop_dict[c_node]
                        aa_mat[node_count, 2] = c_node
                    else:
                        aa_mat[node_count, :] = [np.nan, np.nan, np.nan]

                aa_mat = aa_mat[~np.isnan(aa_mat)]
                aa_mat.shape = (len(aa_mat)/3, 3)

                r, p = stats.pearsonr(aa_mat[:, 0], aa_mat[:, 1])

                if p < self.signif_value:
                    print(line['probe'])
                    # plot graph
                    out = open(out_file, "a")
                    out.writelines(', '.join([str(v) for v in [line['probe'],
                                                               '"' + probe_dict[line['probe']][1]+'"',
                                                               r, p]])+'\node_count')
                    out.close()
                    if scatter_plot:
                        plt.scatter(aa_mat[:, 1], aa_mat[:, 0])
                        plt.savefig(out_file.replace('.txt', line['probe']+'.png'), dpi=300)
                        plt.close()

                    # save data
                    print("Saving data in :"+out_file.replace('.txt',
                                                              line['probe']+graph_metric_name+'.txt'))
                    dat_file = open(out_file.replace('.txt',
                                                     line['probe']+graph_metric_name+'.txt'), "wb")
                    dat_file.writelines(' '.join([line['probe'],
                                                  graph_metric_name,
                                                  "node", "subj"])+'\node_count')
                    dat_file.writelines('\node_count'.join([' '.join([str(aa_mat[node_count, 0]),
                                                                      str(aa_mat[node_count, 1]),
                                                                      str(aa_mat[node_count, 2]),
                                                                      self.allen_subj]) for node_count in range(len(aa_mat[:, 1]))]))
                    dat_file.close()
        else:
            pass

    def probesubt(self, line, prop_dict, probe_dict, graph_metric_name, out_file, probe_numbers=None):
        '''
        line is a line from the probe file.
        The purpose of this function is to write thresholded data to a datafile
        eg for use in ANOVA
        '''
        probe = line['probe']
        dat_file = None
        if probe_numbers:
            if probe in probe_numbers:
                # assign probe values to sample numbers
                for node in self.expr_brain.G.nodes():
                    if line[str(node)]:
                        self.expr_brain.G.node[node][probe] = line[str(node)]
                    else:
                        self.expr_brain.G.node[node][probe] = None

                    out_dict = {probe:probe, 'subj':self.allen_subj}
                    for probe in list(prop_dict.keys()):
                        out_dict[probe] = prop_dict[probe]

                    if not dat_file:
                        headers = [probe, "subj"]
                        gm_subjs = list(prop_dict[list(probe_dict.keys())[0]].keys())
                        gm_subjs.sort()
                        headers.extend(gm_subjs)

                        dat_file = open(out_file.replace('.txt', probe+graph_metric_name+'.txt'), "wb")
                        writer = csv.DictWriter(dat_file, fieldnames=headers, delimiter=" ")
                        writer.writeheader()

                    writer.writerow(out_dict)

class multisubj:
    """
    This object contains two 'brain' network objects, one for the imaging data and
    one for the Allen data. Embedded functions make a comparison between the imaging
    and Allen data by pairing nodes between the two network objects.
    """
    def __init__(self, assoc_matrix, allen_dir="allen_data", nodesToExclude=[], delim=" ",
                 subj_list=None, spatial_file="parcel_500.txt", symmetrise=False,
                 convertMNI=False, mirror=True):

        self.allen_dir = allen_dir

        if subj_list:
            self.allen_subj_list = subj_list
        else:
            self.allen_subj_list = [v for v in glob("17823*") if path.isdir(path.join(self.allen_dir, v))]

        self.annotation_file = "SampleAnnot.csv"
        self.microarray_expression_file = "MicroarrayExpression.csv"
        self.probe_file = "Probes.csv"
        self.mirror = mirror

        # if symmetrise is true then regions are identified by the structure name,
        # if symmetrise is false, then regions are identified by the structural acronym
        # these refer to headers in the SampleAnnot.csv file
        if symmetrise:
            self.s_lab = "structure_acronym"
        else:
            self.s_lab = "structure_name"

        # set up brain for expression data
        self.expr_brain = brain.Brain()

        node_counter = 0

        self.headers = {}

        # dictionary storing the list of nodes for each structural sa_id by subject
        # for use later in averaging across all subjects
        self.s_id_dict = {}
        for subj in self.allen_subj_list:
            self.s_id_dict[subj] = {}
            # import probe data
            file = open(path.join(self.allen_dir, subj, self.annotation_file), "r")
            reader = csv.DictReader(file, delimiter=",", quotechar='"')

            self.headers[subj] = ['probe']
            for line in reader:
                # node_count = int(line["structure_id"])
                s_id = line[self.s_lab]
                node_count = node_counter # GIVE NODES UNIQUE INCREMENTAL sa_id (across all subjects)
                if not self.expr_brain.G.has_node(node_count):
                    self.expr_brain.G.add_node(node_count)
                    self.expr_brain.G.node[node_count] = line

                    # store the structure_acronym/structure_name for the node
                    self.expr_brain.G.node[node_count]['s_id'] = s_id
                    node_counter += 1
                    #STORE structure_acronym or structure_name depending on symmetrise
                    self.headers[subj].append(s_id)

                    if not s_id in list(self.s_id_dict[subj].keys()):
                        self.s_id_dict[subj][s_id] = [node_count]
                    else:
                        self.s_id_dict[subj][s_id].append(node_count)

            file.close()

            if convertMNI:
                # convert location data for Allen brain from MNI space
                for node_count in self.expr_brain.G.nodes():
                    x = 45 - (float(self.expr_brain.G.node[node_count]['mni_x'])/2)
                    y = 63 + (float(self.expr_brain.G.node[node_count]['mni_y'])/2)
                    z = 36 + (float(self.expr_brain.G.node[node_count]['mni_z'])/2)
                    self.expr_brain.G.node[node_count]['xyz'] = (x, y, z)
            else:
                for node_count in self.expr_brain.G.nodes():
                    x = float(self.expr_brain.G.node[node_count]['mni_x'])
                    y = float(self.expr_brain.G.node[node_count]['mni_y'])
                    z = float(self.expr_brain.G.node[node_count]['mni_z'])
                    self.expr_brain.G.node[node_count]['xyz'] = (x, y, z)

            if self.mirror and len(self.expr_brain.G.nodes()) < 600:
                self.expr_brain.copy_hemisphere()

        # set up brain with graph properties
        self.imaging_brain = brain.Brain()
        self.imaging_brain.import_adj_file(assoc_matrix,
                                           delimiter=delim,
                                           nodes_to_exclude=nodesToExclude)
        self.imaging_brain.import_spatial_info(spatial_file)

    def comparison(self):
        """
        set up dictionary to link nodes from probe data and graph
        keys are the mri nodes and values are dictionaries containing two keys:
        key 1= allen, value= (node_count=sa_id of closest allen node, dist=distance to closest allen node)
        key 2= mri, value= (node_count=sa_id of closest other mri node, dist=distance to closest mri node)
        """
        imaging_node_dict = {}

        # iterate through the imaging nodes and find the closest Allen brain nodes
        for node in self.imaging_brain.G.nodes():
            d_other = (None, 999.) # dummy length of 999
            d_own = (None, 999.)

            # iterate through the Allen nodes
            for node_count in self.expr_brain.G.nodes():
                # find the distance between the imaging and Allen nodes
                dist = np.linalg.norm(np.array(self.imaging_brain.G.node[node]['xyz']
                                               - np.array(self.expr_brain.G.node[node_count]['xyz'])))
                if dist < d_other[1]:
                    d_other = (node_count, dist)

            # now iterate through the imaging nodes to find the closest imaging node
            for node_count in [v for v in self.imaging_brain.G.nodes() if not v == node]:
                dist = np.linalg.norm(np.array(self.imaging_brain.G.node[node]['xyz']
                                               - np.array(self.imaging_brain.G.node[node_count]['xyz'])))
                if dist < d_own[1]:
                    d_own = (node_count, dist)
            imaging_node_dict[node] = {"allen":d_other, "MRIs":d_own}

        # iterate through the allen nodes to find the closest imaging node
        allen_node_dict = {}
        for node in self.expr_brain.G.nodes():
            d_other = (None, 999.)
            d_own = (None, 999.)

            # iterate through the imaging nodes
            for node_count in self.imaging_brain.G.nodes():
                dist = np.linalg.norm(np.array(self.expr_brain.G.node[node]['xyz']
                                               - np.array(self.imaging_brain.G.node[node_count]['xyz'])))
                if dist < d_other[1]:
                    d_other = (node_count, dist)

            # iterate through the allen nodes to find the closes node
            for node_count in [v for v in self.expr_brain.G.nodes() if not v == node]:
                dist = np.linalg.norm(np.array(self.expr_brain.G.node[node]['xyz']
                                               - np.array(self.expr_brain.G.node[node_count]['xyz'])))
                if dist < d_own[1]:
                    d_own = (node_count, dist)
            allen_node_dict[node] = {"allen":d_own, "MRIs":d_other}

        node_pairs = [] # create a list of node pairs

        # for each MRI node
        for node in list(imaging_node_dict.keys()):
            # find closest allen node 'node_count'
            node_count = imaging_node_dict[node]['allen'][0]
            self.imaging_brain.G.node[node]['pair'] = node_count
            self.expr_brain.G.node[node_count]['pair'] = node
            node_pairs.append((node, node_count))

        for node in self.expr_brain.G.nodes():
            if 'pair' not in list(self.expr_brain.G.node[node].keys()):
                self.expr_brain.G.remove_node(node)

    def comparisonaveraged(self):
        """
        This function should generate sets of nodes in the imaging data associated
        with single nodes in the Allen data, ie all the closest imaging data nodes will
        be associated with any specific Allen node.
        """
        for node_count in self.expr_brain.G.nodes():
            self.expr_brain.G.node[node_count]['pairNodes'] = []

        # iterate through imaging nodes to find closes Allen node
        for node in self.imaging_brain.G.nodes():
            d_other = (None, 999.) # dummy length of 999

            for node_count in self.expr_brain.G.nodes():
                dist = np.linalg.norm(np.array(self.imaging_brain.G.node[node]['xyz']
                                               - np.array(self.expr_brain.G.node[node_count]['xyz'])))
                if dist < d_other[1]:
                    d_other = (node_count, dist)

            self.expr_brain.G.node[d_other[0]]['pairNodes'].append(node)

        for node in self.expr_brain.G.nodes():
            if not self.expr_brain.G.node[node]['pairNodes']:
                self.expr_brain.G.remove_node(node)

    def probedata(self, probe_numbers=[], mean_vals=True):
        """
        If mean_vals is specified, this takes the mean probe value across subjects
        """
        for subj in self.allen_subj_list:
            # import probe data
            file = open(path.join(self.allen_dir, subj, self.microarray_expression_file), "r")
            reader = csv.DictReader(file, delimiter=",",
                                    fieldnames=self.headers[subj],
                                    quotechar='"')
            for line in reader:
                probe = line['probe']
                if probe in probe_numbers:
                    # assign probe values to sample numbers
                    for c_node in self.imaging_brain.G.nodes():
                        node = self.imaging_brain.G.node[c_node]['pair']
                        s_id = self.expr_brain.G.node[node][self.s_lab]
                        if not probe in list(self.expr_brain.G.node[node].keys()):
                            self.expr_brain.G.node[node][probe] = {}

                        if s_id in list(line.keys()):
                            self.expr_brain.G.node[node][probe][subj] = line[s_id]
            file.close()
            del reader

        # self.expr_brain.G.nodes is a dict containing every
        # UNIQUE structure sa_id (across all subjects)
        if mean_vals:
            for node_count in self.expr_brain.G.nodes():
                for probe in probe_numbers:
                    self.expr_brain.G.node[node_count][probe] = np.mean([float(v) for v in list(self.expr_brain.G.node[node_count][probe].values())])

    def xmatrix(self, out_file="Xmatrix.csv", allen_dir="allen_data",
                probe_numbers=None, temp_mat_name="tempMat.txt", st_dev=False,
                st_dev_file="NodesSd.txt"):
        """
        Needs writing
        """
        # get all probes if otherwise unspecified
        if not probe_numbers:
            file = open(path.join(allen_dir, self.allen_subj_list[0], self.probe_file))
            reader = csv.DictReader(file, delimiter=",", quotechar='"')
            probe_numbers = [line['probe_id'] for line in reader]
            del reader
            file.close()

        # set up out file
        out = open(out_file, "wb")
        headers = ["Gene"]
        headers.extend([str(v) for v in self.imaging_brain.G.nodes()])
        writer = csv.DictWriter(out, fieldnames=headers, delimiter=" ")
        writer.writeheader()

        # set up matrix
        x = len(self.allen_subj_list)
        y = len(probe_numbers)
        z = len(self.imaging_brain.G.nodes())
        
        # max numbers of nodes for any region
        st = np.max([np.max([len(v) for v in list(self.s_id_dict[subj].values())]) for subj in list(self.s_id_dict.keys())])

        probe_mat = np.memmap(temp_mat_name,
                              dtype="float64",
                              mode="w+",
                              shape=(y, z, st*x))

        # set up gene list
        gene_flag = True
        probe_file = open(path.join(allen_dir, self.allen_subj_list[0], self.probe_file))
        probe_reader = csv.DictReader(probe_file, delimiter=",", quotechar='"')
        probe_dict = {line['probe_id']:line['gene_symbol'] for line in probe_reader}
        gene_list = list(probe_dict.values())
        set(gene_list)
        gene_list = {gene:[] for gene in gene_list}

        s_id_list = []
        for subj in self.allen_subj_list:
            s_id_list.extend([v for v in self.headers[subj] if not v == 'probe'])
        s_id_list = set(s_id_list)

        # get the corresponding node names in the MRI graph
        # c_nodes is a dict whose keys are all the MRI nodes and values are the matched alen nodes
        #### PV modified line below which constructed c_nodes by looping through allen nodes
        # but with PV's lax matching criteria several mri nodes can be matched to same allen node
        # the mri pair of these allen nodes gets overwritten in self.expr_brain.G.nodes and so not
        # all mri nodes will appear as pairs of allen nodes in this dict... need to look up pairs
        # in self.imaging_brain.G.nodes instead, where each mri node is matched to an allen region
        # c_nodes = {str(self.expr_brain.G.node[v]['pair']):v for v in self.expr_brain.G.nodes()}

        c_nodes = {str(v):self.imaging_brain.G.node[v]['pair'] for v in self.imaging_brain.G.nodes()}

        # assign values to matrix
        for x, subj in enumerate(self.allen_subj_list):
            print((str(subj)))
            print('\node_count')
            # Generate custom fieldnames list for DictReader which doesn'thresh rely on structure_id
            # *************************************
            fieldnames_pv = ['probe']
            my_node_dict = {}
            temp_headers = self.headers[subj]
            temp_headers.remove('probe')

            for p, q in enumerate(self.headers[subj]):
                my_node_dict[p] = q
                fieldnames_pv.append(str(p))   #####
            # *************************************

            # import probe data
            file = open(path.join(self.allen_dir, subj, self.microarray_expression_file), "r")
            reader = csv.DictReader(file, delimiter=",",
                                    fieldnames=fieldnames_pv,
                                    quotechar='"')
            y = 0
            for line in reader:
                probe = line['probe']
                if probe in probe_numbers:
                    # assign probe values to sample numbers
                    for z, c_node in enumerate(self.imaging_brain.G.nodes()):
                        a_node = c_nodes[str(c_node)]
                        # *************************************
                        # Find structure_acronym corresponding to the matched allen node
                        acronym = self.expr_brain.G.node[a_node][self.s_lab]
                        # Initialise list for summing expression values for
                        # allen nodes with same structure_acronym
                        total_expression = []
                        # Loop over allen nodes to find all those with the correct acronym
                        # for the current MRI node and add their expression values to
                        # the list for averaging
                        s = 0
                        for sa_id, struct_ac in list(my_node_dict.items()):
                            if struct_ac == acronym:
                                # print(sa_id, struct_ac, line[str(sa_id)])
                                # print('\node_count')
                                # if line[str(sa_id)]:
                                total_expression.append(float(line[str(sa_id)]))
                                probe_mat[y, z, st*x + s] = float(line[str(sa_id)])
                                s += 1

                    if gene_flag:
                        # records the position of the probe in a dictionary with genes as a key
                        gene_list[probe_dict[probe]].append(y)
                    y += 1
            file.close()
            del reader

            # get values to normalise expression levels for each probe within subject
            for y in range(probe_mat.shape[0]):
                # create a masked array removing the 0. values
                subj_mat = np.ma.array(probe_mat[y, :, st*x:st*x+st],
                                       mask=probe_mat[y, :, st*x:st*x+st] == 0.,
                                       dtype="float64")

                subj_mat = (subj_mat - np.mean(np.ma.array(subj_mat, mask=subj_mat == 0.))) / np.std(np.ma.array(subj_mat, mask=subj_mat == 0.))
                probe_mat[y, :, st*x:st*x+st] = subj_mat

            gene_flag = False

        # collapse across subjects and probes by gene
        gene_names = list(gene_list.keys())
        gene_names.sort() # sort in to alphabetical order

        # collapse across nodes within regions (averaging across all subjects)
        sh = probe_mat.shape
        probe_mat_temp = np.memmap("probe_mat_temp.txt", mode="w+", dtype="float64", shape=sh[:2])

        # write out the standard deviation for each probe if specified
        if st_dev:
            st_dev_out = open(st_dev_file, "wb")
            st_dev_out.writelines("Probe Node st_dev\node_count")

        for y in range(sh[0]):
            for z in range(sh[1]):
                # mask out unused values, ie where there are less than the maximum number
                # of homolous nodes in a structural region
                probe_mat_temp[y, z] = np.mean(np.ma.array(probe_mat[y, z],
                                                           mask=probe_mat[y, z] == 0.))
                if st_dev:
                    std = np.std(np.ma.array(probe_mat[y, z], mask=probe_mat[y, z] == 0.))
                    st_dev_out.writelines(' '.join([str(int(probe_numbers[y])),
                                                    str(self.imaging_brain.G.nodes()[z]),
                                                    "{:2.5f}".format(float(std))])+'\node_count')

        if st_dev:
            st_dev_out.close()

        # reassign the probe matrix and delete temporary memory-mapped file
        probe_mat = probe_mat_temp
        del probe_mat_temp
        remove(temp_mat_name)

        for gene in gene_names:
            if gene_list[gene]:
                x = probe_mat.shape[1] # number of nodes
                y = len(gene_list[gene]) # number of probes

                gene_mat = np.zeros(shape=(x, y), dtype="float64")

                # nb: probe is the position of the probe recorded above
                for node_count, probe in enumerate(gene_list[gene]):
                    gene_mat[:, node_count] = probe_mat[probe, :]

                # collapse values across probes for each gene
                mean_gene = np.mean(np.ma.array(gene_mat, mask=np.isnan(gene_mat)), axis=1)

                out_dict = dict(list(zip([str(v) for v in self.imaging_brain.G.nodes()],
                                         ["{:10.20f}".format(v) for v in mean_gene])))
                out_dict["Gene"] = gene
                writer.writerow(out_dict)

        out.close()
        remove("probe_mat_temp.txt") # delete memory map file

    def ymatrixgroup(self, metric_dict, subj="Control",
                     allen_dir="allen_data", out_file="YmatrixGroup.csv"):
        '''
        Collates metrics in to a matrix for use in partial least squares analysis.
        Note, the metric_dict contains the metric name as a key and filename as
        the value. Takes group level measures
        '''
        out = open(out_file, "wb")
        headers = ["Metric"]
        headers.extend([str(v) for v in self.imaging_brain.G.nodes()])
        writer = csv.DictWriter(out, fieldnames=headers)
        writer.writeheader()

        # iterate through the metrics
        for m in list(metric_dict.keys()):
            file = open(path.join(allen_dir, subj, metric_dict[m]), "r")
            reader = csv.DictReader(file, delimiter=" ")
            line = next(reader)

            # remove non-numeric keys
            for v in list(line.keys()):
                try:
                    int(v)
                except ValueError:
                    del line[v]

            m_dict = {v:line[v] for v in list(line.keys()) if int(v) in self.imaging_brain.G.nodes()}
            file.close()

            # normalise within metric
            mean_metric = np.mean([float(v) for v in list(m_dict.values())])
            sd_metric = np.std([float(v) for v in list(m_dict.values())])

            m_dict = {str(v):(float(m_dict[str(v)])-mean_metric)/sd_metric for v in list(m_dict.keys())}

            m_dict["Metric"] = m

            writer.writerow(m_dict)

    def ymatrixindividuals(self, metric_dict, subj_list,
                           allen_dir="allen_data", out_file="YmatrixInd.csv"):
        '''
        Collates metrics in to a matrix for use in partial least squares analysis.
        Note, the metric_dict contains the metric name as a key and filename as
        the value. Takes metrics for individual subjects defined in the subject list.
        '''
        out = open(out_file, "wb")
        headers = ["Metric", "Subject"]
        headers.extend([str(v) for v in self.imaging_brain.G.nodes()])
        writer = csv.DictWriter(out, fieldnames=headers)
        writer.writeheader()

        # iterate through the metrics
        for m in list(metric_dict.keys()):
            for subj in subj_list:
                file = open(path.join(allen_dir, subj, metric_dict[m]), "r")
                reader = csv.DictReader(file, delimiter=" ")
                line = next(reader)

                # remove non-numeric keys
                for v in list(line.keys()):
                    try:
                        int(v)
                    except ValueError:
                        del line[v]

                m_dict = {v:line[v] for v in list(line.keys()) if int(v) in self.imaging_brain.G.nodes()}
                file.close()

                # normalise within metric
                mean_metric = np.mean([float(v) for v in list(m_dict.values())])
                sd_metric = np.std([float(v) for v in list(m_dict.values())])

                m_dict = {str(v):(float(m_dict[str(v)]) - mean_metric) / sd_metric for v in list(m_dict.keys())}

                m_dict["Metric"] = m
                m_dict["Subject"] = subj

                writer.writerow(m_dict)
