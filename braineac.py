# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 21:59:59 2015

@author: tim
"""

import maybrain.brainObjs as mbo

import csv
from os import path,rename
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from glob import glob

import rpy2.robjects as ro

data = open('t3723687_rda/3723687.rda').read()