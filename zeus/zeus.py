import parse
import numpy as np
import math
import pandas as pd

from collections import namedtuple, defaultdict
from sklearn.neighbors import DistanceMetric

def hl_distance(p,q):
    """
    Compute distance between two distribution
    :param p: [mean, std]
    :param q: [mean, std]
    :return:
    """
    mu1, sigma1 = p[0], p[1]
    mu2, sigma2 = q[0], q[1]
    return (1 - math.sqrt( (2*math.sqrt(sigma1)*math.sqrt(sigma2))/(sigma1+sigma2) )*math.exp(-0.25*((mu1-mu2)**2)/(sigma1+sigma2)))

def batch_parser(fileBatch):
    """
    :param fileBatch: file contain paths of splicing events
    :return:
    """
    mainFrame = defaultdict(dict)
    batch = pd.read_csv(fileBatch, sep="\t", header=None).to_dict()
    batch = batch[0]
    for sample, file in batch.items():
        schizo = pd.read_csv(file, sep="\t").set_index("#ID").transpose().to_dict()
        map(lambda x: mainFrame[x].update({sample:[schizo[x]["PSI_bootstrap"], schizo[x]["PSI_bootstrap_std"]]}), schizo.keys())
    return mainFrame

def coordGet(distMat, size_1):
    """
    get x,y of coordinate of certain gene
    :param distMat:
    :param size_1:
    :param size_2:
    :return:
    """
    block1  = reduce(list.__add__,[list(x)[:size_1] for x in [y for y in list(distMat)[:size_1]]], [])
    block2  = reduce(list.__add__,[list(x)[size_1:] for x in [y for y in list(distMat)[size_1:]]], [])
    return (np.mean(block1), np.mean(block2))

def main(fcases, fcontrol):
    """
    main instance for different splicing
    :param fcases:
    :param fcontrol:
    :return:
    """
    mainFrame1 = batch_parser(fcases)
    mainFrame2 = batch_parser(fcontrol)
    def _pairwiseCompute(id):
        """
        compute pairwise distances of all distributions, which were represented as means and stds
        :param disMat:
        :return:
        """
        mat1, mat2 = np.array(mainFrame1[id].values()), np.array(mainFrame2[id].values())
        distMat = np.concatenate(mat1, mat2)
        dist = DistanceMetric("pyfunc", func=hl_distance)
        p_dist = dist.pairwise(distMat)

        return coordGet(p_dist, len(mat1))
    map(_pairwiseCompute, mainFrame1.keys())
    return zip(mainFrame1.keys(), _pairwiseCompute)