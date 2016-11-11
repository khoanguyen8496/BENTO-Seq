import numpy as np
import math
import pandas as pd
import json
import logger as Logger

from sklearn.neighbors import DistanceMetric
from multiprocessing import Pool
from contextlib import closing
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from time import time

diagonal_d = lambda x, y: abs(x-y)/math.sqrt(2)
oLogger = Logger.LogInstance("splicing_log.txt")

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

def _read_file(file):
    ret = pd.read_csv(file, sep="\t").set_index("#ID").transpose().to_dict()
    return ret

def batch_parser(fileBatch, oLogger=oLogger, threads=20):
    """
    :param fileBatch: file contain paths of splicing result for all samples
    :return:
    """
    start=time()
    oLogger.info("Parsing splicing results file..")
    batch = pd.read_csv(fileBatch, sep="\t", header=None).to_dict()[0]
    with closing(Pool(processes=threads)) as pool:
        dummies = pool.map(_read_file, batch.values())
        pool.close()
        pool.terminate()
    splicing_batches = dict(zip(batch.keys(), dummies))
    dummy_keys = dummies[0].keys()
    arrPSIs = reduce(list.__add__, [ [(batch[key]["PSI_bootstrap"],batch[key]["PSI_bootstrap_std"])]\
                                                  for batch in (splicing_batches.values()) for key in dummy_keys])
    print time()-start
    return dict(zip(dummy_keys, arrPSIs))

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

def main(fcontrol, fcases, control, cases, base_data="/raid/ptdtan/base_data.json", threads=20, oLogger=oLogger):
    """
    main instance for different splicing
    :param fcases:
    :param fcontrol:
    :return:
    """
    oLogger.info("Parsing splicing result file: %s" %(fcases))
    mainFrame1 = batch_parser(fcases)
    oLogger.info("Parsing splicing result file: %s" %(fcontrol))
    mainFrame2 = batch_parser(fcontrol)

    keys = mainFrame1.keys()
    base_Data = json.load(open(base_data))

    def _pairwiseCompute(id):
        """
        compute pairwise distances of all distributions, which were represented as means and stds
        :param disMat:
        :return:
        """
        #oLogger.info("Computing exons ID: %s" %id)
        mat1, mat2 = mainFrame1[id].values(), mainFrame2[id].values()
        distMat = mat1 + mat2
        dist = DistanceMetric.get_metric("pyfunc", func=hl_distance)
        p_dist = dist.pairwise(distMat)
        return coordGet(p_dist, len(mat1))

    with closing(Pool(processes=threads)) as pool:
        coords = pool.map(_pairwiseCompute, mainFrame1.keys())
        pool.close()
        pool.terminate()

    d = np.array(map(diagonal_d, coords)) #distances
    z = (d-np.mean(d))/np.std(d) #z-scores
    p_vals = stats.norm.sf(abs(z))
    p_adjs = multipletests(p_vals, method="fdr_tsbh") #False discovery rate, two-stage Benjamini/Hochberg

    def updateJSON(i):
        """
        update JSON data
        :param idx:
        :return:
        """
        uKey = keys[i]
        base_Data[uKey].update({
            "PSI_mean":{
                cases:np.mean([x for x,y in mainFrame1[uKey].values()]),
                control:np.mean([x for x,y in mainFrame2[uKey].values()]),
            },
            "distance": d[i],
            "z-score": z[i],
            "p_val": p_vals[i],
            "p_adjs": p_adjs[1][i],
            "x": coords[i][0],
            "y": coords[i][1]
        })
    map(updateJSON, range(len(keys)))
    return base_Data

def write_data(base_Data, fileOut="splicing_result.tsv"):
    """
    fout = json.dump
    :param base_Data:
    :return:
    """
    fout = open(fileOut, "w")
    fout.write("gene\ttriplet\tstrand\tC1_start\tC1_end\tA_start\tA_end\tC2_start\tC2_end\tPSI_mean %s\tPSI_mean %s\tdistance\tz-score\tp_val\tp_adjs\tx\ty\n" %(cases, control))
    for key, value in base_Data.items():
        fout.write("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" %(value["gene"], value["triplet"], value["strand"], value["C1_start"], value["C1_end"], value["A_start"], value["A_end"], value["C2_start"], value["C2_end"], value["PSI_mean"][cases], value["PSI_mean"][control], value["distance"], value["z-score"], value["p_val"], value["p_adjs"], value["x"], value["y"]))



