import numpy as np
import math
import pandas as pd
import json
import logger as Logger
import os
import pickle
import argparse
import sys, traceback

from sklearn.neighbors import DistanceMetric
from multiprocessing import Pool
from contextlib import closing
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from collections import defaultdict

diagonal_d = lambda x: abs(x[0]-x[1])/math.sqrt(2)
oLogger = Logger.LogInstance("splicing_log.txt")

def index(falls, serializedObject="/raid/ptdtan/splicing_db.bin", threads=20):
    """
    index splicing database for fast computing
    :param fall:
    :param serializedObject:
    :return:
    """
    Oall = batch_parser(open(falls).read().strip().split("\n"), threads=threads)
    pickle.dump(Oall, open(serializedObject, "wb"))
    pass

def debugException(exc_type, exc_value, exc_traceback):
    formatted_lines = traceback.format_exc().splitlines()
    oLogger.debug(formatted_lines[0])
    oLogger.debug(formatted_lines[-1])
    oLogger.debug(repr(traceback.format_exception(exc_type, exc_value,
                                          exc_traceback)))
    oLogger.debug(repr(traceback.extract_tb(exc_traceback)))
    oLogger.debug(repr(traceback.format_tb(exc_traceback)))
    pass

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

def batch_parser(batches, oLogger=oLogger, threads=20):
    """
    :param fileBatch: file contain paths of splicing result for all samples
    :return:
    """
    oLogger.info("Parsing splicing results files..")
    with closing(Pool(processes=threads)) as pool:
        dummies = pool.map(_read_file, batches)
        pool.close()
        pool.terminate()
    oLogger.info("Done reading files. Reducing")
    splicing_batches = dict(zip([file.split("/")[-1] for file in batches], dummies))
    dummy_keys = dummies[0].keys()
    mainFrame = defaultdict(list)
    for key in dummy_keys:
        mainFrame[key] = dict(zip(splicing_batches.keys(), [(batch[key]["PSI_bootstrap"],batch[key]["PSI_bootstrap_std"]) \
                                                            for batch in splicing_batches.values()]))
    oLogger.info("Done reducing..")
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

def _pairwiseCompute(arg):
    """
    compute pairwise distances of all distributions, which were represented as means and stds
    :param disMat:
    :return
    """
    mat1, mat2, idx = arg
    oLogger.info("Computing exons ID: %s" %idx)
    distMat = mat1 + mat2
    dist = DistanceMetric.get_metric("pyfunc", func=hl_distance)
    p_dist = dist.pairwise(distMat)
    return coordGet(p_dist, len(mat1))

def main(fcases, fcontrols,serializedObject=None,
         base_data="/raid/ptdtan/base_data.json",
         threads=20,
         oLogger=oLogger,
         fileOut="splicing_result.tsv"):
    """
    main instance for different splicing
    :param fcases:
    :param fcontrol:
    :return:
    """
    def _parser(fcases, fcontrols):
        """
        parse list of file or file of files
        :param fcontrols:
        :param fcases:
        :return:
        """
        if not any([os.path.exists(f) for f in [fcontrols, fcases]]):
            fcases, fcontrols = fcontrols.strip().split(","), fcases.strip().split(",")
            return fcases, fcontrols
        fcases, fcontrols = open(fcases).read().strip().split("\n"), open(fcontrols).read().strip().split("\n")
        return fcases, fcontrols
    def _isSerialied(serializedObject, cases, controls):
        """
        load serialized splicing results if presented
        :param serializedObject:
        :return:
        """
        try:
            if not serializedObject:
                oLogger.info("Parsing splicing result file: %s" % (fcases))
                mainFrame1 = batch_parser(cases)
                oLogger.info("Parsing splicing result file: %s" % (fcontrols))
                mainFrame2 = batch_parser(controls)
            else:
                mainFrame = pickle.load(open(serializedObject))
                k_cases, k_controls = [x.strip().split("/")[-1] for x in cases], [x.strip().split("/")[-1] for x in controls]
                mainFrame1 = {key:[mainFrame[key][sample] for sample in k_cases] for key in mainFrame.keys()}
                mainFrame2 = {key: [mainFrame[key][sample] for sample in k_controls] for key in mainFrame.keys()}
        except Exception as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            debugException(exc_type, exc_obj, exc_tb)
        return mainFrame1, mainFrame2
    try:
        fcases, fcontrols = _parser(fcases, fcontrols)
        mainFrame1, mainFrame2 = _isSerialied(serializedObject, fcases, fcontrols)
    except Exception as err:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        debugException(exc_type, exc_obj, exc_tb)
    keys = mainFrame1.keys()
    base_Data = json.load(open(base_data))
    with closing(Pool(processes=threads)) as pool:
        coords = pool.map(_pairwiseCompute, [(mainFrame1[idx], mainFrame2[idx], idx) \
                                                for idx in mainFrame1.keys()])
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
                1:np.mean([x for x,y in mainFrame1[uKey]]),
                2:np.mean([x for x,y in mainFrame2[uKey]]),
            },
            "distance": d[i],
            "z-score": z[i],
            "p_val": p_vals[i],
            "p_adjs": p_adjs[1][i],
            "x": coords[i][0],
            "y": coords[i][1]
        })
    map(updateJSON, range(len(keys)))
    write_data(base_Data, fileOut=fileOut)
    pass
def write_data(base_Data, fileOut="splicing_result.tsv"):
    """
    fout = json.dump
    :param base_Data:
    :return:
    """
    fout = open(fileOut, "w")
    fout.write("gene\ttriplet\tstrand\tC1_start\tC1_end\tA_start\tA_end\tC2_start\tC2_end\tPSI_mean_1\tPSI_mean_2\tdistance\tz-score\tp_val\tp_adjs\tx\ty\n")
    for key, value in base_Data.items():
        fout.write("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%e\t%e\t%.5f\t%.5f\n" %(value["gene"], value["triplet"], value["strand"], value["C1_start"], value["C1_end"], value["A_start"], value["A_end"], value["C2_start"], value["C2_end"], value["PSI_mean"][1], value["PSI_mean"][2], value["distance"], value["z-score"], value["p_val"], value["p_adjs"], value["x"], value["y"]))
if __name__ == "__main__":
    main_parser = argparse.ArgumentParser(description='Zeus - Different exons splicing comparison')
    subparsers = main_parser.add_subparsers(help='Zeus Utilities', dest="cmd")
    main_parser.add_argument("--ThreadsN", "-t",
                             type=int,
                             help="Number of threads",
                             default=20)

    index_parser = subparsers.add_parser('index', help='Index all splicing results')
    index_parser.add_argument("--files", "-i",
                              type=str,
                              required=True,
                              help="file contained path of all splicing files in workspace")
    index_parser.add_argument("--outfile", "-o",
                              type=str,
                              help="Output binary file of all index splicing result",
                              default="./splicing_result.bin")
    index_parser.set_defaults(which='index')

    compute_parser = subparsers.add_parser('compute', help='Compute different splicing events')
    compute_parser.add_argument("--outfile", "-o",
                                type=str,
                                help="Output binary file of all index splicing result",
                                default="./splicing_result_2.bin")
    compute_parser.add_argument('--controls', '-co', type=str,
                                help='Input list of file, separated by "," or file of files controls',
                                required=True)
    compute_parser.add_argument('--cases', '-ca', type=str,
                                help='Input list of file, separated by "," or file of files cases',
                                required=True)
    compute_parser.add_argument('--base-data', '-b', type=str, help='Annotation DB as JSON', required=True)
    compute_parser.add_argument('--serialized-bin', '-sb', type=str, help='serialized bin', required=True)
    compute_parser.set_defaults(which='compute')

    args = main_parser.parse_args()
    if args.cmd =="index":
        exit(index(args.files, serializedObject=args.outfile, threads=args.ThreadsN))
    if args.cmd == "compute":
        exit(main(fcases=args.cases, fcontrols=args.controls, serializedObject=args.serialized_bin, threads=args.ThreadsN,
                  base_data=args.base_data, fileOut=args.outfile))
