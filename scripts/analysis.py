#!/usr/bin/python3

import sys
import os
import argparse
import re
import numpy as np
import pandas as pd
# import MDAnalysis as mda
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt
import h5py
import pprint
import shutil
import subprocess
import sklearn
import time
import analysis_helper_functions as helper

from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from MDAnalysis.lib.distances import distance_array
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, default=None, help="path to trajectory file")
parser.add_argument("--config", type=str, default="config.ini", help="path to config file")
parser.add_argument("--solvent", type=str, nargs='*', default=["OSMOT"], help="solvent selection rule")
parser.add_argument("--nonsolvent", type=str, nargs='*', default=["MOBIL","FRAME"], help="nonsolvent selection rule")
parser.add_argument("--clstr_eps", type=float, default=1.2, help="max distance for cluster algorithm")
parser.add_argument("--start", type=int, default=0, help="starting time of analysis")
parser.add_argument("--stop", type=int, default=10e10, help="starting time of analysis")
parser.add_argument("--forcenew", action='store_true', help="force new hdf5 file")
parser.add_argument("--lowmem", action='store_true', help="dont save resname, saves memory BIG TIME")
parser.add_argument("--timestats", action='store_true', help="show timer statistics")
# parser.add_argument("--reanalyze", action='store_true', help="reanalze from data.h5 file instead of trajectory")
args = parser.parse_args()


pp = pprint.PrettyPrinter(indent=4, compact=False)
np.set_printoptions(suppress=True)

try:
    trajfile = h5py.File(args.file,'r')
except Exception as e:
    print(e)
    print("unable to open ", args.file)
# else:
#     raise EnvironmentError("unable to find ", args.file)


# creating the storage file object
if args.forcenew and os.path.exists("data.h5"):
    os.remove("data.h5")
try:
    # datafile = pd.HDFStore("data.h5", "a")
    datafile = pd.HDFStore('data.h5', 'a', complevel=9, complib='zlib')
except Exception as e:
    print(e)
    sys.exit()


resname_map = defaultdict(lambda : -1, {
    "MOBIL" : 0,
    "FRAME" : 1,
    "OSMOT" : 2
})

attributes_setup_done = False

t_start = time.perf_counter()
for key in sorted([s for s in trajfile.keys() if s.startswith("snapshot")], key=lambda x: int(re.findall('\d+', x )[0])):

    hdfgroup = trajfile.get(key)
    attributes = hdfgroup.attrs
    actual_time = int(re.findall('\d+', key )[0])
    
    if actual_time < args.start:
        continue
    elif actual_time > args.stop:
        continue
    # print("time ", actual_time)

    dimensions = [
        attributes.get("system.box.x"),
        attributes.get("system.box.y"),
        attributes.get("system.box.z"),
        90,90,90
    ]

    if not attributes_setup_done:
        _attributes = pd.DataFrame(helper.getAttributeDict(args.config, dimensions[:3]))
        datafile["attributes"] = _attributes
        attributes_setup_done = True

    epot_calc = helper.EpotCalculator(attributes)

    t_prep = time.perf_counter()
    
    positions = pd.DataFrame(hdfgroup.get("position")[()], columns=['x','y','z'])
    orientations = pd.DataFrame(hdfgroup.get("orientation")[()], columns=['ux','uy','uz'])
    particledata = pd.concat([positions,orientations], axis=1).astype(np.float32) 



    """
    add particle (residue) names
    """
    particledata["resname"] = pd.Series(hdfgroup.get("type")[()]).astype(np.int16)

    relevant_cond = particledata["resname"] != 2
    relevant_positions = particledata[relevant_cond].filter(["x","y","z"])
    relevant_orientations = particledata[relevant_cond].filter(["ux","uy","uz"])
    particledata["ux"] = np.where(relevant_cond, particledata["ux"], np.nan)
    particledata["uy"] = np.where(relevant_cond, particledata["uy"], np.nan)
    particledata["uz"] = np.where(relevant_cond, particledata["uz"], np.nan)



    """
    scan for clusters
    """
    distances_array = distance_array(relevant_positions.values, relevant_positions.values, box=dimensions)
    dbscan = DBSCAN(min_samples=2, eps=args.clstr_eps, metric="precomputed", n_jobs=-1).fit(distances_array)
    labels = pd.DataFrame(np.append(dbscan.labels_, np.full(np.count_nonzero(~relevant_cond), -2)), columns=['cluster'])
    # add to data and sort for cluster id
    particledata["cluster"] = labels

    unique, counts = np.unique(labels, return_counts=True)
    particledata["clustersize"] = particledata["cluster"].apply( lambda x: counts[np.where(unique == x)][0] if x >= 0 else 1 )
    particledata.loc[particledata["cluster"] == -1, "clustersize"] = 1

    if args.timestats: print(f"prep took     {time.perf_counter()-t_prep:.4f} seconds")



    """
    subcluster identification
    """
    t_sub = time.perf_counter()
    particledata["subcluster"] = -1
    for ID, group in particledata[relevant_cond].groupby(["cluster"], as_index=False):
        subclusters = helper.getSubclusterLabels(ID, group, args.clstr_eps)
        particledata.loc[group.index, "subcluster"] = subclusters

    if args.timestats: print(f"subclstr took {time.perf_counter()-t_sub:.4f} seconds")



    """
    shift subclusters towards largest subcluster
    """
    t_shift = time.perf_counter()
    particledata["shiftx"] = np.where(relevant_cond, particledata["x"], np.nan)
    particledata["shifty"] = np.where(relevant_cond, particledata["y"], np.nan)
    particledata["shiftz"] = np.where(relevant_cond, particledata["z"], np.nan)
    for ID, group in particledata[relevant_cond].groupby("cluster"):
        newx, newy, newz = helper.getShiftedCoordinates(ID, group, args.clstr_eps, dimensions[:3])
        particledata.loc[newx.index, "shiftx"] = newx.values
        particledata.loc[newy.index, "shifty"] = newy.values
        particledata.loc[newz.index, "shiftz"] = newz.values
    
    if args.timestats: print(f"shift took    {time.perf_counter()-t_shift:.4f} seconds")



    """
    get the order of particle in cluster
    """
    t_order = time.perf_counter()
    particledata["order"] = np.nan
    for ID, group in particledata[relevant_cond].groupby("cluster"):
        orders = helper.getOrder(ID, group)
        particledata.loc[group.index, "order"] = orders
    
    if args.timestats: print(f"order took    {time.perf_counter()-t_order:.4f} seconds")



    """
    get the volume per cluster DBSCAN
    """
    t_volume = time.perf_counter()
    particledata["volume"] = np.nan
    for ID, group in particledata[relevant_cond].groupby(["cluster"]):
        volume = helper.getClusterVolume(ID, group, args.clstr_eps, 5)
        particledata.loc[group.index, "volume"] = volume
        if volume / np.cumprod(dimensions[:3])[-1] > 0.9:
            raise Exception(f"volume of cluster {ID} is {volume / np.cumprod(dimensions[:3])[-1]} of box volume")
    
    if args.timestats: print(f"volume took   {time.perf_counter()-t_volume:.4f} seconds")


    
    """
    calculate the potential energy per particle
    """
    t_epot = time.perf_counter()
    epot, chi = epot_calc.get(relevant_positions, relevant_orientations, dimensions, ret="epot+chi")
    particledata["epot"] = np.append(epot, np.full(np.count_nonzero(~relevant_cond), np.nan))
    particledata["chi"] = np.append(chi, np.full(np.count_nonzero(~relevant_cond), np.nan))
    if args.timestats: print(f"epot and chi took     {time.perf_counter()-t_epot:.4f} seconds")



    """
    calculate the curvature of the structure for every particle
    """
    t_curvature = time.perf_counter()
    particledata["curvature"] = np.append(helper.getCurvature(particledata[relevant_cond], dimensions, cutoff=1.3), np.full(np.count_nonzero(~relevant_cond), np.nan))
    if args.timestats: print(f"curvature took        {time.perf_counter()-t_curvature:.4f} seconds")



    if args.lowmem:
        particledata = particledata.drop(columns=["shiftx","shifty","shiftz"])

    t_write = time.perf_counter()
    datafile[f"time{actual_time}"] = particledata
    if args.timestats: print(f"write took    {time.perf_counter()-t_write:.4f} seconds")
    
    t_end = time.perf_counter()
    print(f"time {actual_time} took {t_end-t_start:.4f} seconds")
    t_start = time.perf_counter()

    # print(datafile[f"time{actual_time}"])

trajfile.close()
datafile.close()