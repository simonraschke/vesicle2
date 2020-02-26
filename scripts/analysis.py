#!/usr/bin/python3

import sys
import os
import argparse
import re
import numpy as np
import pandas as pd
import h5py
import time
import analysis_helper_functions as helper

from sklearn.cluster import DBSCAN
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

ge_restriction_map = defaultdict(lambda : -1, {
    "inplace" : 0,
    "structure" : 1
})


LAST_DF = pd.DataFrame()
origin_positions = pd.DataFrame()
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

    dimensions = np.array([
        attributes.get("system.box.x"),
        attributes.get("system.box.y"),
        attributes.get("system.box.z"),
        90,90,90
    ])


    t_prep = time.perf_counter()
    
    positions = pd.DataFrame(hdfgroup.get("position")[()], columns=['x','y','z'])
    orientations = pd.DataFrame(hdfgroup.get("orientation")[()], columns=['ux','uy','uz'])
    particledata = pd.concat([positions,orientations], axis=1).astype(np.float32) 


    if not attributes_setup_done:
        origin_positions = pd.DataFrame(hdfgroup.get("position")[()], columns=['x','y','z'])
        LAST_DF = particledata
        LAST_DF = LAST_DF.assign(xnopbc=LAST_DF["x"], ynopbc=LAST_DF["y"], znopbc=LAST_DF["z"])
        _attributes = pd.DataFrame(helper.getAttributeDict(args.config, dimensions[:3]))
        fga_mode = _attributes['fga_mode'].values[0]
        simulation_mode = _attributes['simulation_mode'].values[0]
        if fga_mode == "plane":
            plane_edge = _attributes['plane_edge'].values[0]
            guiding_elements = _attributes['guiding_elements_each'].values[0]
            guiding_elements_per_dim = int(np.sqrt(guiding_elements))
            print(plane_edge, guiding_elements_per_dim, dimensions[:3]/2)
            domains = helper.generateDomains(plane_edge, guiding_elements_per_dim, dimensions[:3]/2)
            domains_internal = helper.generateDomainsInternal(plane_edge, guiding_elements_per_dim, dimensions[:3]/2, attributes.get('system.ljsigma') )
            # for d in domains_internal:
            #     print(d)
            _attributes["domain_volume"] = domains[0].volume if guiding_elements > 0 else 0.0
        elif fga_mode == "sphere":
            guiding_elements = _attributes['guiding_elements_each'].values[0]
        datafile["attributes"] = _attributes
        attributes_setup_done = True

    epot_calc = helper.EpotCalculator(attributes)


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
    particledata["neighbours"] = np.count_nonzero(distances_array <= attributes["system.ljsigma"]*1.5, axis=1).astype(np.int16) - 1
    particledata["cluster"] = labels

    unique, counts = np.unique(labels, return_counts=True)
    particledata["clustersize"] = particledata["cluster"].apply( lambda x: counts[np.where(unique == x)][0] if x >= 0 else 1 )
    particledata.loc[particledata["cluster"] == -1, "clustersize"] = 1

    if args.timestats: print(f"\nprep took             {time.perf_counter()-t_prep:.4f} seconds")



    """
    subcluster identification
    """
    t_sub = time.perf_counter()
    particledata.loc[relevant_cond.index, "subcluster"] = particledata[relevant_cond].groupby("cluster", group_keys=False).apply(lambda g: helper.getSubclusterLabels(g, args.clstr_eps))["subcluster"]

    if args.timestats: print(f"subclstr took         {time.perf_counter()-t_sub:.4f} seconds")



    """
    shift subclusters towards largest subcluster
    """
    t_shift = time.perf_counter()
    particledata = particledata.assign(shiftx=np.nan, shifty=np.nan, shiftz=np.nan)
    particledata.loc[relevant_cond.index, ["shiftx","shifty","shiftz"]] = particledata[relevant_cond].groupby("cluster", group_keys=False).apply(lambda g: helper.getShiftedCoordinates(g, args.clstr_eps, dimensions[:3]))[["shiftx","shifty","shiftz"]]
    
    if args.timestats: print(f"shift took            {time.perf_counter()-t_shift:.4f} seconds")



    """
    special analysis of fga structure
    """
    if simulation_mode == "FGA" and guiding_elements > 0:
        t_fga_plane = time.perf_counter()
        particledata["in_structure"] = False
        particledata["in_structure_cluster"] = False
        in_structure, structure_volume = helper.isParticleInStructure(particledata[relevant_cond], attributes, dimensions, fga_mode)
        particledata.loc[relevant_cond, "in_structure"] = in_structure
        particledata.loc[in_structure, "structure_volume"] = structure_volume
        particledata.loc[relevant_cond, "in_structure_cluster"] = helper.isParticleInStructureCluster(particledata[relevant_cond], attributes)
        particledata.loc[relevant_cond, "in_structure_env"] = helper.isParticleInStructureEnvironment(particledata[relevant_cond], attributes, dimensions, fga_mode)
        
        if args.timestats: print(f"plane fga case took   {time.perf_counter()-t_fga_plane:.4f} seconds")
        # print(particledata)


    """
    get the order of particle in cluster
    """
    t_order = time.perf_counter()
    particledata.loc[relevant_cond.index, "order"] = particledata[relevant_cond].groupby("cluster", group_keys=False).apply(lambda g: helper.getOrder(g, attributes))["order"]
    
    if args.timestats: print(f"order took            {time.perf_counter()-t_order:.4f} seconds")



    """
    get the volume per cluster DBSCAN
    """
    t_volume = time.perf_counter()
    particledata.loc[relevant_cond.index, "volume"] = particledata[relevant_cond].groupby("cluster", group_keys=False).apply(lambda g: helper.getClusterVolume(g, args.clstr_eps, 5))["volume"]
    
    if args.timestats: print(f"volume took           {time.perf_counter()-t_volume:.4f} seconds")


    
    """
    calculate the potential energy per particle
    """
    t_epot = time.perf_counter()
    epot, chi = epot_calc.get(relevant_positions, relevant_orientations, dimensions, ret="epot+chi", cutoff=attributes.get("system.ljsigma")*3, distances_array=distances_array)
    neigbours_chi = epot_calc.get(relevant_positions, relevant_orientations, dimensions, ret="chi", cutoff=attributes.get("system.ljsigma")*1.5, distances_array=distances_array)
    particledata.loc[relevant_cond.index, "epot"] = epot
    particledata.loc[relevant_cond.index, "chi"] = chi
    particledata.loc[relevant_cond.index, "neighbours_chi"] = neigbours_chi

    if args.timestats: print(f"epot and chi took     {time.perf_counter()-t_epot:.4f} seconds")



    """
    calculate the curvature of the structure for every particle
    """
    t_curvature = time.perf_counter()
    particledata["curvature"] = np.append(helper.getCurvature(particledata[relevant_cond], dimensions, cutoff=attributes.get("system.ljsigma")*1.5), np.full(np.count_nonzero(~relevant_cond), np.nan))
    particledata.loc[particledata["order"].isna(), "curvature"] = np.nan

    if args.timestats: print(f"curvature took        {time.perf_counter()-t_curvature:.4f} seconds")



    """
    calculate particle structure domain and the domain volume
    """
    if fga_mode == "plane" and guiding_elements > 0:
        t_population = time.perf_counter()
        particledata["structure_domain"] = np.int8(-1)
        particledata.loc[particledata["in_structure_env"], "structure_domain"] = helper.getStructureDomainID(particledata[particledata["in_structure_env"]], domains)
        particledata["structure_domain_internal"] = np.int8(-1)
        particledata.loc[particledata["in_structure_env"], "structure_domain_internal"] = helper.getStructureDomainIDInternal(particledata[particledata["in_structure_env"]], domains_internal)

        if args.timestats: print(f"structure domain took {time.perf_counter()-t_population:.4f} seconds")



    """
    Mean Square Displacement
    """
    t_msd = time.perf_counter()
    particledata = particledata.assign(xnopbc=np.nan, ynopbc=np.nan, znopbc=np.nan)
    particledata[["xnopbc","ynopbc","znopbc"]] = helper.correctForPBC(LAST_DF, particledata, dimensions)
    particledata["MSD"] = helper.getMSD(origin_positions, particledata[["xnopbc","ynopbc","znopbc"]].values, dimensions)
    if args.timestats: print(f"MSD took              {time.perf_counter()-t_msd:.4f} seconds")
    
    

    """
    calculate the surface tension for every particle
    """
    # t_tension = time.perf_counter()
    # particledata["tension"] = np.append(helper.getSurfaceTension(particledata[relevant_cond], dimensions, epot_calc, cutoff=3.0), np.full(np.count_nonzero(~relevant_cond), np.nan))
    # if args.timestats: print(f"surface tension took  {time.perf_counter()-t_tension:.4f} seconds")



    if args.lowmem:
        particledata = particledata.drop(columns=["shiftx","shifty","shiftz"])

    t_write = time.perf_counter()
    datafile[f"time{actual_time}"] = particledata
    if args.timestats: print(f"write took            {time.perf_counter()-t_write:.4f} seconds")
    
    t_end = time.perf_counter()
    print(f"time {actual_time} took {t_end-t_start:.4f} seconds")
    t_start = time.perf_counter()

    LAST_DF = particledata
    # with pd.option_context('display.max_rows', None):  # more options can be specified also
    #     print(particledata)\
    # print(particledata.head(30))
    # np.set_printoptions(precision=2, linewidth=200, floatmode="fixed")
    # vals.append(particledata["MSD"].mean())
    # sys.exit()

    # print(datafile[f"time{actual_time}"])
trajfile.close()
datafile.close()