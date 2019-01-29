import sys
import numpy as np
import pandas as pd
import scipy
import MDAnalysis as mda 

from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.lib.distances import distance_array
from sklearn.cluster import DBSCAN
from scipy import ndimage
from sklearn.preprocessing import normalize



# iterate a file linewise
# check line for keyword
# if keyword is found return value after seperator
def fileValueFromKeyword(filepath, keyword, seperator='='):
    found_counter = 0
    assert(filepath)
    value = [np.NaN]
    try:
        with open(filepath) as FILE:
            for line in FILE:
                if keyword in line:
                    found_counter += 1
                    value = [float(line.split(seperator,1)[1].rstrip())]
                    break
    except:
        pass
    if found_counter == 0:
        print(Warning("unable to find keyword "+keyword+" in file "+filepath+"\n"))
        # raise Warning("unable to find keyword "+keyword+" in file "+filepath+"\n")
    # if found_counter >= 2:
        # raise Exception("multiple keywords "+keyword+" in file "+filepath+"\n")
    return value



def getAttributeDict(filename, dimensions):
    return {
        'mobile': fileValueFromKeyword(filename, 'mobile', '='),
        'density': fileValueFromKeyword(filename, 'density', '='),
        'boxx': dimensions[0],
        'boxy': dimensions[1],
        'boxz': dimensions[2],
        'temperature': fileValueFromKeyword(filename, 'temperature', '='),
        'time_max': fileValueFromKeyword(filename, 'time_max', '='),
        'kappa': fileValueFromKeyword(filename, 'kappa', '='),
        'gamma': fileValueFromKeyword(filename, 'gamma', '='),
        'guiding_elements_each': fileValueFromKeyword(filename, 'guiding_elements_each', '='),
        'frame_guides_grid_edge': fileValueFromKeyword(filename, 'frame_guides_grid_edge', '='),
        'plane_edge': fileValueFromKeyword(filename, 'plane_edge', '='),
        'osmotic_density_inside': fileValueFromKeyword(filename, 'osmotic_density_inside', '='),
        'sw_position_min': fileValueFromKeyword(filename, 'sw_position_min', '='),
        'sw_position_max': fileValueFromKeyword(filename, 'sw_position_max', '='),
        'sw_position_target': fileValueFromKeyword(filename, 'sw_position_target', '='),
        'sw_orientation_min': fileValueFromKeyword(filename, 'sw_orientation_min', '='),
        'sw_orientation_max': fileValueFromKeyword(filename, 'sw_orientation_max', '='),
        'sw_orientation_target': fileValueFromKeyword(filename, 'sw_orientation_target', '='),
        'cell_min_edge': fileValueFromKeyword(filename, 'cell_min_edge', '='),
        'max_cells_dim': fileValueFromKeyword(filename, 'max_cells_dim', '='),
        'skip': fileValueFromKeyword(filename, 'skip', '=')
    }



def getSubclusterLabels(ID, group, eps):
    #skip if noise
    if ID == -1:
        return np.zeros( len(group.index), dtype=int )
    else:
        # arange a DBSCAN without PBC to get subclusters
        coms_subcluster = pd.concat([group['x'], group['y'], group['z']], axis=1)
        distances_array_subcluster = distance_array(coms_subcluster.values, coms_subcluster.values, box=None)
        subclusterd = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=-1).fit(distances_array_subcluster)
        return subclusterd.labels_



def getShiftedCoordinates(ID, group, eps, dimensions, alg="dbscan"):
    # get largest subcluster
    if alg == "dbscan":
        unique, counts = np.unique(group["subcluster"], return_counts=True)
    elif alg == "hdbscan":
        unique, counts = np.unique(group["h_subcluster"], return_counts=True)
    if len(unique) == 1:
        return group["x"], group["y"], group["z"]
    max_subclusterID = unique[counts == np.max(counts)][0]
    # calculate shifts per subcluster
    if alg == "dbscan":
        centers = group.groupby("subcluster")['x','y','z'].mean()
    elif alg == "hdbscan":
        centers = group.groupby("h_subcluster")['x','y','z'].mean()
    shifts = np.round(( -centers + centers.loc[max_subclusterID] )/dimensions[:3]).astype(int)
    shifts *= dimensions[:3]

    coms_subcluster = pd.concat([group['x'], group['y'], group['z']], axis=1)
    distances_array_subcluster = distance_array(coms_subcluster.values, coms_subcluster.values, box=None)
    dbscan_subcluster = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=-1).fit(distances_array_subcluster)
    np.set_printoptions(threshold=np.nan, linewidth=np.nan, precision=1)
    # calculate new coordinates based on shift
    if alg == "dbscan":
        newx = np.add(group["x"], shifts.loc[group["subcluster"]]["x"])
        newy = np.add(group["y"], shifts.loc[group["subcluster"]]["y"])
        newz = np.add(group["z"], shifts.loc[group["subcluster"]]["z"])
    elif alg == "hdbscan":
        newx = np.add(group["x"], shifts.loc[group["h_subcluster"]]["x"])
        newy = np.add(group["y"], shifts.loc[group["h_subcluster"]]["y"])
        newz = np.add(group["z"], shifts.loc[group["h_subcluster"]]["z"])
    # assign to main data
    return newx, newy, newz



def getClusterVolume(ID, group, eps, pps):
    # volume = 0
    if ID == -1:
        # single particles as spheres
        return np.pi * (eps**3) * 4/3
    else:
        # generate meshgrid around cluster particles plus threshold
        x_vector = np.arange(np.min(group["shiftx"])-eps, np.max(group["shiftx"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)
        y_vector = np.arange(np.min(group["shifty"])-eps, np.max(group["shifty"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)
        z_vector = np.arange(np.min(group["shiftz"])-eps, np.max(group["shiftz"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)

        # make it recursive if meshgrid becomes too large
        NNN = 80
        if len(x_vector) > NNN or len(y_vector) > NNN or len(z_vector) > NNN:
            if pps-1 < 1:
                return getClusterVolume(ID, group, eps, pps/2)
            else:
                return getClusterVolume(ID, group, eps, pps-1)

        xx, yy, zz = np.meshgrid(x_vector, y_vector, z_vector)
        # stack them together as array of 3D points
        meshgrid = np.stack((xx.ravel(), yy.ravel(), zz.ravel()), axis=1)
        #calculate the distance array with centres of masses of particles
        coms_cluster = pd.concat([group['shiftx'], group['shifty'], group['shiftz']], axis=1)
        getClusterVolume.distances_array_volume = distance_array(meshgrid, coms_cluster.values, box=None).astype(np.float32)
        # check if any point in distance array row is close enough, then reshape to meshgrid
        # result is a binary meshgrid with 1 for the cluster shell region
        isclose = np.where(getClusterVolume.distances_array_volume <= eps, True, False).any(axis=1).reshape(xx.shape[0], yy.shape[1], zz.shape[2])
        # fill hole inside the shell region
        isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
        # calc volum from all points inside cluster
        return ((1.0/pps)**3)*np.count_nonzero(isclose)


        # z,x,y = isclose.nonzero()
        # ax.scatter(x+,y,z, s=2)



def getOrder(ID, group, alg="dbscan"):
    if ID == -1:
        if alg == "dbscan":
            return group["order"].values
        elif alg == "hdbscan":
            return group["h_order"].values
    else:
        if alg == "dbscan":
            shifted_coms = pd.concat([group['shiftx'], group['shifty'], group['shiftz']], axis=1)
        elif alg == "hdbscan":
            shifted_coms = group.filter(['h_shiftx','h_shifty','h_shiftz'])
        normalized_orientations = pd.concat([group['ux'], group['uy'], group['uz']], axis=1)
        return (normalized_orientations.values*normalize(shifted_coms.sub(shifted_coms.mean().values))).sum(axis=1)



def getPairs(distances_array, max_cutoff, same=False):
    valid = distances_array < max_cutoff
    np.fill_diagonal(valid, same)
    pairs = np.column_stack(np.where(valid))
    pairs = np.sort(pairs, axis=1)
    unique, index = np.unique(pairs, axis=0, return_index=True)
    return pairs[index]



def getNormedPairDistanceVectors(coms, pairs, dimensions):
    dist_vecs = np.subtract(coms.loc[pairs[:,1]].values, coms.loc[pairs[:,0]].values)
    dist_vecs = np.subtract(dist_vecs, np.multiply(dimensions[:3], np.round(dist_vecs/dimensions[:3])))
    return pd.DataFrame((dist_vecs.T / np.linalg.norm(dist_vecs, axis=1)).T, columns=["x","y","z"]), np.linalg.norm(dist_vecs, axis=1)



class EpotCalculator(object):
    def __init__(self, attributes):
        self.sigma = attributes["system.ljsigma"]
        self.epsilon = attributes["system.ljepsilon"]
        self.kappa = attributes["system.kappa"]
        self.gamma = attributes["system.gamma"]
        self.a = 1.0 + self.kappa*np.sin(self.gamma*np.pi/180)
        self.b = 1.0 - self.kappa*np.sin(self.gamma*np.pi/180)
        self.c = np.linalg.norm(np.array([self.a,0,0]) + np.array([self.a-1.0, self.kappa*np.cos(self.gamma*np.pi/180),0]))



    def getChi(self, coms, orientations, dimensions, cutoff=3):
        try:
            assert(len(coms) == len(orientations))
            assert(len(dimensions) == 6)
            assert(len(coms) == len(df))
        except:
            print(f"{len(coms)} coms")
            print(f"{len(orientations)} orientations")
            print(f"{len(dimensions)} dimensions")
            return

        distances_array = distance_array(coms.values, coms.values, box=dimensions)
        pairs = getPairs(distances_array, 3)
        normed_dist_vecs, dist_norms = getNormedPairDistanceVectors(coms, pairs, dimensions)
        res1_u = np.multiply(np.take(orientations.values, pairs[:,0], axis=0), self.kappa/2)
        res2_u = np.multiply(np.take(orientations.values, pairs[:,1], axis=0), self.kappa/2)
        chi =  np.power((np.linalg.norm(-res1_u + normed_dist_vecs + res2_u, axis=1) - self.a), 2)
        chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs - res2_u, axis=1) - self.b), 2)
        chi += np.power((np.linalg.norm(-res1_u + normed_dist_vecs - res2_u, axis=1) - self.c), 2)
        chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs + res2_u, axis=1) - self.c), 2)
        epot = 4.0 * self.epsilon * ( np.power(self.sigma/dist_norms, 12) - (1.0 - chi)*np.power(self.sigma/dist_norms, 6) )
        epot_array = np.zeros_like(distances_array)
        pairs_t = pairs.T
        epot_array[tuple(pairs_t)] = epot
        epot_array[tuple([pairs_t[1], pairs_t[0]])] = epot
        return np.sum(epot_array, axis=1)


    
    def get(self, coms, orientations, dimensions, cutoff=3, ret="epot", distances_array=None):
        try:
            assert(len(coms) == len(orientations))
            assert(len(dimensions) == 6)
        except:
            print(f"{len(coms)} coms")
            print(f"{len(orientations)} orientations")
            print(f"{len(dimensions)} dimensions")
            return

        if isinstance(distances_array, type(None)):
            distances_array = distance_array(coms.values, coms.values, box=dimensions)

        if ret == "chi":
            pairs = getPairs(distances_array, 1.3*self.sigma)
            normed_dist_vecs, dist_norms = getNormedPairDistanceVectors(coms, pairs, dimensions)
            res1_u = np.multiply(np.take(orientations.values, pairs[:,0], axis=0), self.kappa/2)
            res2_u = np.multiply(np.take(orientations.values, pairs[:,1], axis=0), self.kappa/2)
            chi =  np.power((np.linalg.norm(-res1_u + normed_dist_vecs + res2_u, axis=1) - self.a), 2)
            chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs - res2_u, axis=1) - self.b), 2)
            chi += np.power((np.linalg.norm(-res1_u + normed_dist_vecs - res2_u, axis=1) - self.c), 2)
            chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs + res2_u, axis=1) - self.c), 2)
            chi_array = np.zeros_like(distances_array)
            pairs_t = pairs.T
            chi_array[tuple(pairs_t)] = chi
            chi_array[tuple([pairs_t[1], pairs_t[0]])] = chi
            return np.sum(chi_array, axis=1)
        elif ret == "epot":
            pairs = getPairs(distances_array, cutoff)
            normed_dist_vecs, dist_norms = getNormedPairDistanceVectors(coms, pairs, dimensions)
            res1_u = np.multiply(np.take(orientations.values, pairs[:,0], axis=0), self.kappa/2)
            res2_u = np.multiply(np.take(orientations.values, pairs[:,1], axis=0), self.kappa/2)
            chi =  np.power((np.linalg.norm(-res1_u + normed_dist_vecs + res2_u, axis=1) - self.a), 2)
            chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs - res2_u, axis=1) - self.b), 2)
            chi += np.power((np.linalg.norm(-res1_u + normed_dist_vecs - res2_u, axis=1) - self.c), 2)
            chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs + res2_u, axis=1) - self.c), 2)
            epot = 4.0 * self.epsilon * ( np.power(self.sigma/dist_norms, 12) - (1.0 - chi)*np.power(self.sigma/dist_norms, 6) )
            epot_array = np.zeros_like(distances_array)
            pairs_t = pairs.T
            epot_array[tuple(pairs_t)] = epot
            epot_array[tuple([pairs_t[1], pairs_t[0]])] = epot
            return np.sum(epot_array, axis=1)

        elif ret == "epot+chi":
            return self.get(coms, orientations, dimensions, cutoff=3, ret="epot", distances_array=distances_array), \
                   self.get(coms, orientations, dimensions, ret="chi",  distances_array=distances_array)

        elif ret == "chi+epot":
            return self.get(coms, orientations, dimensions, ret="chi",  distances_array=distances_array), \
                   self.get(coms, orientations, dimensions, cutoff=3, ret="epot", distances_array=distances_array)



def getCurvature(particledata, dimensions, cutoff=13):
    # np.set_printoptions(threshold=np.nan, linewidth=np.nan, precision=4)
    # b = np.array([1,0,0])
    # a = np.array([-0.2,1,0])
    # a1 = np.dot(a, b) / np.linalg.norm(b)
    # print(a1)
    coms = particledata.filter(['shiftx','shifty','shiftz'])
    orientations = particledata.filter(['ux','uy','uz'])
    distances_array = distance_array(coms.values, coms.values, box=dimensions)
    pairs = getPairs(distances_array, cutoff)
    origin_orientations = orientations.values[pairs[:,0]]
    origin_connections = np.subtract(coms.values[pairs[:,1]], coms.values[pairs[:,0]])
    # origin_connections = np.divide(origin_connections, 10)
    projections = np.einsum('ij,ij->i', origin_connections, origin_orientations) # same as (origin_connections * origin_orientations).sum(axis=1) BUT FASTER
    projections_array = np.zeros_like(distances_array)
    pairs_t = pairs.T
    projections_array[tuple(pairs_t)] = projections
    projections_array[tuple([pairs_t[1], pairs_t[0]])] = projections
    sums = np.sum(projections_array, axis=1)
    nums = np.count_nonzero(projections_array, axis=1)
    averages = np.zeros_like(sums)
    averages[np.where(nums>0)] = sums[np.where(nums>0)]/nums[np.where(nums>0)]

    coms = particledata.filter(['x','y','z'])
    distances_array = distance_array(coms.values, coms.values, box=None)
    pairs = getPairs(distances_array, cutoff)
    origin_orientations = orientations.values[pairs[:,0]]
    origin_connections = np.subtract(coms.values[pairs[:,1]], coms.values[pairs[:,0]])
    # origin_connections = np.divide(origin_connections, 10)
    projections = np.einsum('ij,ij->i', origin_connections, origin_orientations) # same as (origin_connections * origin_orientations).sum(axis=1) BUT FASTER
    projections_array = np.zeros_like(distances_array)
    pairs_t = pairs.T
    projections_array[tuple(pairs_t)] = projections
    projections_array[tuple([pairs_t[1], pairs_t[0]])] = projections
    _sums = np.sum(projections_array, axis=1)
    _nums = np.count_nonzero(projections_array, axis=1)
    _averages = np.zeros_like(_sums)
    _averages[np.where(_nums>0)] = _sums[np.where(_nums>0)]/nums[np.where(_nums>0)]

    condition = np.where(np.logical_or(averages<-1, averages>1))
    # df = pd.DataFrame({"old":averages[condition]})
    averages[condition] = _averages[condition]
    # df["new"] = averages[condition]
    # print(df)

    return np.nan_to_num(averages)