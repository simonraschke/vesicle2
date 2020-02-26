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



# iterate a file linewise
# check line for keyword
# if keyword is found return value after seperator
def fileStringFromKeyword(filepath, keyword, seperator='='):
    found_counter = 0
    assert(filepath)
    value = [np.NaN]
    try:
        with open(filepath) as FILE:
            for line in FILE:
                if keyword in line:
                    found_counter += 1
                    value = [str(line.split(seperator,1)[1].rstrip())]
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
        'fga_mode': fileStringFromKeyword(filename, 'fga_mode', '='),
        'simulation_mode': fileStringFromKeyword(filename, 'simulation_mode', '='),
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
        'acceptance_position_target': fileValueFromKeyword(filename, 'acceptance_position_target', '='),
        'sw_orientation_min': fileValueFromKeyword(filename, 'sw_orientation_min', '='),
        'sw_orientation_max': fileValueFromKeyword(filename, 'sw_orientation_max', '='),
        'acceptance_orientation_target': fileValueFromKeyword(filename, 'acceptance_orientation_target', '='),
        'cell_min_edge': fileValueFromKeyword(filename, 'cell_min_edge', '='),
        'max_cells_dim': fileValueFromKeyword(filename, 'max_cells_dim', '='),
        'skip': fileValueFromKeyword(filename, 'skip', '=')
    }



# def getSubclusterLabels(ID, group, eps):
#     #skip if noise
#     if ID == -1:
#         return np.zeros( len(group.index), dtype=int )
#     else:
#         # arange a DBSCAN without PBC to get subclusters
#         coms_subcluster = pd.concat([group['x'], group['y'], group['z']], axis=1)
#         distances_array_subcluster = distance_array(coms_subcluster.values, coms_subcluster.values, box=None)
#         subclusterd = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=-1).fit(distances_array_subcluster)
#         return subclusterd.labels_



def getSubclusterLabels(group, eps):
    group["subcluster"] = 0
    if group["cluster"].unique() == -1:
        return group
    #skip if noise
    # if ID == -1:
    #     return np.zeros( len(group.index), dtype=int )
    else:
        # arange a DBSCAN without PBC to get subclusters
        coms = group.filter(["x","y","z"])
        distances_array = distance_array(coms.values, coms.values, box=None)
        group["subcluster"] = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=-1).fit(distances_array).labels_
        # print(group)
        return group

        # coms_subcluster = pd.concat([group['x'], group['y'], group['z']], axis=1)
        # distances_array_subcluster = distance_array(coms_subcluster.values, coms_subcluster.values, box=None)
        # subclusterd = DBSCAN(min_samples=1, eps=eps, metric="precomputed", n_jobs=-1).fit(distances_array_subcluster)
        # return subclusterd.labels_



def getShiftedCoordinates(group, eps, dimensions):
    group[["shiftx","shifty","shiftz"]] = group[["x","y","z"]] 
    if group["cluster"].unique() == -1:
        return group

    # get largest subcluster
    unique, counts = np.unique(group["subcluster"], return_counts=True)
    if len(unique) == 1:
        return group
    max_subclusterID = unique[counts == np.max(counts)][0]
    # calculate shifts per subcluster
    centers = group.groupby("subcluster")[['x','y','z']].mean()
    shifts = np.round(( -centers + centers.loc[max_subclusterID] )/dimensions[:3]).astype(int)
    shifts *= dimensions[:3]
    group[["shiftx","shifty","shiftz"]] += shifts.loc[group["subcluster"]].values
    # print(group)
    return group



def getClusterVolume(group, eps, pps):
    # volume = 0
    group["volume"] = np.pi * (eps**3) * 4/3
    if group["cluster"].unique() == -1:
        return group
    else:
        # generate meshgrid around cluster particles plus threshold
        x_vector = np.arange(np.min(group["shiftx"])-eps, np.max(group["shiftx"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)
        y_vector = np.arange(np.min(group["shifty"])-eps, np.max(group["shifty"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)
        z_vector = np.arange(np.min(group["shiftz"])-eps, np.max(group["shiftz"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)

        # make it recursive if meshgrid becomes too large
        NNN = 80
        if len(x_vector) > NNN or len(y_vector) > NNN or len(z_vector) > NNN:
            if pps-1 < 1:
                return getClusterVolume(group, eps, pps/2)
            else:
                return getClusterVolume(group, eps, pps-1)

        xx, yy, zz = np.meshgrid(x_vector, y_vector, z_vector)
        # stack them together as array of 3D points
        meshgrid = np.stack((xx.ravel(), yy.ravel(), zz.ravel()), axis=1)
        #calculate the distance array with centres of masses of particles
        coms_cluster = group.filter(['shiftx','shifty','shiftz'])
        getClusterVolume.distances_array_volume = distance_array(meshgrid, coms_cluster.values, box=None).astype(np.float32)
        # check if any point in distance array row is close enough, then reshape to meshgrid
        # result is a binary meshgrid with 1 for the cluster shell region
        isclose = np.where(getClusterVolume.distances_array_volume <= eps, True, False).any(axis=1).reshape(xx.shape[0], yy.shape[1], zz.shape[2])
        # fill hole inside the shell region
        isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
        # calc volum from all points inside cluster
        group["volume"] = ((1.0/pps)**3)*np.count_nonzero(isclose)
        return group


# def getClusterVolume(ID, group, eps, pps):
#     # volume = 0
#     if ID == -1:
#         # single particles as spheres
#         return np.pi * (eps**3) * 4/3
#     else:
#         # generate meshgrid around cluster particles plus threshold
#         x_vector = np.arange(np.min(group["shiftx"])-eps, np.max(group["shiftx"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)
#         y_vector = np.arange(np.min(group["shifty"])-eps, np.max(group["shifty"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)
#         z_vector = np.arange(np.min(group["shiftz"])-eps, np.max(group["shiftz"])+eps+1.0/pps, 1.0/pps, dtype=np.float16)

#         # make it recursive if meshgrid becomes too large
#         NNN = 80
#         if len(x_vector) > NNN or len(y_vector) > NNN or len(z_vector) > NNN:
#             if pps-1 < 1:
#                 return getClusterVolume(ID, group, eps, pps/2)
#             else:
#                 return getClusterVolume(ID, group, eps, pps-1)

#         xx, yy, zz = np.meshgrid(x_vector, y_vector, z_vector)
#         # stack them together as array of 3D points
#         meshgrid = np.stack((xx.ravel(), yy.ravel(), zz.ravel()), axis=1)
#         #calculate the distance array with centres of masses of particles
#         coms_cluster = pd.concat([group['shiftx'], group['shifty'], group['shiftz']], axis=1)
#         getClusterVolume.distances_array_volume = distance_array(meshgrid, coms_cluster.values, box=None).astype(np.float32)
#         # check if any point in distance array row is close enough, then reshape to meshgrid
#         # result is a binary meshgrid with 1 for the cluster shell region
#         isclose = np.where(getClusterVolume.distances_array_volume <= eps, True, False).any(axis=1).reshape(xx.shape[0], yy.shape[1], zz.shape[2])
#         # fill hole inside the shell region
#         isclose = ndimage.morphology.binary_fill_holes(isclose).astype(bool)
#         # calc volum from all points inside cluster
#         return ((1.0/pps)**3)*np.count_nonzero(isclose)


        # z,x,y = isclose.nonzero()
        # ax.scatter(x+,y,z, s=2)

        


# def getOrder(ID, group):
    # if ID == -1:
    #     return group["order"].values
    # else:
    #     shifted_coms = group.filter(['shiftx','shifty','shiftz'])
    #     normalized_orientations = group.filter(['ux','uy','uz'])
    #     return (normalized_orientations.values*normalize(shifted_coms.sub(shifted_coms.mean().values))).sum(axis=1)



def getPairs(distances_array, max_cutoff, same=False):
    valid = distances_array < max_cutoff
    np.fill_diagonal(valid, same)
    pairs = np.column_stack(np.where(valid))
    if len(pairs) == 0:
        return []
    pairs = np.sort(pairs, axis=1)
    unique, index = np.unique(pairs, axis=0, return_index=True)
    return pairs[index]



def getNormedPairDistanceVectors(coms, pairs, dimensions):
    dist_vecs = np.subtract(coms.loc[pairs[:,1]].values, coms.loc[pairs[:,0]].values)
    dist_vecs = np.subtract(dist_vecs, np.multiply(dimensions[:3], np.round(dist_vecs/dimensions[:3])))
    return pd.DataFrame((dist_vecs.T / np.linalg.norm(dist_vecs, axis=1)).T, columns=["x","y","z"]), np.linalg.norm(dist_vecs, axis=1)



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
    if len(pairs) == 0:
        return np.full((len(particledata.index),1), np.nan, dtype=np.float32)
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

    # nonan = np.where(~np.isnan(averages))
    # condition = np.where(np.logical_or(averages[nonan]<-1, averages[nonan]>1))
    condition = np.where(np.logical_or(averages<-1, averages>1))
    
    averages[condition] = _averages[condition]
    return np.nan_to_num(averages)



def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:,-1]



def normalVector(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return svd(M)[0][:,-1]



def getOrder(group, attributes):
    group["order"] = np.nan
    if group["cluster"].unique() == -1:
        return group
    else:
        if attributes.get("system.gamma") < 1e-3:
            if len(group.index) <= 3:
                return group
            shifted_coms = group.filter(['shiftx','shifty','shiftz'])
            normalized_orientations = group.filter(['ux','uy','uz'])
            center, normal = planeFit(shifted_coms.values.T)
            data = np.sum(normalized_orientations.values*normal, axis=1)
            if data.mean() < 0:
                data_compare = data * (-1)
                group["order"] = [data, data_compare][np.argmax([data.mean(), data_compare.mean()])]
            else:
                group["order"] = data
            return group
        else :
            shifted_coms = group.filter(['shiftx','shifty','shiftz'])
            normalized_orientations = group.filter(['ux','uy','uz'])
            data = np.sum(normalized_orientations.values*normalize(shifted_coms.sub(shifted_coms.mean()).values, copy=False), axis=1)
            group["order"] = data
            return group



def getSurfaceTension(particledata, dimensions, epot_calc, cutoff=3.0, dr=0.1):
    """
    find closest neighbours closer 3 sigma
    calculate energy of particle
    set particle away from local plane in upwards direction
    calculate energy of particle
    calculate energy difference
    """
    from pprint import pprint

    # cluster_centers = particledata.groupby("cluster")[['shiftx','shifty','shiftz']].mean()
    # # print(cluster_centers)
    # particledata["cluster_x"] = particledata.cluster.apply(lambda x: cluster_centers.loc[x][0])
    # particledata["cluster_y"] = particledata.cluster.apply(lambda x: cluster_centers.loc[x][1])
    # particledata["cluster_z"] = particledata.cluster.apply(lambda x: cluster_centers.loc[x][2])
    # # print(cluster_centers.loc[-1])
    # # particledata.loc[cluster_centers.index, "cluster_x"] = cluster_centers.loc[cluster_centers.index, "shiftx"]

    # coms = particledata.filter(['shiftx','shifty','shiftz'])
    # orientations = particledata.filter(['ux','uy','uz'])

    # distances_array = distance_array(coms.values, coms.values, box=dimensions)
    # pairs = getPairs(distances_array, cutoff)
    
    coms = particledata.filter(['shiftx','shifty','shiftz'])
    orientations = particledata.filter(['ux','uy','uz'])
    distances_array = distance_array(coms.values, coms.values, box=dimensions)
    pairs = getPairs(distances_array, cutoff)

    for origin_id in np.unique(pairs.T[0]):
        # get all neighbour ids
        origin_neighbours = pairs.T[1][np.where(pairs.T[0] == origin_id)]
        # if origin_neighbours.size < 10: continue
        # all ids of plane 
        all_ids = np.append(origin_neighbours, origin_id)
        print(all_ids)
        point_cloud = coms.loc[all_ids].values.T
        #if more dimenstions than points
        if point_cloud.shape[0] >= point_cloud.shape[1]: continue
        
        normal = normalVector(point_cloud)
        print(normal)
        orthogonal = np.random.rand(3)
        orthogonal /= np.linalg.norm(orthogonal)
        orthogonal_xx = np.cross(normal, orthogonal)
        orthogonal_yy = np.cross(normal, orthogonal_xx)
        print(orthogonal_xx)
        print(orthogonal_yy)
        
        # initial
        epot_coms = coms.loc[all_ids].reset_index(drop=True)
        epot_orientations = orientations.loc[all_ids].reset_index(drop=True)
        epot_before = epot_calc.get(epot_coms, epot_orientations, dimensions, cutoff=cutoff*10, ret="epot")[-1]
        print(particledata.loc[origin_id]["epot"], epot_before)
        assert np.abs(particledata.loc[origin_id]["epot"] - epot_before) < 5e-1
        # print(epot_before)

        # pxx
        epot_coms.iloc[-1] += orthogonal_xx*dr
        d_epot_xx = epot_calc.get(epot_coms, epot_orientations, dimensions, cutoff=cutoff*10, ret="epot")[-1] - epot_before
        epot_coms.iloc[-1] -= orthogonal_xx*dr

        # pyy
        epot_coms.iloc[-1] += orthogonal_yy*dr
        d_epot_yy = epot_calc.get(epot_coms, epot_orientations, dimensions, cutoff=cutoff*10, ret="epot")[-1] - epot_before
        epot_coms.iloc[-1] -= orthogonal_yy*dr

        # pzz
        epot_coms.iloc[-1] += normal*dr
        d_epot_zz = epot_calc.get(epot_coms, epot_orientations, dimensions, cutoff=cutoff*10, ret="epot")[-1] - epot_before
        epot_coms.iloc[-1] -= normal*dr

        epot_coms.index = all_ids

        tensor = np.zeros((3,3))
        tensor[0,0], tensor[1,1], tensor[2,2] = d_epot_xx, d_epot_yy, d_epot_zz
        print(tensor)
        print()

    sys.exit()
    

    return np.nan_to_num(np.zeros_like(particledata.index))



# plane : [[-x,-y,-z],[+x,+y,+z]]
def isParticleInStructure(df, attributes, dimensions, fga_mode):
    xyz = np.array(dimensions[:3])
    if fga_mode == "plane":
        plane = np.array([xyz/2 - attributes["system.plane_edge"]/2, xyz/2 + attributes["system.plane_edge"]/2])
        plane[0,2] = xyz[2]/2 - attributes["system.ljsigma"]
        plane[1,2] = xyz[2]/2 + attributes["system.ljsigma"]
        xcond = np.logical_and(df["x"] >= plane[0,0], df["x"] <= plane[1,0])
        ycond = np.logical_and(df["y"] >= plane[0,1], df["y"] <= plane[1,1])
        zcond = np.logical_and(df["z"] >= plane[0,2], df["z"] <= plane[1,2])
        return np.logical_and.reduce((xcond,ycond,zcond)), np.abs(np.cumprod(np.subtract(plane[1], plane[0])))[-1]
    elif fga_mode == "sphere":
        center = np.array([xyz/2])
        radius = float(attributes["system.ljsigma"])**(1.0/6) / (2.0*np.sin(attributes["system.gamma"])) + attributes["system.ljsigma"]
        return (distance_array(center, df.filter(["x","y","z"]).values, box=dimensions) <= radius).ravel(), np.pi*4/3*radius**3*attributes["system.frame_guides_grid_edge"]**3
    elif fga_mode == "pair":
        plane = np.array([xyz/2 - attributes["system.plane_edge"]/2, xyz/2 + attributes["system.plane_edge"]/2])
        plane[0,1] = xyz[1]/2 - attributes["system.ljsigma"]
        plane[1,1] = xyz[1]/2 + attributes["system.ljsigma"]
        plane[0,2] = xyz[2]/2 - attributes["system.ljsigma"]
        plane[1,2] = xyz[2]/2 + attributes["system.ljsigma"]
        xcond = np.logical_and(df["x"] >= plane[0,0], df["x"] <= plane[1,0])
        ycond = np.logical_and(df["y"] >= plane[0,1], df["y"] <= plane[1,1])
        zcond = np.logical_and(df["z"] >= plane[0,2], df["z"] <= plane[1,2])
        return np.logical_and.reduce((xcond,ycond,zcond)), np.abs(np.cumprod(np.subtract(plane[1], plane[0])))[-1]



def isParticleInStructureEnvironment(df, attributes, dimensions, fga_mode):
    xyz = np.array(dimensions[:3])
    if fga_mode == "plane":
        plane = np.array([xyz/2 - attributes["system.plane_edge"]/2, xyz/2 + attributes["system.plane_edge"]/2])
        plane[0,2] = xyz[2]/2 - attributes["system.ljsigma"]*3
        plane[1,2] = xyz[2]/2 + attributes["system.ljsigma"]*3
        xcond = np.logical_and(df["x"] >= plane[0,0], df["x"] <= plane[1,0])
        ycond = np.logical_and(df["y"] >= plane[0,1], df["y"] <= plane[1,1])
        zcond = np.logical_and(df["z"] >= plane[0,2], df["z"] <= plane[1,2])
        return np.logical_and.reduce((xcond,ycond,zcond))
    elif fga_mode == "sphere":
        center = np.array([xyz/2])
        radius = float(attributes["system.ljsigma"])**(1.0/6) / (2.0*np.sin(attributes["system.gamma"])) + attributes["system.ljsigma"]*3
        return (distance_array(center, df.filter(["x","y","z"]).values, box=dimensions) <= radius).ravel()
    elif fga_mode == "pair":
        plane = np.array([xyz/2 - attributes["system.plane_edge"]/2, xyz/2 + attributes["system.plane_edge"]/2])
        plane[0,1] = xyz[1]/2 - attributes["system.ljsigma"]*3
        plane[1,1] = xyz[1]/2 + attributes["system.ljsigma"]*3
        plane[0,2] = xyz[2]/2 - attributes["system.ljsigma"]*3
        plane[1,2] = xyz[2]/2 + attributes["system.ljsigma"]*3
        xcond = np.logical_and(df["x"] >= plane[0,0], df["x"] <= plane[1,0])
        ycond = np.logical_and(df["y"] >= plane[0,1], df["y"] <= plane[1,1])
        zcond = np.logical_and(df["z"] >= plane[0,2], df["z"] <= plane[1,2])
        return np.logical_and.reduce((xcond,ycond,zcond))



def isParticleInStructureCluster(df, fga_mode):

    in_structure = df[df["in_structure"]]
    # print(in_structure)
    uniques, counts = np.unique(in_structure[in_structure["cluster"] != -1]["cluster"], return_counts=True)
    # print(uniques, counts)
    # largest_cluster_size = in_structure[in_structure["cluster"] != -1].groupby(["cluster"])["cluster"].count().max()
    if len(counts) == 0 or max(counts) < 2:
        return False
    largest_cluster_ID = uniques[np.argmax(counts)]
    # print(largest_cluster_ID)
    # sys.exit()
    return df["cluster"] == largest_cluster_ID




class EpotCalculator(object):
    def __init__(self, attributes):
        self.sigma = attributes["system.ljsigma"]
        self.epsilon = attributes["system.ljepsilon"]
        self.kappa = attributes["system.kappa"]
        self.gamma = attributes["system.gamma"]
        self.a = 1.0 + self.kappa*np.sin(self.gamma*np.pi/180)
        self.b = 1.0 - self.kappa*np.sin(self.gamma*np.pi/180)
        self.c = np.sqrt(1+self.kappa**2*np.cos(self.gamma*np.pi/180))  

    
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
            pairs = getPairs(distances_array, cutoff)
            normed_dist_vecs, dist_norms = getNormedPairDistanceVectors(coms, pairs, dimensions)
            res1_u = np.multiply(np.take(orientations.values, pairs[:,0], axis=0), self.kappa/2)
            res2_u = np.multiply(np.take(orientations.values, pairs[:,1], axis=0), self.kappa/2)
            chi =  np.power((np.linalg.norm(-res1_u + normed_dist_vecs + res2_u, axis=1) - self.a), 2)
            chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs - res2_u, axis=1) - self.b), 2)
            chi += np.power((np.linalg.norm(-res1_u + normed_dist_vecs - res2_u, axis=1) - self.c), 2)
            chi += np.power((np.linalg.norm( res1_u + normed_dist_vecs + res2_u, axis=1) - self.c), 2)
            chi = chi * np.power(self.sigma/dist_norms, 6)
            chi_array = np.zeros_like(distances_array)
            pairs_t = pairs.T
            chi_array[tuple(pairs_t)] = chi
            chi_array[tuple([pairs_t[1], pairs_t[0]])] = chi
            # chi_array = np.where(chi_array > 1e-5, chi_array * np.power(self.sigma/distances_array, 6), 0)
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
            return self.get(coms, orientations, dimensions, cutoff=cutoff, ret="epot", distances_array=distances_array), \
                   self.get(coms, orientations, dimensions, cutoff=cutoff, ret="chi",  distances_array=distances_array)

        elif ret == "chi+epot":
            return self.get(coms, orientations, dimensions, cutoff=cutoff, ret="chi",  distances_array=distances_array), \
                   self.get(coms, orientations, dimensions, cutoff=cutoff, ret="epot", distances_array=distances_array)



# from dataclasses import dataclass

class Plane:
    ID: int
    x: int
    xmin: float
    xmax: float
    y: int
    ymin: float
    ymax: float
        
    def __init__(
            self,
            ID: int,
            x: int,
            xmin: float,
            xmax: float,
            y: int,
            ymin: float,
            ymax: float,
        ) -> None:
        self.ID = ID
        self.x = x
        self.xmin = xmin
        self.xmax = xmax
        self.y = y
        self.ymin = ymin
        self.ymax = ymax

    def __repr__(self) -> str:
        return f"Plane: ID {self.ID},  x {self.x}, xmin {self.xmin},  xmax {self.xmax},  y {self.y},  ymin {self.ymin},  ymax {self.ymax}"

    def __hash__(self) -> int:
        return hash((self.ID, self.x, self.xmin, self.xmax, self.y, self.ymin, self.ymax))

    def __eq__(self, other) -> bool:
        if not isinstance(other, Plane):
            return NotImplemented
        return (
            (self.ID, self.x, self.xmin, self.xmax, self.y, self.ymin, self.ymax)== 
            (other.ID, other.x, other.xmin, other.xmax, other.y, other.ymin, other.ymax))
    
    @property
    def dx(self):
        return self.xmax-self.xmin

    @property
    def dy(self):
        return self.ymax-self.ymin

    @property
    def area(self):
        return self.dx*self.dy
    
    def contains(self, point):
        return self.xmin <= point[0] <= self.xmax and self.ymin <= point[1] <= self.ymax

    def isNeighboursOf(self, other):
        return abs(self.x - other.x) == 1 and abs(self.y - other.y) == 0 or abs(self.y - other.y) == 1 and abs(self.x - other.x) == 0



class Cuboid(Plane):
    z: int
    zmin: float
    zmax: float
    
    def __init__(
        self,
        ID: int,
        x: int,
        xmin: float,
        xmax: float,
        y: int,
        ymin: float,
        ymax: float,
        z: int,
        zmin: float,
        zmax: float
        ) -> None:
        self.z = z
        self.zmin = zmin
        self.zmax = zmax
        super(Cuboid, self).__init__(ID, x, xmin, xmax, y, ymin, ymax)

    def __repr__(self) -> str:
        return f"Coboid: {super(Cuboid, self).__repr__()}, z{self.z}, zmin{self.zmin}, xzmax{self.zmax}"

    def __hash__(self) -> int:
        return hash((super(Cuboid, self), self.z, self.zmin, self.zmax))

    def __eq__(self, other) -> bool:
        if not isinstance(other, Cuboid):
            return NotImplemented
        return (
            (super(Cuboid, self), self.z, self.zmin, self.zmax) == 
            (super(Cuboid, self), other.z, other.zmin, other.zmax))
    
    
    @property
    def dz(self):
        return self.zmax-self.zmin

    @property
    def volume(self):
        return self.dx*self.dy*self.dz

    @property
    def center(self):
        return [
            self.xmin + self.dx/2,
            self.ymin + self.dy/2,
            self.zmin + self.dz/2
        ]
    
    def contains(self, point):
        return self.xmin < point[0] <= self.xmax and self.ymin < point[1] <= self.ymax and self.zmin < point[2] <= self.zmax
    
    
def generateDomains(plane_edge, num_dim, center, ):
    domains = []
    domain_edge = plane_edge/num_dim
    for x in range(num_dim):
        for y in range(num_dim):
            domains.append(
                Cuboid(
                    # (num_dim*num_dim+1) - (x*num_dim+y+1),
                    (num_dim*num_dim-1) - (x*num_dim+y),
                    x, 
                    center[0] - plane_edge/2 + domain_edge*x,
                    center[0] - plane_edge/2 + domain_edge*(x+1),
                    y, 
                    center[1] - plane_edge/2 + domain_edge*y,
                    center[1] - plane_edge/2 + domain_edge*(y+1),
                    0, 
                    # FIXME: hardocded == BAD
                    center[2] - 1.5,
                    center[2] + 1.5
                )
            )
    return domains
    
    
def generateDomainsInternal(plane_edge, num_dim, center, sigma):
    domains = []
    domain_edge = plane_edge/num_dim
    count = (num_dim-1)**2-1
    for x in range(num_dim-1):
        for y in range(num_dim-1):
            domains.append(
                Cuboid(
                    # (num_dim*num_dim+1) - (x*num_dim+y+1),
                    count,
                    x,
                    center[0] - (num_dim-1)*(domain_edge)/2 + domain_edge*x,
                    center[0] - (num_dim-1)*(domain_edge)/2 + domain_edge*(x+1),
                    y, 
                    center[1] - (num_dim-1)*(domain_edge)/2 + domain_edge*y,
                    center[1] - (num_dim-1)*(domain_edge)/2 + domain_edge*(y+1),
                    0, 
                    # FIXME: hardocded == BAD
                    center[2] - 1.5,
                    center[2] + 1.5
                )
            )
            count -= 1
    return domains




def getDomainIDofPoint(point, domains):
    for d in domains:
        if d.contains(point.values):
            return d.ID
    return -1



def getStructureDomainID(particledata, domains):
    return particledata.filter(["shiftx","shifty","shiftz"]).apply(axis=1, func=lambda point: getDomainIDofPoint(point, domains)).astype(np.int8)

def getStructureDomainIDInternal(particledata, domains):
    return particledata.filter(["shiftx","shifty","shiftz"]).apply(axis=1, func=lambda point: getDomainIDofPoint(point, domains)).astype(np.int8)


def _internal_squared_distance(x0, x1, dims):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dims, delta - dims, delta)
    # return np.linalg.norm(delta, axis=1)
    return (delta ** 2).sum(axis=1)



def _internal_distance(squared_distances):
    return np.linalg.norm(squared_distances, axis=1)



def getMSD(origin, compare, dimensions):
    return _internal_squared_distance(origin, compare, dimensions[:3])



def correctForPBC(last, compare, dimensions):
    delta = compare[["x","y","z"]] - last[["x","y","z"]]
    shifts = np.round(delta/dimensions[:3]).astype(int)
    corrected_delta = delta - shifts*dimensions[:3]
    data = last[["xnopbc","ynopbc","znopbc"]] + corrected_delta.values
    return pd.DataFrame.from_records(data)