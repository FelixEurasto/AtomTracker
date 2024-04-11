import numpy as np
import scipy

def make_cylindrical_grid(R_min, R_max, z_min, z_max, n_R, n_theta, n_z):
    """
    Creates cylindrical grids of given dimensions
    Parameters:
        R_min -- float: Minimum radial distance from origin in xy-plane
        R_max -- float: Maximum radial distance from origin in xy-plane
        z_min -- float: Minimum value of z-coordinate
        z_max -- float: Maximum value of z-coordinate
        n_R -- int: Number of grids in radial direction
        n_theta -- int: Number of grids in theta direction
        n_z -- int: Number of grids in z direction
    Returns:
        R_grids -- np.ndarray: Values of grid points in R-direction
        theta_grids -- np.ndarray: Values of grid points in theta-direction
        z_grids -- np.ndarray: Values of grid points in z-direction
        volumes -- np.ndarray: Volumes of voxels corresponding to each R-value
    """
    R_grids = np.linspace(R_min, R_max, n_R)
    z_grids = np.linspace(z_min, z_max, n_z)
    theta_grids = np.linspace(0, 2*np.pi, n_theta)
    theta_diff = theta_grids[1] - theta_grids[0]
    areas = [1, theta_diff/2 * R_grids[1]**2]
    for i in range(2, n_R - 1):
        areas.append(theta_diff/2 * (R_grids[i]**2 - R_grids[i - 1]**2)) 
    volumes = (z_grids[1] - z_grids[0]) * np.array(areas)

    return R_grids, theta_grids, z_grids, volumes


def cartesian_to_cylindrical(X, Y):
    """
    Converts cartesian coordinates to cylindrical coordinates
    Parameters:
        X -- np.ndarray: Values of x-coordinates
        Y -- np.ndarray: Values of y-coordinates
    Returns:
        R -- np.ndarray: Values of R-coordinates
        thetas -- np.ndarray: Values of theta coordinates
    """
    R = np.linalg.norm(np.concatenate([X.reshape(-1,1),Y.reshape(-1,1)], axis=1), axis=1)
    thetas = np.zeros(R.shape[0])

    thetas[np.where((X > 0) & (Y > 0))[0]] = np.arctan(Y[np.where((X > 0) & (Y > 0))[0]] / X[np.where((X > 0) & (Y > 0))[0]])
    thetas[np.where((X > 0) & (Y < 0))[0]] = 2*np.pi + np.arctan(Y[np.where((X > 0) & (Y < 0))[0]] / X[np.where((X > 0) & (Y < 0))[0]])
    thetas[np.where((X < 0) & (Y > 0))[0]] = np.pi + np.arctan(Y[np.where((X < 0) & (Y > 0))[0]] / X[np.where((X < 0) & (Y > 0))[0]])
    thetas[np.where((X < 0) & (Y < 0))[0]] = np.pi + np.arctan(Y[np.where((X < 0) & (Y < 0))[0]] / X[np.where((X < 0) & (Y < 0))[0]])

    return R, thetas

def cylindrical_to_cartesian(R, theta, z):
    """
    Convert cylindrical coordinates (R, theta, z) to Cartesian coordinates (x, y, z).

    Parameters:
        R (array_like): Array of radial distances.
        theta (array_like): Array of angles in radians.
        z (array_like): Array of z-coordinates.

    Returns:
        ndarray: Array of Cartesian coordinates (x, y, z).
    """
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    return np.column_stack((x, y, z))


def map_to_grid_points(coordinates, R_vals, theta_vals, z_vals):
    """
    Maps every coordinate in 'coordinates' to the closest grid point
    Parameters:
        coordinates -- np.ndarray: Coordinates of the atoms
        R_vals -- np.ndarray: Grid coordinates in R-dimension
        theta_vals -- np.ndarray: Grid coordinates in theta-dimension
        z_vals -- np.ndarray: Grid coordinates in z-dimension
    Returns:
        R_closest -- np.ndarray: Indeces of closest grid points in R-dimension
        theta_closest -- np.ndarray: Indeces of closest grid points in theta-dimension
        z_closest -- np.ndarray: Indeces of closest grid points in z-dimension
    """
    R_distances = scipy.spatial.distance.cdist(coordinates[:,0].reshape(-1,1), R_vals)
    theta_distances = scipy.spatial.distance.cdist(coordinates[:,1].reshape(-1,1), theta_vals)
    z_distances = scipy.spatial.distance.cdist(coordinates[:,2].reshape(-1,1), z_vals)
    R_closest = np.argmin(R_distances, axis=1)
    theta_closest = np.argmin(theta_distances, axis=1)
    z_closest = np.argmin(z_distances, axis=1)    
    
    return R_closest, theta_closest, z_closest

def create_maps(X, R_grids, theta_grids, z_grids, volumes, normalization="within"):
    
    if X.size != 0:
        maps = []
        for coords in X:
            frame_map = np.zeros((1, len(R_grids), len(theta_grids), len(z_grids)))
            R_coordinates, theta_coordinates = cartesian_to_cylindrical(coords[:,0], coords[:,1])
            cylindrical_coordinates = np.concatenate([R_coordinates.reshape(-1,1),
                                                    theta_coordinates.reshape(-1,1),
                                                    coords[:,2].reshape(-1,1)], axis=1)
            R_closest, theta_closest, z_closest = map_to_grid_points(cylindrical_coordinates, R_grids.reshape(-1,1),
                                                                    theta_grids.reshape(-1,1), z_grids.reshape(-1,1))
            R_within_range = np.where(R_closest != n_R)[0]
            np.add.at(frame_map, (0, R_closest[R_within_range],
                                theta_closest[R_within_range],
                                z_closest[R_within_range]), 1)
            if normalization == "within":
                if len(R_within_range) != 0:
                    frame_map /= len(R_within_range)
            elif normalization == "all":
                frame_map /= len(coords)
            else:
                print(f"Normalization {normalization} not known! Use either 'within' or 'all'.")
                break
            frame_map /= volumes.reshape(1,-1,1,1)
            maps.append(frame_map)
        maps = np.concatenate(maps, axis=0)
        return maps
    else:
        return None


