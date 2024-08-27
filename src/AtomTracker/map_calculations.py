import numpy as np
import scipy
from ..AtomTracker import position_calculations
from ..AtomTracker import plotting
import MDAnalysis as mda
from MDAnalysis.analysis import align
import os
from collections import defaultdict

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

def create_maps(X, R_grids, theta_grids, z_grids, volumes, normalization="within", verbose=False):
    
    if X.size != 0:
        maps = []
        for INDEX, coords in enumerate(X):
            print(f"At map {INDEX+1}/{len(X)}", end="\r")
            frame_map = np.zeros((1, len(R_grids) - 1, len(theta_grids), len(z_grids)))
            R_coordinates, theta_coordinates = cartesian_to_cylindrical(coords[:,0], coords[:,1])
            cylindrical_coordinates = np.concatenate([R_coordinates.reshape(-1,1),
                                                    theta_coordinates.reshape(-1,1),
                                                    coords[:,2].reshape(-1,1)], axis=1)
            R_closest, theta_closest, z_closest = map_to_grid_points(cylindrical_coordinates, R_grids.reshape(-1,1),
                                                                    theta_grids.reshape(-1,1), z_grids.reshape(-1,1))
            R_within_range = np.where(R_closest != len(R_grids) - 1)[0] # Filter out atoms that are not within radius
            np.add.at(frame_map, (0, R_closest[R_within_range],
                                theta_closest[R_within_range],
                                z_closest[R_within_range]), 1)
            if normalization == "within":
                if len(R_within_range) != 0:
                    frame_map /= len(R_within_range)
            elif normalization == "all":
                frame_map /= len(coords)
            elif not normalization:
                pass
            else:
                print(f"Normalization {normalization} not known! Use either 'within' or 'all'.")
                break
            frame_map /= volumes.reshape(1,-1,1,1)
            maps.append(frame_map)
        maps = np.concatenate(maps, axis=0)
        return maps
    else:
        return None


class Mapper():
    def __init__(
                self,
                save_path,
                map_selections,
                other_selections,
                ref_structure,
                alignment_selection,
                fit_structures,
                centering_selection,
                function,
                R_min,
                R_max,
                n_R=30,
                n_theta=30,
                n_z=30,
                normalization=None,
                skip=1
                ):
        self.map_selections = map_selections
        self.save_path = save_path
        self.other_selections = other_selections
        self.ref_structure = ref_structure
        self.alignment_selection = alignment_selection
        self.fit_structures = fit_structures
        self.centering_selection = centering_selection
        self.function = function
        self.R_min = R_min
        self.R_max = R_max
        self.n_R = n_R
        self.n_theta = n_theta
        self.n_z = n_z
        self.normalization = normalization
        self.skip = skip
        self.maps = None
        self.fvals = None
        self.other_coordinates = None
        self.R_grids = None
        self.theta_grids = None
        self.z_grids = None
    
    def run(self, universe, verbose=True):
        if not os.path.isdir(self.save_path):
            print(f"Directory {self.save_path} not found, creating it now ...")
            os.makedirs(self.save_path)
        if verbose:
            print(f"Calculating coordinate info ...")
        if self.ref_structure:
            ref_struct = mda.Universe(self.ref_structure)
        else:
            ref_struct = self.universe.copy()
        ref_struct.atoms.positions -= ref_struct.select_atoms(self.centering_selection).center_of_mass()
        all_selections = " or ".join(self.map_selections)
        atoms = universe.select_atoms(all_selections)
        if atoms.n_atoms > 0:
            z_min_max = []
            for ts in universe.trajectory[::int(universe.trajectory.n_frames / 100)]:
                if self.fit_structures:
                    align.alignto(universe, ref_struct, select=self.alignment_selection)
                else:
                    universe.atoms.positions -= universe.select_atoms(self.centering_selection).center_of_mass()
                z_min_max.append(atoms.positions[:,2].min())
                z_min_max.append(atoms.positions[:,2].max())
            z_min = np.min(z_min_max)
            z_max = np.max(z_min_max)
            R_grids, theta_grids, z_grids, volumes = make_cylindrical_grid(R_min=self.R_min, R_max=self.R_max,
                                                                                z_min=z_min, z_max=z_max,
                                                                                n_R=self.n_R + 1, n_theta=self.n_theta,
                                                                                n_z=self.n_z)
            self.R_grids = R_grids
            self.theta_grids = theta_grids
            self.z_grids = z_grids
            maps = defaultdict(np.ndarray)

            coordinates = position_calculations.calculate_positions(universe=universe,
                                                map_selections=self.map_selections,
                                                skip=self.skip, fit_structures=self.fit_structures,
                                                reference_structure=self.ref_structure,
                                                other_selections=self.other_selections,
                                                centering_selection=self.centering_selection,
                                                alignment_selection=self.alignment_selection,
                                                verbose=verbose)

            for sel in self.map_selections:
                if verbose:
                    print(f"Calculating density maps for selection '{sel}' ...")
                X = coordinates[sel]
                sel_maps = create_maps(X=X,
                                        R_grids=R_grids,
                                        theta_grids=theta_grids,
                                        z_grids=z_grids,
                                        volumes=volumes,
                                        normalization=self.normalization,
                                        verbose=verbose)
                if sel_maps is not None: # Check if map for this selection exists
                    maps[sel] = sel_maps

            if self.other_selections:
                self.other_coordinates = dict([(sel, coordinates[sel]) for sel in self.other_selections if coordinates[sel].size > 0])
            
            if verbose: 
                print(f"Calculating function values ...")
            if self.function:
                fvals = np.array([self.function(universe) for ts in universe.trajectory[::self.skip]])
                self.fvals = fvals

            if len(maps.keys()) > 0:
                self.maps = maps
            if self.save_path:
                if verbose:
                    print(f"Saving data into {self.save_path} ...")
                self.save_data()
            print("Done!")
        else:
            print(f"None of {','.join(self.map_selections)} found in this Universe!")
        return self

    def save_data(self):
        if self.maps is not None:
            for selection, data in self.maps.items():
                np.save(f"{self.save_path}/{selection.replace(' ', '_')}_maps.npy", data)
        if self.fvals is not None:
            if self.fvals.size > 0:
                np.save(f"{self.save_path}/function_values.npy", self.fvals)
        if self.other_coordinates is not None:
            for selection, data in self.other_coordinates.items():
                np.save(f"{self.save_path}/{selection.replace(' ', '_')}_coordinates.npy", data)

    def interactive_3d(self, density_sel, coordinates, fval_range=(-10,10), n_frames=21, threshold=0.4, res=(10,10,10)):
        densities = self.maps[density_sel]
        fig = plotting._3d_plot(maps=densities,
                                coordinates=coordinates,
                                fvals=self.fvals,
                                R_grids=self.R_grids,
                                theta_grids=self.theta_grids,
                                z_grids=self.z_grids,
                                res=res,
                                fval_range=fval_range,
                                n_frames=n_frames,
                                threshold=threshold,
                                title=density_sel)
        fig.show()

