# A tool to build density maps of atoms from MD simulations


## Required Python Modules
- NumPy
- matplotlib
- MDAnalysis
- scikit-learn

## Usage

All of the required settings are defined in /src/AtomTracker/config.py. Variables with `*` cannot be `None`.

First, define the directory paths from which the simulation trajectories and corresponding .gro/.pdb files are fetched.

| Variable | Description | Default value |
| --- | - | - |
| `traj_paths*` | List of paths from which to fetch .xtc files |  - |
| `gro_paths*` | List of paths from which to fetch corresponding .gro/.pdb files | - |
| `save_paths*` | List of paths to which the corresponding results are saved | - |


Next, define all settings needed for the calculations of the density maps.

| Variable | Description | Default value |
| --- | - | - |
| `map_selections*` | List of MDAnalysis selection strings for which density maps are calculated |  - |
| `other_selections` | List of MDAnalysis selection strings for which coordinates are saved at every frame. | `None` |
| `reference_structure` | Path to .gro/.pdb file which will be used as reference structure for fitting. If `None` and `fit_structures` is `True`, each frame will be aligned to the first frame in its trajectory. | `None` |
| `alignment_selection` | MDAnalysis selection string that is used to align each frame to `reference_structure`, if `fit_structures` is `True` | `None` |
| `fit_structures*` | Whether to use `alignment_selection` to fit every frame to `reference_structure` | `False` |
| `centering_selection*` | MDAnalysis string of structure whose center of mass (COM) is shifted to the origin (center of the grid) on each frame | - |
| `R_min*` | Minimum radial distance (in Angstroms) from the COM of `centering_selection` for which atoms are considered in the density maps | `0` |
| `R_max*` | Maximum radial distance (in Angstroms) from the COM of `centering_selection` for which atoms are considered in the density maps | - |
| `n_R*` | Number of grid points in the radial direction | - |
| `n_theta*` | Number of grid points in the (azimuthal) $\theta$-direction | - |
| `n_z*` | Number of grid points in the z-direction | - |
| `normalization` | How to normalize density values. `"all"` normalizes by dividing by the total number of atoms in the map selection. `"within"` divides by the number of atoms seen within the radius limits. If `None`, no normalization is done.  | `"all"` |
| `skip*` | Number of frames to skip when computing the density maps | 10 |


Along with the density maps, you may wish to evaluate some function related to the system at each frame. This is done by assigning the `function` variable to some function. If `function` is `None`, no function is evaluated when computing the maps. 


If `function` is not `None`, you can optionally create movies showing the trends in how atom densities shift, as the function takes different values. This is done by fitting a partial least squares (PLS) model using the density maps to predict the function value. PLS then allows the interpolation of density maps from particular function values, producing best estimates of what the corresponding density maps should look like.

The movies are structured such that the z-dimension of the system is split in 3, producing a "lower section", "central section", and "upper section". The density shifts in all of these regions are animated separately. These movies are calculated separately for each selection in `map_selections`.

To make density, movies you need to define the following variables:

| Variable | Description | Default value |
| --- | - | - |
| `create_movies*` | Whether to create movies | `False` |
| `movie_n_frames*` | Number of frames in the movie | `100` |
| `movie_length*` | Length of movie in seconds | `5` |
| `function_min_and_max*` | Tuple containing the range of function values for which atom densities are interpolated (there will be `movie_n_frames` interpolations) | (0,1) |





