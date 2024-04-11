# A tool to build density maps of atoms from MD simulations


## Required Python Modules
- NumPy
- matplotlib
- MDAnalysis
- scikit-learn

## Usage

All of the required settings are defined in /src/AtomTracker/config.py. Variables with `*` are required.

First, define the directory paths from which the simulation trajectories and corresponding .gro/.pdb files are fetched.

| Variable | Description | Default value |
| --- | - | - |
| `traj_paths` | List of paths from which to fetch .xtc files |  - |
| `gro_paths` | List of paths from which to fetch corresponding .gro/.pdb files | - |
| `save_paths` | List of paths to which the corresponding results are saved | - |


Next, define all settings needed for the calculations of the density maps.

| Variable | Description | Default value |
| --- | - | - |
| `map_selections`* | List of MDAnalysis selection strings for which density maps are calculated |  - |
| `other_selections` | List of MDAnalysis selection strings for which coordinates are saved at every frame | None |
| `reference_structure` | Path to .gro/.pdb file which will be used as reference structure | None |
| `alignment_selection` | MDAnalysis selection string that is used to align each frame to `reference_structure` | None |
| `fit_structures` | Whether to use `alignment_selection` to fit every frame to `reference_structure` | False |
| `centering_selection`* | MDAnalysis string of structure whose center of mass is shifted to the origin on each frame | - |
| `R_min` |  |  |
| `R_max` |  |  |
| `n_R` |  |  |
| `n_theta` |  |  |
| `n_z` |  |  |
| `normalization` |  |  |
| `skip` |  |  |








