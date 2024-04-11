# A tool to build density maps of atoms from MD simulations


## Required Python Modules
- NumPy
- matplotlib
- MDAnalysis
- scikit-learn

## Usage

All of the required settings are defined in /src/AtomTracker/config.py

First, define the directory paths from which the simulation trajectories and corresponding .gro/.pdb files are fetched.

| Variable | Description | Default value |
| --- | - | - |
| `traj_paths` | List of paths from which to fetch .xtc files |  - |
| `gro_paths` | List of paths from which to fetch corresponding .gro/.pdb files | - |
| `save_paths` | List of paths to which the corresponding results are saved | - |


Next, define all settings needed for the calculations of the density maps.

| Variable | Description | Default value |
| --- | - | - |
| `map_selections` | List of paths from which to fetch .xtc files |  - |
| `other_selections` | List of paths from which to fetch corresponding .gro/.pdb files | - |
| `reference_structure` | List of paths to which the corresponding results are saved | - |
| `alignment_selection` |  |  |
| `fit_structures` |  |  |
| `centering_selection` |  |  |
| `` |  |  |







