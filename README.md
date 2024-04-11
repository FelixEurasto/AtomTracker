# A tool to build density maps of atoms from MD simulations


## Python Modules
- NumPy
- matplotlib
- MDAnalysis
- scikit-learn

## Usage

All of the required settings are defined in /src/AtomTracker/config.py

First, define the directory paths from which the simulation trajectories and corresponding gro/pdb files are fetched.

| Variable | Description | Default value |
| --- | - | - |
| `traj_paths` | List of paths from which to fetch .xtc files |  - |
| `gro_paths` | List of paths from which to fetch corresponding .gro/.pdb files | - |
| `save_paths` | List of paths to which the corresponding results are saved | - |
