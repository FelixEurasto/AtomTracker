########### DIRECTORY DEFINITIONS ###########

# Define list of paths where the trajectories are. Can use '*' to select recursively.
traj_paths = [
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol-sdpc/a100/epoch*/rep*/solu_memb_centered.xtc",
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol/a100/epoch*/rep*/solu_memb_centered.xtc",
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/sdpc/a100/epoch*/rep*/solu_memb_centered.xtc",
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/popc/a100/epoch*/rep*/solu_memb_centered.xtc"
]

# Define list of gro/pdb file paths in corresponding order with trajectory dirs.
gro_paths = [
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol-sdpc/a100/epoch01/rep01/solu_memb.gro",
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol/a100/epoch01/rep01/solu_memb.gro",
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/sdpc/a100/epoch01/rep01/solu_memb.gro",
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/popc/a100/epoch01/rep01/solu_memb.gro"
]

# Define list of paths to which position data will be saved.
save_paths = [
    "/wrk/eurastof/mapper/data/chol-sdpc/",
    "/wrk/eurastof/mapper/data/chol/",
    "/wrk/eurastof/mapper/data/sdpc/",
    "/wrk/eurastof/mapper/data/popc/"
]

########### POSITION MAP SETTINGS ###########

map_selections = ["resname CHL1", "resname SDPC", "resname POPC"] # Selections for which maps are calculated
other_selections = ["name CA"] # Selections for which coordinates corresponding to maps are calculated. Set to None if not needed.

# Reference structure used for alignment at each timestep of each simulation. If None, each trajectoyr will be aligned with it's first frame.
reference_structure = "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol-sdpc/a100/epoch01/rep01/solu_memb.gro"
alignment_selection = "name CA" # selection used to align each ts to reference structure
centering_selection = "protein" # selection whose COM is centered to origin at each frame

R_min = 0 # Minimum radial distance from origin in xy-plane
R_max = 30 # Maximum radial distance from origin in xy-plane
n_R = 30 # Number of grid points in R-direction
n_z = 30 # Number of grid points in z-direction
n_theta = 30 # Number of grid points in theta-direction
use_com = False # Whether to use the COMs of selections when calculating position maps
# Normalization method. 'within' normalizes each frame's map to sum to 1. 'all' normalizes each frame's map by the total number of atoms in the map selection.
normalization = "within" 

skip = 30 # Number of frames to skip 


########### FUNCTION DEFINITION ###########

# Define function that is calculated for the corresponding position maps.
# If you do not need this, just set function = None below. 
def calculate_function(universe):
    import scipy
    pairs = [(23, 297), (48, 87), (92, 119), (196, 241), (265, 277)]
    coefs = [-14.43, -7.62, 9.11, -6.32, -5.22]
    a100 = 0
    for i, pair in enumerate(pairs):
        resid1 = pair[0]
        resid2 = pair[1]
        coef = coefs[i]
        p1 = universe.select_atoms(f"name CA and resid {resid1}")
        p2 = universe.select_atoms(f"name CA and resid {resid2}")
        R = scipy.spatial.distance.euclidean(p1.positions.reshape(-1), p2.positions.reshape(-1))
        a100 += coef * R
    a100 += 278.88
    return a100


function = calculate_function # Set function = None if you do not wish to calculate any function

