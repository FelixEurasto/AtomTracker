########### DIRECTORY DEFINITIONS ###########

# Define list of paths where the trajectories are. Can use '*' to select recursively.
traj_paths = [
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol-sdpc/a100/epoch*/rep*/mdrun.xtc",
    #"/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol/a100/epoch*/rep*/mdrun.xtc",
    #"/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/sdpc/a100/epoch*/rep*/mdrun.xtc",
    #"/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/popc/a100/epoch*/rep*/mdrun.xtc"
]

# Define list of gro/pdb file paths in corresponding order with trajectory dirs.
gro_paths = [
    "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol-sdpc/a100/epoch01/rep01/mdrun.gro",
    #"/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol/a100/epoch01/rep01/mdrun.gro",
    #"/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/sdpc/a100/epoch01/rep01/mdrun.gro",
    #"/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/popc/a100/epoch01/rep01/mdrun.gro"
]

# Define list of paths to which position data will be saved.
save_paths = [
    "/wrk/eurastof/AtomTracker/data/chol-sdpc-full-system/",
    #"/wrk/eurastof/AtomTracker/data/chol-full-system/",
    #"/wrk/eurastof/AtomTracker/data/sdpc-full-system/",
    #"/wrk/eurastof/AtomTracker/data/popc-full-system/"
]

########### DENSITY MAP SETTINGS ###########

map_selections = ["resname TIP3", "resname SOD", "resname CLA"] # Selections for which maps are calculated
other_selections = ["name CA"] # Selections for which coordinates corresponding to maps are calculated. Set to None if not needed.

# Reference structure used for alignment at each timestep of each simulation. If None, each trajectory will be aligned with its first frame.
reference_structure = "/wrk/eurastof/binding_spots_project/gpcr_sampling/b2ar-fst/chol-sdpc/a100/epoch01/rep01/solu_memb.gro"
alignment_selection = "name CA and not resid 199-235 and not resid 298-311 and not resid 265-275 and not resid 135-165 and not resid 105-115 and not resid 60-75 and not resid 29-36" # selection used to align each ts to reference structure
fit_structures = True # Whether to align frames to reference_structure, using alignment_selection.
centering_selection = "protein and not name H and not type H" # selection whose COM is centered to origin at each frame

R_min = 0 # Minimum radial distance from origin in xy-plane
R_max = 30 # Maximum radial distance from origin in xy-plane
n_R = 30 # Number of grid points in R-direction
n_theta = 30 # Number of grid points in theta-direction
n_z = 30 # Number of grid points in z-direction
# Normalization method. 'within' normalizes each frame's map by the number of atoms within R_max.
# 'all' normalizes each frame's map by the total number of atoms in the map selection.
# If None, no normalization based on atom count is performed
normalization = "all"
skip = 50 # Number of frames to skip 

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

########### PLOTTING & PLS MOVIE CREATION ###########

# iterable containing z-dimension cut-offs used for partitioning the density maps for the plots and movies.
# For example, z_bins = [10, 20] will use the density in grids below index 10 as the lower section,
# the density between indices 10 and 20 as the middle section, and the density above index 20 as the upper section.
# If None, the z-direction will be partitioned into 3 (approximately) equal sections.
z_bins = [13, 19]

create_movies = True # Whether to create PLS movies
movie_n_frames = 500 # Number of frames in created movies
movie_length = 5 # Length of created movies in seconds
function_min_and_max = (-50, 50) # Range of function values at which atom densities are interpolated