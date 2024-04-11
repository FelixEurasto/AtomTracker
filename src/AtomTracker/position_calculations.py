import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align


def calculate_positions(universe, map_selections, skip=1, fit_structures=False,
                        reference_structure=None, other_selections=None,
                        centering_selection=None, alignment_selection=None):
    coordinates = dict([(sel, []) for sel in map_selections])
    atoms = dict([(sel, universe.select_atoms(sel)) for sel in map_selections])
    
    if other_selections:
        for sel in other_selections:
            coordinates[sel] = []
            atoms[sel] = universe.select_atoms(sel)
    
    if reference_structure:
        ref_struct = mda.Universe(reference_structure)
    else:
        ref_struct = universe.copy()
    ref_struct.atoms.positions -= ref_struct.select_atoms(centering_selection).center_of_mass()

    for ts in universe.trajectory[::skip]:
        if fit_structures:
            align.alignto(universe, ref_struct, select=alignment_selection)
        else:
            universe.atoms.positions -= universe.select_atoms(centering_selection).center_of_mass()
        for sel, sel_atoms in atoms.items():
            coordinates[sel].append(sel_atoms.positions)
    for sel, coords in coordinates.items():
        coordinates[sel] = np.array(coords)

    return coordinates