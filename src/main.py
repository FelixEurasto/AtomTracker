from AtomTracker import map_calculations, position_calculations, config, plotting
import glob
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
import os

def main():
    for system_ind, traj_path in enumerate(config.traj_paths):
        save_path = config.save_paths[system_ind]
        if not os.path.isdir(save_path):
            print(f"Directory {save_path} not found, creating it now...")
            os.mkdir(save_path)
        print(f"Extracting data from trajectories in '{traj_path}' ...")
        try:
            gro = glob.glob(config.gro_paths[system_ind])[0]
            trajs = glob.glob(traj_path)
            cosmos = mda.Universe(gro, trajs)
            print(f"Calculating z-coordinate info ...")
            if config.reference_structure:
                ref_struct = mda.Universe(config.reference_structure)
            else:
                ref_struct = cosmos.copy()
            ref_struct.atoms.positions -= ref_struct.select_atoms(config.centering_selection).center_of_mass()
            selection = " or ".join(config.map_selections)
            atoms = cosmos.select_atoms(selection)
            z_min_max = []
            for ts in cosmos.trajectory[::int(cosmos.trajectory.n_frames / 100)]:
                if config.fit_structures:
                    align.alignto(cosmos, ref_struct, select=config.alignment_selection)
                else:
                    cosmos.atoms.positions -= cosmos.select_atoms(config.centering_selection).center_of_mass()
                z_min_max.append(atoms.positions[:,2].min())
                z_min_max.append(atoms.positions[:,2].max())
            z_min = np.min(z_min_max)
            z_max = np.max(z_min_max)
            R_grids, theta_grids, z_grids, volumes = map_calculations.make_cylindrical_grid(R_min=config.R_min, R_max=config.R_max,
                                                                                            z_min=z_min, z_max=z_max,
                                                                                            n_R=config.n_R + 1, n_theta=config.n_theta,
                                                                                            n_z=config.n_z)
        except:
            print(f"Unable to load Universe from provided gro & xtc(s)!")
            continue
        
        maps = dict([(sel, []) for sel in config.map_selections])
        if config.function:
            fvals = []
        traj_inds = []
        if config.other_selections:
            other_coordinates = dict([(sel, []) for sel in config.other_selections])
        
        for traj_ind, traj in enumerate(trajs):
            if len(trajs) > 1:
                print(f"At trajectory {traj_ind + 1} / {len(trajs)} ...", end="\r")
            cosmos = mda.Universe(gro, traj)
            coordinates = position_calculations.calculate_positions(universe=cosmos,
                                                                    map_selections=config.map_selections,
                                                                    skip=config.skip, fit_structures=config.fit_structures,
                                                                    reference_structure=config.reference_structure,
                                                                    other_selections=config.other_selections,
                                                                    centering_selection=config.centering_selection,
                                                                    alignment_selection=config.alignment_selection)
            for sel in config.map_selections:
                X = coordinates[sel]
                traj_sel_maps = map_calculations.create_maps(X=X,
                                                             R_grids=R_grids,
                                                             theta_grids=theta_grids,
                                                             z_grids=z_grids,
                                                             volumes=volumes,
                                                             normalization=config.normalization)
                if traj_sel_maps is not None: # Check if map for this selection exists
                    maps[sel].append(traj_sel_maps)
            if config.other_selections:
                for sel in config.other_selections:
                    other_coordinates[sel].append(coordinates[sel])
            
            if config.function:
                traj_fvals = [config.calculate_function(cosmos) for ts in cosmos.trajectory[::config.skip]]
                fvals += traj_fvals
            traj_inds += [traj_ind for _ in range(len(traj_fvals))]
        
        if config.function:
            fvals = np.array(fvals)
        traj_inds = np.array(traj_inds)

        for sel in config.map_selections:
            sel_in_file_name = sel.replace(" ", "_")
            if len(maps[sel]) > 0:
                maps[sel] = np.concatenate(maps[sel], axis=0)
                np.save(f"{save_path}{sel_in_file_name}_maps.npy", maps[sel])
                fig = plotting.make_radial_plot(maps=maps[sel], R_min=config.R_min, R_max=config.R_max, vmax=maps[sel].mean(axis=0).max()/10, title=f"Mean densities for {sel}")
                fig.savefig(f"{save_path}{sel_in_file_name}_densities.png", bbox_inches="tight", dpi=300)
    

        if config.function:
            print("\nFitting PLS models ...")
            # PLS models for every selection in 'map_selections'
            for sel in config.map_selections:
                sel_in_file_name = sel.replace(" ", "_")
                if np.array(maps[sel]).size != 0:
                    train_maps, test_maps, train_fvals, test_fvals = train_test_split(maps[sel].reshape(maps[sel].shape[0],
                                                                                                        maps[sel].shape[1]*maps[sel].shape[2]*maps[sel].shape[3]),
                                                                                        fvals, test_size=.25)
                    pls_1d = PLSRegression(n_components=1).fit(train_maps, train_fvals)
                    pls_1d_score = pls_1d.score(test_maps, test_fvals)
                    pls_1d_projected = pls_1d.transform(test_maps)

                    pls_2d = PLSRegression(n_components=2).fit(train_maps, train_fvals)
                    pls_2d_score = pls_2d.score(test_maps, test_fvals)
                    pls_2d_projected = pls_2d.transform(test_maps)

                    fig, ax = plt.subplots(1,2,figsize=(6,3), layout="constrained")
                    ax[0].hist(pls_1d_projected, bins=30)
                    s = ax[1].scatter(pls_2d_projected[:,0], pls_2d_projected[:,1], c=test_fvals, s=5)
                    ax[0].set_title(f"1 Component PLS (R^2: {pls_1d_score:.3})")
                    ax[1].set_title(f"2 Component PLS (R^2: {pls_2d_score:.3})")
                    fig.colorbar(mappable=s, label="Function value")
                    fig.savefig(f"{save_path}{sel_in_file_name}_pls_projections.png", bbox_inches="tight", dpi=300)
                else:
                    continue
            # Combined PLS model
            sel_in_file_name = "_".join([sel.replace(" ", "_") for sel in config.map_selections if np.array(maps[sel]).size != 0])
            combined_maps = np.concatenate([maps[sel] for sel in config.map_selections if np.array(maps[sel]).size != 0], axis=1)
            combined_maps = combined_maps.reshape(combined_maps.shape[0], combined_maps.shape[1]*combined_maps.shape[2]*combined_maps.shape[3])
            train_maps, test_maps, train_fvals, test_fvals = train_test_split(combined_maps, fvals, test_size=.5)

            pls_1d = PLSRegression(n_components=1).fit(train_maps, train_fvals)
            pls_1d_score = pls_1d.score(test_maps, test_fvals)
            pls_1d_projected = pls_1d.transform(test_maps)

            pls_2d = PLSRegression(n_components=2).fit(train_maps, train_fvals)
            pls_2d_score = pls_2d.score(test_maps, test_fvals)
            pls_2d_projected = pls_2d.transform(test_maps)

            fig, ax = plt.subplots(1,2,figsize=(8,3),layout="constrained")
            ax[0].hist(pls_1d_projected, bins=30)
            s = ax[1].scatter(pls_2d_projected[:,0], pls_2d_projected[:,1], c=test_fvals, s=5)
            fig.colorbar(mappable=s, label="Function value")
            ax[0].set_title(f"1 Component PLS (R^2: {pls_1d_score:.3})")
            ax[1].set_title(f"2 Component PLS (R^2: {pls_2d_score:.3})")
            fig.savefig(f"{save_path}{sel_in_file_name}_pls_projections.png", bbox_inches="tight", dpi=300)
        
        print(f"Saving arrays into {save_path} ...")
        if config.other_selections:
            for sel in config.other_selections:
                sel_in_file_name = sel.replace(" ", "_")
                other_coordinates[sel] = np.concatenate(other_coordinates[sel], axis=0)
                np.save(f"{save_path}{sel_in_file_name}_coordinates.npy", other_coordinates[sel])
        if config.function:
            np.save(f"{save_path}function_values.npy", fvals)
        np.save(f"{save_path}traj_inds.npy", traj_inds) 


if __name__ == "__main__":
    main()