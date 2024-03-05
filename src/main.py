from mapper import map_calculations, position_calculations, config, plotting
import glob
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def main():
    for system_ind, traj_path in enumerate(config.traj_paths):
        print(f"Extracting data from trajectories in '{traj_path}'")
        try:
            gro = glob.glob(config.gro_paths[system_ind])[0]
            trajs = glob.glob(traj_path)
            cosmos = mda.Universe(gro, trajs)
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
            cosmos = mda.Universe(gro, traj)
            coordinates = position_calculations.calculate_positions(universe=cosmos,
                                                                    map_selections=config.map_selections,
                                                                    use_com=config.use_com,
                                                                    skip=config.skip,
                                                                    reference_structure=config.reference_structure,
                                                                    other_selections=config.other_selections,
                                                                    centering_selection=config.centering_selection,
                                                                    alignment_selection=config.alignment_selection)
            for sel in config.map_selections:
                X = coordinates[sel]
                traj_sel_maps = map_calculations.create_maps(X=X,
                                                             R_min=config.R_min,
                                                             R_max=config.R_max,
                                                             n_R=config.n_R,
                                                             n_theta=config.n_theta,
                                                             n_z=config.n_z,
                                                             normalization=config.normalization)
                maps[sel].append(traj_sel_maps)
            if config.other_selections:
                for sel in config.other_selections:
                    other_coordinates[sel].append(coordinates[sel])
            
            traj_fvals = [config.calculate_function(cosmos) for ts in cosmos.trajectory[::config.skip]]
            if config.function:
                fvals += traj_fvals
            traj_inds += [traj_ind for _ in range(len(traj_fvals))]
        
        for sel in config.map_selections:
            sel_in_file_name = sel.replace(" ", "_")
            maps[sel] = np.concatenate(maps[sel], axis=0)
            np.save(f"{config.save_paths[system_ind]}{sel_in_file_name}_maps.npy", maps[sel])
            fig = plotting.make_radial_plot(maps=maps[sel], R_min=config.R_min, R_max=config.R_max, vmax=maps[sel].mean(axis=0).max(), title=f"Mean densities for {sel}")
            plt.savefig(f"{config.save_paths[system_ind]}{sel_in_file_name}_densities.png", bbox_inches="tight", dpi=300)
        if config.other_selections:
            for sel in config.other_selections:
                sel_in_file_name = sel.replace(" ", "_")
                other_coordinates[sel] = np.concatenate(other_coordinates[sel], axis=0)
                np.save(f"{config.save_paths[system_ind]}{sel_in_file_name}_coordinates.npy", other_coordinates[sel])
        if config.function:
            np.save(f"{config.save_paths[system_ind]}function_values.npy", np.array(fvals))
        np.save(f"{config.save_paths[system_ind]}traj_inds.npy", np.array(traj_inds)) 


if __name__ == "__main__":
    main()