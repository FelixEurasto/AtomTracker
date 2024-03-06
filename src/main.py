from mapper import map_calculations, position_calculations, config, plotting
import glob
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression

def main():
    for system_ind, traj_path in enumerate(config.traj_paths):
        print(f"Extracting data from trajectories in '{traj_path}' ...")
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
            if len(trajs) > 1:
                print(f"At trajectory {traj_ind + 1} / {len(trajs)} ...", end="\r")
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
            try:
                maps[sel] = np.concatenate(maps[sel], axis=0)
                np.save(f"{config.save_paths[system_ind]}{sel_in_file_name}_maps.npy", maps[sel])
                fig = plotting.make_radial_plot(maps=maps[sel], R_min=config.R_min, R_max=config.R_max, vmax=maps[sel].mean(axis=0).max(), title=f"Mean densities for {sel}")
                fig.savefig(f"{config.save_paths[system_ind]}{sel_in_file_name}_densities.png", bbox_inches="tight", dpi=300)
            except:
                print(f"Selection {sel} not found in this system!")
                continue

        if config.function:
            print("Fitting PLS models ...")
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
                    s = ax[1].scatter(pls_2d_projected[:,0], pls_2d_projected[:,1], c=test_fvals, edgecolors="k")
                    ax[0].set_title(f"1 Component PLS (R^2: {pls_1d_score:.3})")
                    ax[1].set_title(f"2 Component PLS (R^2: {pls_2d_score:.3})")
                    fig.colorbar(mappable=s, label="Function value")
                    fig.savefig(f"{config.save_paths[system_ind]}{sel_in_file_name}_pls_projections.png", bbox_inches="tight", dpi=300)
                else:
                    continue
            # Combined PLS model
            sel_in_file_name = "_".join([sel.replace(" ", "_") for sel in config.map_selections if np.array(maps[sel]).size != 0])
            combined_maps = np.concatenate([maps[sel] for sel in config.map_selections if np.array(maps[sel]).size != 0], axis=1)
            combined_maps = combined_maps.reshape(combined_maps.shape[0], combined_maps.shape[1]*combined_maps.shape[2]*combined_maps.shape[3])
            train_maps, test_maps, train_fvals, test_fvals = train_test_split(combined_maps, fvals, test_size=.25)

            pls_1d = PLSRegression(n_components=1).fit(train_maps, train_fvals)
            pls_1d_score = pls_1d.score(test_maps, test_fvals)
            pls_1d_projected = pls_1d.transform(test_maps)

            pls_2d = PLSRegression(n_components=2).fit(train_maps, train_fvals)
            pls_2d_score = pls_2d.score(test_maps, test_fvals)
            pls_2d_projected = pls_2d.transform(test_maps)

            fig, ax = plt.subplots(1,2,figsize=(6,3),layout="constrained")
            ax[0].hist(pls_1d_projected, bins=30)
            s = ax[1].scatter(pls_2d_projected[:,0], pls_2d_projected[:,1], c=test_fvals, edgecolors="k")
            fig.colorbar(mappable=s, label="Function value")
            ax[0].set_title(f"1 Component PLS (R^2: {pls_1d_score:.3})")
            ax[1].set_title(f"2 Component PLS (R^2: {pls_2d_score:.3})")
            fig.savefig(f"{config.save_paths[system_ind]}{sel_in_file_name}_pls_projections.png", bbox_inches="tight", dpi=300)
        
        print(f"Saving arrays into {config.save_paths[system_ind]} ...")
        if config.other_selections:
            for sel in config.other_selections:
                sel_in_file_name = sel.replace(" ", "_")
                other_coordinates[sel] = np.concatenate(other_coordinates[sel], axis=0)
                np.save(f"{config.save_paths[system_ind]}{sel_in_file_name}_coordinates.npy", other_coordinates[sel])
        if config.function:
            np.save(f"{config.save_paths[system_ind]}function_values.npy", fvals)
        np.save(f"{config.save_paths[system_ind]}traj_inds.npy", traj_inds) 


if __name__ == "__main__":
    main()