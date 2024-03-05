import numpy as np
import matplotlib.pyplot as plt


def make_radial_plot(maps, R_min, R_max, vmax, title):
    """
    Plots mean position maps for 3 sections along the z-direction
    Parameters:
        lower_maps -- np.ndarray: Maps from lowest section as returned by 'three_maps'
        middle_maps -- np.ndarray: Maps from middle section as returned by 'three_maps'
        upper_maps -- np.ndarray: Maps from upper section as returned by 'three_maps'
        R_max -- float: Maximum radial distance from origin in xy-plane
        vmax -- float: vmax given to 'pcolormesh'
        title -- string: Title given to figure
    Returns:
        fig -- plt.figure: Mean radial position maps
    """
    
    fig = plt.figure(figsize=(9, 3), layout="constrained")
    ax = [fig.add_subplot(1,3,i,projection="polar") for i in [1,2,3]]

    z_bins = np.linspace(0, maps.shape[-1], 4, dtype=int)
    mean_map = maps.mean(axis=0)
    lower_map = mean_map[:,:,0:z_bins[1]].mean(axis=2).T
    middle_map = mean_map[:,:,z_bins[1]:z_bins[2]].mean(axis=2).T
    upper_map = mean_map[:,:,z_bins[2]:].mean(axis=2).T
    n_R, n_theta = maps.shape[1], maps.shape[2]
    rad = np.linspace(R_min, R_max, n_R)
    azm = np.linspace(0, 2 * np.pi, n_theta)
    r, th = np.meshgrid(rad, azm)
    ax[0].pcolormesh(th, r, lower_map, cmap="Reds", vmin=0, vmax=vmax)
    ax[1].pcolormesh(th, r, middle_map, cmap="Reds", vmin=0, vmax=vmax)
    ax[2].pcolormesh(th, r, upper_map, cmap="Reds", vmin=0, vmax=vmax)

    ax[0].set_title("Lower")
    ax[1].set_title("Central")
    ax[2].set_title("Upper")
    
    fig.suptitle(title, fontsize=15)
    return fig