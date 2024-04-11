import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"

def make_radial_plot(maps, R_min, R_max, vmax, title, cmap="Reds"):
    """
    Plots mean position maps for 3 sections along the z-direction
    Parameters:
        lower_maps -- np.ndarray: Maps from lowest section as returned by 'three_maps'
        middle_maps -- np.ndarray: Maps from middle section as returned by 'three_maps'
        upper_maps -- np.ndarray: Maps from upper section as returned by 'three_maps'
        R_max -- float: Maximum radial distance from origin in the xy-plane
        vmax -- float: vmax given to 'pcolormesh'
        title -- string: Title given to figure
    Returns:
        fig -- plt.figure: Mean radial position maps
    """
    
    fig = plt.figure(figsize=(9, 3), layout="constrained")
    ax = [fig.add_subplot(1,3,i,projection="polar") for i in [1,2,3]]

    total_density = maps.sum()
    z_low_ind = [i for i in range(maps.shape[-1]) if maps[:,:,:,:i].sum() >= total_density/3][0]
    z_middle_ind = [i for i in range(maps.shape[-1]) if maps[:,:,:,:i].sum() >= 2*total_density/3][0]
    z_bins = [z_low_ind, z_middle_ind]

    mean_map = maps.mean(axis=0)
    lower_map = mean_map[:,:,0:z_bins[0]].mean(axis=2).T
    middle_map = mean_map[:,:,z_bins[0]:z_bins[1]].mean(axis=2).T
    upper_map = mean_map[:,:,z_bins[1]:].mean(axis=2).T

    n_R, n_theta = maps.shape[1], maps.shape[2]
    rad = np.linspace(R_min, R_max, n_R)
    azm = np.linspace(0, 2 * np.pi, n_theta)
    r, th = np.meshgrid(rad, azm)

    ax[0].contourf(th, r, lower_map, cmap=cmap, vmin=0, vmax=vmax)
    ax[1].contourf(th, r, middle_map, cmap=cmap, vmin=0, vmax=vmax)
    ax[2].contourf(th, r, upper_map, cmap=cmap, vmin=0, vmax=vmax)

    for a in ax:
        a.set_thetagrids([], [])
        a.set_rticks([])
    ax[0].set_title("Lower")
    ax[1].set_title("Central")
    ax[2].set_title("Upper")
    
    fig.suptitle(title, fontsize=15)
    return fig