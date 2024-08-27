import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.rcParams["font.family"] = "serif"

def density_movie(density_interpolations, function_values, R_grids, theta_grids,
                  vmax=1e-5, save_file="./movie.mp4", interval=10, title="densities", z_bins=None):


    if not z_bins: # If z_bins is None, use approximately equal partition
        z_bins = [density_interpolations.shape[-1]//3, 2*density_interpolations.shape[-1]//3]

    r, th = np.meshgrid(R_grids, theta_grids)

    fig = plt.figure(figsize=(6, 2), dpi=200)
    for i in range(1, 4):
        ax = fig.add_subplot(1,3,i, projection="polar")
        ax.grid()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["polar"].set_color("k")
    
    fig.subplots_adjust(top=0.8)

    frame_titles = [f"Function value = {int(round(val, 0))}" for val in function_values.flatten()]
    max_length = np.max([len(title) for title in frame_titles])
    frame_titles = [s.ljust(max_length) for s in frame_titles]
    fig.text(0.6, 0.95, title, ha="center", va="center", fontsize=12, fontweight="normal")    
    suptitle = fig.suptitle("", fontsize=12, x=0.015, y=0.98, ha="left", color="k")
    bbox_props = dict(boxstyle="round, pad=0.2", fc="lightgray", ec="k", lw=1)
    suptitle.set_bbox(bbox_props)

    def update(frame):
        print(f"At frame {frame + 1} / {len(density_interpolations)}", end="\r")

        for ax in fig.axes:
            ax.clear()
            ax.grid()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines["polar"].set_color("k")

        interp_low = density_interpolations[frame,:,:,:z_bins[0]].mean(axis=2).T
        interp_middle = density_interpolations[frame,:,:,z_bins[0]:z_bins[1]].mean(axis=2).T
        interp_upper = density_interpolations[frame,:,:,z_bins[1]:].mean(axis=2).T


        fig.axes[0].contourf(th, r, interp_low, cmap="Reds", vmin=0, vmax=vmax)
        fig.axes[1].contourf(th, r, interp_middle, cmap="Reds", vmin=0, vmax=vmax)
        fig.axes[2].contourf(th, r, interp_upper, cmap="Reds", vmin=0, vmax=vmax)
        
        suptitle.set_text(frame_titles[frame])
        return suptitle, fig.axes

    ani = FuncAnimation(fig, update, frames=len(density_interpolations), interval=interval)
    fig.tight_layout()
    ani.save(save_file, writer="imagemagick")
    return ani