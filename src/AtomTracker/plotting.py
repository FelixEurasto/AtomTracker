import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
import plotly.graph_objects as go
from scipy.stats import gaussian_kde
from ..AtomTracker import map_calculations
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


def make_radial_plot(maps, R_grids, theta_grids, vmax, title, cmap="Reds", z_bins=None):
    """
    Plots mean position maps for 3 sections along the z-direction
    Parameters:

        vmax -- float: vmax given to 'pcolormesh'
        title -- string: Title given to figure
    Returns:
        fig -- plt.figure: Mean radial position maps
    """
    
    fig = plt.figure(figsize=(9, 3), layout="constrained")
    ax = [fig.add_subplot(1,3,i,projection="polar") for i in [1,2,3]]
    mean_map = maps.mean(axis=0)
    
    if not z_bins: # If z_bins is None, use approximately equal partition
        z_bins = [mean_map.shape[-1]//3, 2*mean_map.shape[-1]//3]
        
    lower_map = maps[:,:,:,:z_bins[0]].mean(axis=0).mean(axis=2).T
    middle_map = maps[:,:,:,z_bins[0]:z_bins[1]].mean(axis=0).mean(axis=2).T
    upper_map = maps[:,:,:,z_bins[1]:].mean(axis=0).mean(axis=2).T
    r, th = np.meshgrid(R_grids, theta_grids)
    
    ax[0].contourf(th, r, lower_map, cmap="Reds", vmax=vmax, vmin=0)
    ax[1].contourf(th, r, middle_map, cmap="Reds", vmax=vmax, vmin=0)
    ax[2].contourf(th, r, upper_map, cmap="Reds", vmax=vmax, vmin=0)

    for a in ax:
        a.set_thetagrids([], [])
        a.set_rticks([])
    ax[0].set_title("Lower")
    ax[1].set_title("Central")
    ax[2].set_title("Upper")
    
    fig.suptitle(title, fontsize=15)

    return fig


def animate(coordinates, densities, points, fvals, threshold=0.4, cmap="hot", title=""):
    
    n_frames = len(densities)
    frames = []
    for i in range(n_frames):
        data = []
        if coordinates is not None:
            if coordinates.shape[0] > 1:
                data.append(go.Scatter3d(
                            visible=True,
                            x=coordinates[i,:,0],
                            y=coordinates[i,:,1],
                            z=coordinates[i,:,2],
                            hoverinfo="none",
                            mode="lines",
                            line=dict(width=12, color="darkgreen"),
                    ))
            else:
                data.append(go.Scatter3d(
                            visible=True,
                            x=coordinates[0,:,0],
                            y=coordinates[0,:,1],
                            z=coordinates[0,:,2],
                            hoverinfo="none",
                            mode="lines",
                            line=dict(width=12, color="darkgreen"),
                    ))
        if np.where(densities[i] >= threshold)[0].size != 0:
            data.append(go.Volume(
                            visible=True,
                            x=points[:,0].flatten(),
                            y=points[:,1].flatten(),
                            z=points[:,2].flatten(),
                            value=densities[i],
                            isomin=threshold,
                            isomax=1,
                            opacity=.7,
                            surface_count=10,
                            colorscale=cmap,
                            showscale=False
                    ))
        frames.append(go.Frame(data=data, name=str(i)))

    fig = go.Figure(
        frames=frames
    )

    for d in fig.frames[0].data:
        fig.add_trace(d)


    def frame_args(duration):
        return {
                "frame": {"duration": duration},
                "mode": "immediate",
                "fromcurrent": True,
                "transition": {"duration": duration, "easing": "linear"},
            }

    sliders = [
                {
                    "pad": {"b": 10, "t": 60},
                    "len": 0.9,
                    "x": 0.1,
                    "y": 0,
                    "steps": [
                        {
                            "args": [[f.name], frame_args(0)],
                            "label": str(int(fvals[k])),
                            "method": "animate",
                        }
                        for k, f in enumerate(fig.frames)
                    ],
                }
            ]
    
    fig.update_layout(
         title=title,
         width=600,
         height=600,
         scene=dict(
                    xaxis=dict(showticklabels=False, showgrid=False, title=""),
                    yaxis=dict(showticklabels=False, showgrid=False, title=""),
                    zaxis=dict(showticklabels=False, showgrid=False, title=""),
                    bgcolor="white",
                    camera_projection_type="orthographic",
                    dragmode="orbit"
                    ),
        template="plotly_white",

        updatemenus = [
            {
                "buttons": [
                    {
                        "args": [None, frame_args(5)],
                        "label": "&#9654;",
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
         ],
        sliders=sliders,
        )
    fig.add_annotation(
        text="Function value",
        x=0.5,
        y=-0.2,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=20)
    )


    return fig

def calculate_kdes(maps, map_points, density_points):
    X = []
    for m in maps:
        m[m < 0] = 0
        kde = gaussian_kde(map_points.T, weights=m)
        density = kde.evaluate(density_points.T)
        X.append(density.reshape(1,-1))
    return np.concatenate(X, axis=0)


def _3d_plot(maps, coordinates, fvals, R_grids, theta_grids, z_grids,
             res=(10,10,10), fval_range=(0,1), n_frames=20, threshold=0.1, title=""):
    x_limits = (-R_grids.max(), R_grids.max())
    y_limits = (-R_grids.max(), R_grids.max())
    z_limits = (z_grids.min(), z_grids.max())

    points = []
    for x in np.linspace(x_limits[0], x_limits[1], res[0]):
        for y in np.linspace(y_limits[0], y_limits[1], res[1]):
            for z in np.linspace(z_limits[0], z_limits[1], res[2]):
                points.append([x,y,z])
    points = np.array(points)

    grid_coordinates = []
    for r in R_grids[:-1]:
        for t in theta_grids:
            for z in z_grids:
                grid_coordinates.append([r, t, z])
    grid_coordinates = np.array(grid_coordinates)
    XYZ = map_calculations.cylindrical_to_cartesian(grid_coordinates[:,0], grid_coordinates[:,1], grid_coordinates[:,2])
    X = StandardScaler().fit_transform(maps.reshape(len(maps), -1))
    y_scaler = StandardScaler().fit(fvals.reshape(-1,1))
    y = y_scaler.transform(fvals.reshape(-1,1))
    train_X, test_X, train_y, test_y = train_test_split(X, y, test_size=.33)

    pls_model = PLSRegression(n_components=1).fit(train_X, train_y)
    model_score = pls_model.score(test_X, test_y)
    fvals_to_interpolate = y_scaler.transform(np.linspace(fval_range[0], fval_range[1], n_frames).reshape(n_frames, 1))
    density_interpolations = pls_model.inverse_transform(fvals_to_interpolate).reshape((len(fvals_to_interpolate),
                                                                                        R_grids.shape[0] - 1,
                                                                                        theta_grids.shape[0],
                                                                                        z_grids.shape[0]))
    density_kde = calculate_kdes(density_interpolations.reshape(len(density_interpolations), -1), XYZ, points)
    density_kde /= density_kde.max()
    
    
    fig = animate(coordinates=coordinates, densities=density_kde,
                  points=points, threshold=threshold,
                  fvals=y_scaler.inverse_transform(fvals_to_interpolate.reshape(1, n_frames)).flatten(),
                  title=title + "<br>" + f"PLS model R2: {model_score:.3}")
    
    return fig
    

