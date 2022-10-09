from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colorbar
import matplotlib.colors
import numpy as np
import pandas as pd
from pylab import *

maxsite_1d = 8
num_grid = maxsite_1d ** 3
chunksize = 2 * num_grid + 1

for chunk in pd.read_csv(
    "dens_profile_N32_A100_B100.txt", header=0, delimiter="\t", chunksize=chunksize
):
    timestep = int(chunk.iat[0, 0].split("= ")[1])

    print("timestep =", timestep, ", step =", timestep // 1000)
    # print(chunk.iloc[1:])
    old_arr = chunk.iloc[1:].to_numpy()
    x = old_arr[:, 1].astype(float)
    y = old_arr[:, 2].astype(float)
    z = old_arr[:, 3].astype(float)
    itype = old_arr[:, 0].astype(int)
    dens_data = old_arr[:, 4].astype(float)
    # dens_matrix = dens_data.reshape((2, maxsite_1d, maxsite_1d, maxsite_1d))
    # For debug:
    # for itype in range(2):
    #     for i in range(maxsite_1d):
    #         for j in range(maxsite_1d):
    #             for k in range(maxsite_1d):
    #                 if (
    #                     dens_matrix[itype][i][j][k]
    #                     != dens_data[itype * 512 + i * 64 + j * 8 + k]
    #                 ):
    #                     print(False)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(
        x[:num_grid],
        y[:num_grid],
        z[:num_grid],
        c=dens_data[:num_grid],
        cmap="hot",
        marker="s",
        s=50,
    )
    color_map = cm.ScalarMappable(cmap="hot")
    plt.colorbar(color_map)
    ax.set_title("3D Heatmap")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    plt.show()

