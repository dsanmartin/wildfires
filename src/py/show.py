import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot(x, y, u, v, w, T, Y, u_lims, v_lims, w_lims, T_lims, Y_lims, label_coor, axis):
    variables = [u, v, w, T, Y]
    variable_name = ['u', 'v', 'w', 'T', 'Y']
    units = ['m/s', 'm/s', 'm/s', 'K', '\%']
    lims = [u_lims, v_lims, w_lims, T_lims, Y_lims]
    cmaps = ['viridis', 'viridis', 'viridis', 'jet', 'Oranges']
    plots = [None, None, None, None, None]
    if 'z' in axis:
        figsize = (12, 12)
        n_rows, n_cols = 5, 1
        sharex, sharey = True, False
    else:
        figsize = (18, 4)
        n_rows, n_cols = 1, 5
        sharex, sharey = False, True
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=sharex, sharey=sharey)
    h_var, v_var = axis[0], axis[1]
    if sharex:
        for i in range(n_rows):
            axes[i].set_ylabel(r'${}$'.format(v_var))
        axes[-1].set_xlabel(r'${}$'.format(h_var))
    if sharey:
        for j in range(n_cols):
            axes[j].set_xlabel(r'${}$'.format(h_var))
        axes[0].set_ylabel(r'${}$'.format(v_var))

    if 'z' in axis:
        for i in range(n_rows):
            axes[i].set_ylim(0, 20)

    for i in range(len(variables)):
        plots[i] = axes[i].contourf(x, y, variables[i], cmap=cmaps[i], vmin=lims[i][0], vmax=lims[i][1])
        m = plt.cm.ScalarMappable(cmap=cmaps[i])
        m.set_clim(lims[i][0], lims[i][1])
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(m, ax=axes[i], cax=cax, label=r'${}$'.format(units[i]))
        axes[i].set_title(r'${}{}$'.format(variable_name[i], label_coor))

    plt.tight_layout()
    plt.show()


if len(sys.argv) < 2:
    print("Usage: python test.py <plot>")
    sys.exit(1)

# Read the data
data_dir = 'data/output/'
# variable = sys.argv[1]
# filename = variable + '.csv.{}'

# print(data[:, 3].min(), data[:, 3].max())
# # Reshape data
Nx, Ny, Nz, Nt = 128, 128, 64, 5
Nz_Y = 4
# if variable == "Y":
#     Nz = 4
#     cmap = 'Oranges'
# else:
#     cmap = 'jet'

# data = data.reshape((Ny, Nx, Nz), order='C')
# T = np.zeros((Ny, Nx, Nz))
# for k in range(Nz):
#     for j in range(Ny):
#         for i in range(Nx):
#             T[j, i, k] = data[k + Nz * (j + Ny * i), 3]
u = np.zeros((Nt, Ny, Nx, Nz))
v = np.zeros((Nt, Ny, Nx, Nz))
w = np.zeros((Nt, Ny, Nx, Nz))
T = np.zeros((Nt, Ny, Nx, Nz))
Y = np.zeros((Nt, Ny, Nx, Nz))
for n in range(Nt):
    data_u = np.loadtxt(data_dir + 'u.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_v = np.loadtxt(data_dir + 'v.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_w = np.loadtxt(data_dir + 'w.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_T = np.loadtxt(data_dir + 'T.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_Y = np.loadtxt(data_dir + 'Y.csv.{}'.format(n), delimiter=',', skiprows=1)
    u[n, :, :, :] = data_u[:, 3].reshape((Ny, Nx, Nz), order='F')
    v[n, :, :, :] = data_v[:, 3].reshape((Ny, Nx, Nz), order='F')
    w[n, :, :, :] = data_w[:, 3].reshape((Ny, Nx, Nz), order='F')
    T[n, :, :, :] = data_T[:, 3].reshape((Ny, Nx, Nz), order='F')
    Y[n, :, :, :Nz_Y] = data_Y[:, 3].reshape((Ny, Nx, Nz_Y), order='F')

x = data_u[:, 0].reshape((Ny, Nx, Nz), order='F')
y = data_u[:, 1].reshape((Ny, Nx, Nz), order='F')
z = data_u[:, 2].reshape((Ny, Nx, Nz), order='F')

# T1 = T[:, 0, :]
# T2 = T[:, -1, :]
# print(np.linalg.norm(T1 - T2))

# print(T.min(), T.max())
# Plot the data
axes = sys.argv[1]#'xz'
k = 0
j = Ny // 2
i = Nx // 2
n1 = 0
n2 = -1
# vmin = R.min()
# vmax = R.max()
u_lims = [u.min(), u.max()]
v_lims = [v.min(), v.max()]
w_lims = [w.min(), w.max()]
T_lims = [T.min(), T.max()]
Y_lims = [Y.min(), Y.max()]

# print(vmin, vmax)
# idxs = np.where(T == 0)
# idxsy, idxsx, idxsz = idxs[1], idxs[2], idxs[3]

# print(len(idxsx), len(idxsy), len(idxsz))

# for i in zip(idxsx, idxsy, idxsz):
#     print(i)
# 3D scatter plot of idxs
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(idxsx, idxsy, idxsz, c='r', marker='o')
# ax.set_xlabel(r'$x$')
# ax.set_ylabel(r'$y$')
# ax.set_zlabel(r'$z$')
# plt.show()


# # plt.imshow(T[:, :, k], cmap='jet', origin='lower')
# # plt.contourf(X[:, :, k], Y[:, :, k], np.abs(T[n1, :, :, k] - T[n2, :, :, k]), cmap='viridis')
# plt.figure(figsize=(12, 4))

# if plot == 'xy':
#     plt.subplot(121)
#     plt.contourf(X[:, :, k], Y[:, :, k], T[n1, :, :, k], cmap='jet', vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.subplot(122)
#     plt.contourf(X[:, :, k], Y[:, :, k], T[n2, :, :, k], cmap='jet', vmin=vmin, vmax=vmax)
#     plt.colorbar()
# elif plot == 'xz':
#     plt.subplot(121)
#     plt.contourf(X[:, j, :], Z[:, j, :], T[n1, :, j, :], cmap='jet', vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.subplot(122)
#     plt.contourf(X[:, j, :], Z[:, j, :], T[n2, :, j, :], cmap='jet', vmin=vmin, vmax=vmax)
#     plt.colorbar()
# elif plot == 'yz':
#     plt.subplot(121)
#     plt.contourf(Y[i, :, :], Z[i, :, :], T[n1, i, :, :], cmap='jet', vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.subplot(122)
#     plt.contourf(Y[i, :, :], Z[i, :, :], T[n2, i, :, :], cmap='jet', vmin=vmin, vmax=vmax)
#     plt.colorbar()

# plt.tight_layout()
# plt.show()

# fig = plt.figure(figsize=(12, 4))
# for n in range(Nt):
#     plt.figure()
#     ax = plt.gca()
#     if plot == 'xy':
#         cf = plt.contourf(X[:, :, k], Y[:, :, k], R[n, :, :, k], cmap=cmap, vmin=vmin, vmax=vmax)
#         plt.xlabel(r'$x$')
#         plt.ylabel(r'$y$')
#         label_coor = "(x, y, {})".format(round(Z[:,:,k].min(), 2))
#     elif plot == 'xz':
#         cf = plt.contourf(X[:, j, :], Z[:, j, :], R[n, :, j, :], cmap=cmap, vmin=vmin, vmax=vmax)
#         plt.xlabel(r'$x$')
#         plt.ylabel(r'$z$')
#         label_coor = "(x, {}, z)".format(round(Y[:,j,:].min(), 2))
#     elif plot == 'yz':
#         cf = plt.contourf(Y[i, :, :], Z[i, :, :], R[n, i, :, :], cmap=cmap, vmin=vmin, vmax=vmax)
#         plt.xlabel(r'$y$')
#         plt.ylabel(r'$z$')
#         label_coor = "({}, y, z)".format(round(X[i,:,:].min(), 2))
#     if "z" in plot:
#         plt.ylim(0, 20)
    
#     # Force colorbar to show full range using vmin and vmax
#     m = plt.cm.ScalarMappable(cmap=cmap)
#     # m.set_array(R[0])
#     m.set_clim(vmin, vmax)
    
#     # Set aspect ratio according to plot limits
#     # plt.gca().set_aspect('equal', adjustable='box')
#     # Fit colorbar to plot height
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.1)
#     plt.colorbar(m, ax=ax, cax=cax, label=r'${}{}$'.format(variable, label_coor))
#     # colorbar = plt.colorbar(cf)
#     # plt.clim(vmin, vmax)
#     plt.tight_layout()
#     plt.show()

for n in range(Nt):
    if axes == 'xy':
        label_coor = "(x, y, {})".format(round(z[:,:,k].min(), 2))
        plot(x[:, :, k], y[:, :, k], u[n, :, :, k], v[n, :, :, k], w[n, :, :, k], T[n, :, :, k], Y[n, :, :, k], u_lims, v_lims, w_lims, T_lims, Y_lims, label_coor, axes)
    elif axes == 'xz':
        label_coor = "(x, {}, z)".format(round(y[:,j,:].min(), 2))
        plot(x[:, j, :], z[:, j, :], u[n, :, j, :], v[n, :, j, :], w[n, :, j, :], T[n, :, j, :], Y[n, :, j, :], u_lims, v_lims, w_lims, T_lims, Y_lims, label_coor, axes)
    elif axes == 'yz':
        label_coor = "({}, y, z)".format(round(x[i,:,:].min(), 2))
        plot(y[i, :, :], z[i, :, :], u[n, i, :, :], v[n, i, :, :], w[n, i, :, :], T[n, i, :, :], Y[n, i, :, :], u_lims, v_lims, w_lims, T_lims, Y_lims, label_coor, axes)
