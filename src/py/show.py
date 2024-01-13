import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot(domain, variables, lims, label_coor, axis):
    variable_name = ['u', 'v', 'w', 'T', 'Y']
    units = ['m/s', 'm/s', 'm/s', 'K', '\%']
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
    x, y = domain
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

Nx, Ny, Nz, Nt = 128, 128, 64, 5
Nz_Y = 4

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

axes = sys.argv[1] 
k = 0
j = Ny // 2
i = Nx // 2

u_lims = [u.min(), u.max()]
v_lims = [v.min(), v.max()]
w_lims = [w.min(), w.max()]
T_lims = [T.min(), T.max()]
Y_lims = [Y.min(), Y.max()]
lims = [u_lims, v_lims, w_lims, T_lims, Y_lims]

for n in range(Nt):
    if axes == 'xy':
        label_coor = "(x, y, {})".format(round(z[:,:,k].min(), 2))
        variables = [u[n, :, :, k], v[n, :, :, k], w[n, :, :, k], T[n, :, :, k], Y[n, :, :, k]]
        domain = [x[:, :, k], y[:, :, k]]
    elif axes == 'xz':
        label_coor = "(x, {}, z)".format(round(y[:,j,:].min(), 2))
        variables = [u[n, :, j, :], w[n, :, j, :], T[n, :, j, :], Y[n, :, j, :]]
        domain = [x[:, j, :], z[:, j, :]]
    elif axes == 'yz':
        label_coor = "({}, y, z)".format(round(x[i,:,:].min(), 2))
        variables = [v[n, i, :, :], w[n, i, :, :], T[n, i, :, :], Y[n, i, :, :]]
        domain = [y[i, :, :], z[i, :, :]]
    plot(domain, variables, lims, label_coor, axes)