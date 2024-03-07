import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot(domain, variables, label_coor, axis, t):
    # variable_name = ['u', 'v', 'w', 'T', 'Y']
    # units = ['m/s', 'm/s', 'm/s', 'K', '\%']
    # cmaps = ['viridis', 'viridis', 'viridis', 'jet', 'Oranges']
    # plots = [None, None, None, None, None]
    n_plots = len(variables)
    plots = [None] * n_plots
    x_min_lim, x_max_lim = 0, 200
    if 'z' in axis:
        figsize = (12, 12)
        n_rows, n_cols = n_plots, 1
        sharex, sharey = True, False
    else:
        figsize = (18, 4)
        n_rows, n_cols = 1, n_plots
        sharex, sharey = False, True
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=sharex, sharey=sharey)
    h_var, v_var = axis[0], axis[1]
    # Use t as title
    fig.suptitle(r'Simulation at $t = {}$ s'.format(round(t, 2)), fontsize=16)
    fig.subplots_adjust(top=0.88)
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
    # for i in range(len(variables)):
    i = 0
    for var in variables:
        data = variables[var]['data_plot']
        cmap = variables[var]['cmap']
        lims = variables[var]['lims']
        units = variables[var]['units']
        label = variables[var]['label']
        # plots[i] = axes[i].contourf(x, y, variables[i], cmap=cmaps[i], vmin=lims[i][0], vmax=lims[i][1])
        # m = plt.cm.ScalarMappable(cmap=cmaps[i])
        # m.set_clim(lims[i][0], lims[i][1])
        # divider = make_axes_locatable(axes[i])
        # cax = divider.append_axes("right", size="5%", pad=0.1)
        # plt.colorbar(m, ax=axes[i], cax=cax, label=r'${}$'.format(units[i]))
        # axes[i].set_title(r'${}{}$'.format(variable_name[i], label_coor))
        # axes[i].set_xlim(x_min_lim, x_max_lim)
        plots[i] = axes[i].contourf(x, y, data, cmap=cmap, vmin=lims[0], vmax=lims[1])
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_clim(lims[0], lims[1])
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        plt.colorbar(m, ax=axes[i], cax=cax, label=r'${}$'.format(units))
        axes[i].set_title(r'${}{}$'.format(label, label_coor))
        axes[i].set_xlim(x_min_lim, x_max_lim)
        i += 1
    plt.tight_layout()
    plt.show()


if len(sys.argv) < 3:
    print("Usage: python test.py <axes> <v1,v2,...>")
    sys.exit(1)

# Read the data
data_dir = 'data/output/'
params_path = data_dir + 'parameters.txt'

with open(params_path, 'r') as f:
    # Read lines
    lines = f.readlines()
    # Read the number of cells in each direction
    for line in lines:
        line = line.strip()
        if 'Nx' in line:
            Nx = int(line.split(",")[0].split(":")[1])
        if 'Ny' in line:
            Ny = int(line.split(",")[1].split(":")[1])
        if 'Nz' in line:
            Nz = int(line.split(",")[2].split(":")[1])
        if 'Nt' in line:
            Nt = int(line.split(",")[3].split(":")[1])
        if 'NT' in line:
            NT = int(line.split(":")[1])

# Nx, Ny, Nz, Nt = 128, 128, 64, 21
# Nz_Y = 4
# t_min, t_max = 0, 5
if NT > 1:
    Nt = Nt // NT + 1
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
p = np.zeros((Nt, Ny, Nx, Nz))
t = np.loadtxt(data_dir + 't.csv', delimiter=',', skiprows=1)
for n in range(Nt):
    data_u = np.loadtxt(data_dir + 'u.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_v = np.loadtxt(data_dir + 'v.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_w = np.loadtxt(data_dir + 'w.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_T = np.loadtxt(data_dir + 'T.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_Y = np.loadtxt(data_dir + 'Y.csv.{}'.format(n), delimiter=',', skiprows=1)
    data_p = np.loadtxt(data_dir + 'p.csv.{}'.format(n), delimiter=',', skiprows=1)
    u[n, :, :, :] = data_u[:, 3].reshape((Ny, Nx, Nz), order='F')
    v[n, :, :, :] = data_v[:, 3].reshape((Ny, Nx, Nz), order='F')
    w[n, :, :, :] = data_w[:, 3].reshape((Ny, Nx, Nz), order='F')
    T[n, :, :, :] = data_T[:, 3].reshape((Ny, Nx, Nz), order='F')
    p[n, :, :, :] = data_p[:, 3].reshape((Ny, Nx, Nz), order='F')
    Y[n, :, :, :Nz_Y] = data_Y[:, 3].reshape((Ny, Nx, Nz_Y), order='F')

x = data_u[:, 0].reshape((Ny, Nx, Nz), order='F')
y = data_u[:, 1].reshape((Ny, Nx, Nz), order='F')
z = data_u[:, 2].reshape((Ny, Nx, Nz), order='F')
# t = np.linspace(t_min, t_max, Nt)

axes = sys.argv[1] 
variables_to_plot = sys.argv[2].split(',')
k = 2
j = Ny // 2
i = Nx // 2 - 8

u_lims = [u.min(), u.max()]
v_lims = [v.min(), v.max()]
w_lims = [w.min(), w.max()]
T_lims = [T.min(), T.max()]
Y_lims = [Y.min(), Y.max()]
p_lims = [p.min(), p.max()]
lims = [u_lims, v_lims, w_lims, T_lims, Y_lims]
all_variables = {
    'u': {
        'label': 'u',
        'lims': u_lims,
        'data': u,
        'units': 'm/s',
        'cmap': 'viridis'
    },
    'v': {
        'label': 'v',
        'lims': v_lims,
        'data': v,
        'units': 'm/s',
        'cmap': 'viridis'
    },
    'w': {
        'label': 'w',
        'lims': w_lims,
        'data': w,
        'units': 'm/s',
        'cmap': 'viridis'
    },
    'T': {
        'label': 'T',
        'lims': T_lims,
        'data': T,
        'units': 'K',
        'cmap': 'jet'
    },
    'Y': {
        'label': 'Y',
        'lims': Y_lims,
        'data': Y,
        'units': '\%',
        'cmap': 'Oranges'
    },
    'p': {
        'label': 'p',
        'lims': p_lims,
        'data': p,
        'units': 'Pa',
        'cmap': 'viridis'
    }
}

variables = {k: all_variables[k] for k in variables_to_plot if k in all_variables}

for n in range(0, Nt, 1):
    if axes == 'xy':
        label_coor = "(x, y, {})".format(round(z[:,:,k].min(), 2))
        # variables = [u[n, :, :, k], v[n, :, :, k], w[n, :, :, k], T[n, :, :, k], Y[n, :, :, k]]
        for v in variables_to_plot:
            variables[v]['data_plot'] = variables[v]['data'][n, :, :, k]
        domain = [x[:, :, k], y[:, :, k]]
    elif axes == 'xz':
        label_coor = "(x, {}, z)".format(round(y[:,j,:].min(), 2))
        # variables = [u[n, :, j, :], v[n, :, j, :], w[n, :, j, :], T[n, :, j, :], Y[n, :, j, :]]
        for v in variables_to_plot:
            variables[v]['data_plot'] = variables[v]['data'][n, :, j, :]
        domain = [x[:, j, :], z[:, j, :]]
    elif axes == 'yz':
        label_coor = "({}, y, z)".format(round(x[i,:,:].min(), 2))
        # variables = [u[n, i, :, :], v[n, i, :, :], w[n, i, :, :], T[n, i, :, :], Y[n, i, :, :]]
        for v in variables_to_plot:
            variables[v]['data_plot'] = variables[v]['data'][n, i, :, :]
        domain = [y[i, :, :], z[i, :, :]]
    plot(domain, variables, label_coor, axes, t[n])