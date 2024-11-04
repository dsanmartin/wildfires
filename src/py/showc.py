import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot(domain, variables, label_coor, axis, t, axes_lims=[[0, 200], [0, 200], [0,20]], limits=False, imshow=False):
    # variable_name = ['u', 'v', 'w', 'T', 'Y']
    # units = ['m/s', 'm/s', 'm/s', 'K', '\%']
    # cmaps = ['viridis', 'viridis', 'viridis', 'jet', 'Oranges']
    # plots = [None, None, None, None, None]
    n_plots = len(variables)
    plots = [None] * n_plots
    x_min_lim, x_max_lim = axes_lims[0]
    y_min_lim, y_max_lim = axes_lims[1]
    z_min_lim, z_max_lim = axes_lims[2]
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

    x, y = domain
    # for i in range(len(variables)):
    i = 0
    for var in variables:
        data = variables[var]['data_plot']
        cmap = variables[var]['cmap']
        lims = variables[var]['lims']
        units = variables[var]['units']
        label = variables[var]['label']
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        if limits:
            if imshow:
                plots[i] = axes[i].imshow(data, cmap=cmap, vmin=lims[0], vmax=lims[1], extent=[x_min_lim, x_max_lim, y_min_lim, y_max_lim], origin='lower', interpolation='none')
            else:
                plots[i] = axes[i].contourf(x, y, data, cmap=cmap, vmin=lims[0], vmax=lims[1])
            m = plt.cm.ScalarMappable(cmap=cmap)
            m.set_clim(lims[0], lims[1])
            plt.colorbar(m, ax=axes[i], cax=cax, label=r'${}$'.format(units))
        else:
            if imshow:
                plots[i] = axes[i].imshow(data, cmap=cmap, extent=[x_min_lim, x_max_lim, y_min_lim, y_max_lim], origin='lower', interpolation='none')
            else:
                plots[i] = axes[i].contourf(x, y, data, cmap=cmap)
            plt.colorbar(plots[i], ax=axes[i], cax=cax, label=r'${}$'.format(units))
        axes[i].set_title(r'${}{}$'.format(label, label_coor))
        if axis == "xy":
            axes[i].set_xlim(x_min_lim, x_max_lim)
            axes[i].set_ylim(y_min_lim, y_max_lim)
        elif axis == "xz":
            axes[i].set_xlim(x_min_lim, x_max_lim)
            axes[i].set_ylim(z_min_lim, z_max_lim)
        elif axis == "yz":
            axes[i].set_xlim(y_min_lim, y_max_lim)
            axes[i].set_ylim(z_min_lim, z_max_lim)
        # if 'x' in axis:
        #     axes[i].set_xlim(x_min_lim, x_max_lim)
        # if 'y' in axis:
        #     axes[i].set_ylim(y_min_lim, y_max_lim)
        # if 'z' in axis:
        #     axes[i].set_ylim(z_min_lim, z_max_lim)
        i += 1
    plt.tight_layout()
    plt.show()
    return None

def load_c(c_dir, Nut=None, topo=False):
    with open(c_dir + 'parameters.txt', 'r') as f:
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
                Nt = int(line.split(",")[3].split(":")[1]) + 1
            if 'Time samples' in line:
                NT = int(line.split(":")[1])
            if "NT" in line:
                NT = int(line.split(":")[1])
                Nt -= 1
    if NT > 1:
        Nt = Nt // NT + 1
    if Nut is not None:
        Nt = Nut + 1
    u_c = np.zeros((Nt, Ny, Nx, Nz))
    v_c = np.zeros((Nt, Ny, Nx, Nz))
    w_c = np.zeros((Nt, Ny, Nx, Nz))
    T_c = np.zeros((Nt, Ny, Nx, Nz))
    Y_c = np.zeros((Nt, Ny, Nx, Nz))
    p_c = np.zeros((Nt, Ny, Nx, Nz))
    # t = np.loadtxt(c_dir + 't.csv', delimiter=',', skiprows=1)
    t = np.fromfile(c_dir + 't.bin', dtype=np.float64)
    # Nz_Y = np.fromfile(c_dir + 'Y_mask.bin', dtype=np.int32).reshape((Nx, Ny)).T    
    for n in range(Nt):
        # data_u = np.loadtxt(c_dir + 'u.csv.{}'.format(n), delimiter=',', skiprows=1)
        # data_v = np.loadtxt(c_dir + 'v.csv.{}'.format(n), delimiter=',', skiprows=1)
        # data_w = np.loadtxt(c_dir + 'w.csv.{}'.format(n), delimiter=',', skiprows=1)
        # data_T = np.loadtxt(c_dir + 'T.csv.{}'.format(n), delimiter=',', skiprows=1)
        # data_Y = np.loadtxt(c_dir + 'Y.csv.{}'.format(n), delimiter=',', skiprows=1)
        # data_p = np.loadtxt(c_dir + 'p.csv.{}'.format(n), delimiter=',', skiprows=1)
        # u_tmp = data_u[:, 3].reshape((Nx, Ny, Nz), order='F')
        # v_tmp = data_v[:, 3].reshape((Nx, Ny, Nz), order='F')
        # w_tmp = data_w[:, 3].reshape((Nx, Ny, Nz), order='F')
        # T_tmp = data_T[:, 3].reshape((Nx, Ny, Nz), order='F')
        # p_tmp = data_p[:, 3].reshape((Nx, Ny, Nz), order='F')
        # Y_tmp = data_Y[:, 3].reshape((Nx, Ny, Nz_Y), order='F')
        data_u = np.fromfile(c_dir + 'u.bin.{}'.format(n), dtype=np.float64)
        data_v = np.fromfile(c_dir + 'v.bin.{}'.format(n), dtype=np.float64)
        data_w = np.fromfile(c_dir + 'w.bin.{}'.format(n), dtype=np.float64)
        data_T = np.fromfile(c_dir + 'T.bin.{}'.format(n), dtype=np.float64)
        data_Y = np.fromfile(c_dir + 'Y.bin.{}'.format(n), dtype=np.float64)
        data_p = np.fromfile(c_dir + 'p.bin.{}'.format(n), dtype=np.float64)
        if n == 0: Nz_Y = data_Y.shape[0] // (Ny * Nx)
        u_tmp = data_u.reshape((Nx, Ny, Nz), order='C')
        v_tmp = data_v.reshape((Nx, Ny, Nz), order='C')
        w_tmp = data_w.reshape((Nx, Ny, Nz), order='C')
        T_tmp = data_T.reshape((Nx, Ny, Nz), order='C')
        p_tmp = data_p.reshape((Nx, Ny, Nz), order='C')
        Y_tmp = data_Y.reshape((Nx, Ny, Nz_Y), order='C')
        u_c[n, :, :, :] = np.transpose(u_tmp, (1, 0, 2))
        v_c[n, :, :, :] = np.transpose(v_tmp, (1, 0, 2))
        w_c[n, :, :, :] = np.transpose(w_tmp, (1, 0, 2))
        T_c[n, :, :, :] = np.transpose(T_tmp, (1, 0, 2))
        p_c[n, :, :, :] = np.transpose(p_tmp, (1, 0, 2))
        Y_c[n, :, :, :Nz_Y] = np.transpose(Y_tmp, (1, 0, 2))
    #     u_c[n, :, :, :] = data_u[:, 3].reshape((Ny, Nx, Nz), order='F')
    #     v_c[n, :, :, :] = data_v[:, 3].reshape((Ny, Nx, Nz), order='F')
    #     w_c[n, :, :, :] = data_w[:, 3].reshape((Ny, Nx, Nz), order='F')
    #     T_c[n, :, :, :] = data_T[:, 3].reshape((Ny, Nx, Nz), order='F')
    #     p_c[n, :, :, :] = data_p[:, 3].reshape((Ny, Nx, Nz), order='F')
    #     Y_c[n, :, :, :Nz_Y] = data_Y[:, 3].reshape((Ny, Nx, Nz_Y), order='F')
    # x = data_p[:, 0].reshape((Ny, Nx, Nz), order='F')
    # y = data_p[:, 1].reshape((Ny, Nx, Nz), order='F')
    # z = data_p[:, 2].reshape((Ny, Nx, Nz), order='F')
    # x_tmp = data_u[:, 0].reshape((Nx, Ny, Nz), order='F')
    # y_tmp = data_u[:, 1].reshape((Nx, Ny, Nz), order='F')
    # z_tmp = data_u[:, 2].reshape((Nx, Ny, Nz), order='F')
    # x = np.transpose(x_tmp, (1, 0, 2))[0,:,0]
    # y = np.transpose(y_tmp, (1, 0, 2))[:,0,0]
    # z = np.transpose(z_tmp, (1, 0, 2))[0,0,:]
    x = np.fromfile(c_dir + 'x.bin', dtype=np.float64)
    y = np.fromfile(c_dir + 'y.bin', dtype=np.float64)
    z = np.fromfile(c_dir + 'z.bin', dtype=np.float64)
    c_vars = [u_c, v_c, w_c, T_c, Y_c, p_c]
    domain = [x, y, z, t]
    out = [c_vars, domain]
    if topo:
        topography = np.fromfile(c_dir + 'topography.bin', dtype=np.float64).reshape(Nx, Ny).T
        out.append(topography)
    return out


"""
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
    Nt = Nt // NT + 0
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
    if n == 0: Nz_Y = data_Y.shape[0] // (Ny * Nx)
    u_tmp = data_u[:, 3].reshape((Nx, Ny, Nz), order='F')
    v_tmp = data_v[:, 3].reshape((Nx, Ny, Nz), order='F')
    w_tmp = data_w[:, 3].reshape((Nx, Ny, Nz), order='F')
    T_tmp = data_T[:, 3].reshape((Nx, Ny, Nz), order='F')
    p_tmp = data_p[:, 3].reshape((Nx, Ny, Nz), order='F')
    Y_tmp = data_Y[:, 3].reshape((Nx, Ny, Nz_Y), order='F')
    u[n, :, :, :] = np.transpose(u_tmp, (1, 0, 2))
    v[n, :, :, :] = np.transpose(v_tmp, (1, 0, 2))
    w[n, :, :, :] = np.transpose(w_tmp, (1, 0, 2))
    T[n, :, :, :] = np.transpose(T_tmp, (1, 0, 2))
    p[n, :, :, :] = np.transpose(p_tmp, (1, 0, 2))
    Y[n, :, :, :Nz_Y] = np.transpose(Y_tmp, (1, 0, 2))

x_tmp = data_u[:, 0].reshape((Nx, Ny, Nz), order='F')
y_tmp = data_u[:, 1].reshape((Nx, Ny, Nz), order='F')
z_tmp = data_u[:, 2].reshape((Nx, Ny, Nz), order='F')
x = np.transpose(x_tmp, (1, 0, 2))[0,:,0]
y = np.transpose(y_tmp, (1, 0, 2))[:,0,0]
z = np.transpose(z_tmp, (1, 0, 2))[0,0,:]
# t = np.linspace(t_min, t_max, Nt)

"""