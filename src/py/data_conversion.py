import sys
import numpy as np

def load_c(c_dir, Nut=None):
    # with open(c_dir + 'parameters.txt', 'r') as f:
    #     # Read lines
    #     lines = f.readlines()
    #     # Read the number of cells in each direction
    #     for line in lines:
    #         line = line.strip()
    #         if 'Nx' in line:
    #             Nx = int(line.split(",")[0].split(":")[1])
    #         if 'Ny' in line:
    #             Ny = int(line.split(",")[1].split(":")[1])
    #         if 'Nz' in line:
    #             Nz = int(line.split(",")[2].split(":")[1])
    #         if 'Nt' in line:
    #             Nt = int(line.split(",")[3].split(":")[1]) + 1
    #         if 'Time samples' in line:
    #             NT = int(line.split(":")[1])
    #         if "NT" in line:
    #             NT = int(line.split(":")[1])
    #             Nt -= 1
    # if NT > 1:
    #     Nt = Nt // NT + 1
    # if Nut is not None:
    #     Nt = Nut + 1
    # Load domain
    x = np.fromfile(c_dir + 'x.bin', dtype=np.float64)
    y = np.fromfile(c_dir + 'y.bin', dtype=np.float64)
    z = np.fromfile(c_dir + 'z.bin', dtype=np.float64)
    t = np.fromfile(c_dir + 't.bin', dtype=np.float64)
    # Domain shapes
    Nx, Ny, Nz, Nt = x.shape[0], y.shape[0], z.shape[0], t.shape[0]
    # Load data
    u = np.zeros((Nt, Ny, Nx, Nz))
    v = np.zeros((Nt, Ny, Nx, Nz))
    w = np.zeros((Nt, Ny, Nx, Nz))
    T = np.zeros((Nt, Ny, Nx, Nz))
    Y = np.zeros((Nt, Ny, Nx, Nz))
    p = np.zeros((Nt, Ny, Nx, Nz))
    # Read the data
    for n in range(Nt):
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
        u[n, :, :, :] = np.transpose(u_tmp, (1, 0, 2))
        v[n, :, :, :] = np.transpose(v_tmp, (1, 0, 2))
        w[n, :, :, :] = np.transpose(w_tmp, (1, 0, 2))
        T[n, :, :, :] = np.transpose(T_tmp, (1, 0, 2))
        p[n, :, :, :] = np.transpose(p_tmp, (1, 0, 2))
        Y[n, :, :, :Nz_Y] = np.transpose(Y_tmp, (1, 0, 2))

    variables = [u, v, w, T, Y, p]
    domain = [x, y, z, t]
    return variables, domain


# Read input path
input_dir = sys.argv[1]
output_dir = sys.argv[2]
# Load data
variables, domain = load_c(input_dir)
x, y, z, t = domain
u, v, w, T, Y, p = variables
# Pack the data
data = {
    'u': u, 'v': v, 'w': w, 'T': T, 'Y': Y, 'p': p
}
domain_data = {
    'x': x, 'y': y, 'z': z, 't': t
}
# Save the data
np.savez(output_dir + 'data.npz', **data)
np.savez(output_dir + 'domain.npz', **domain_data)
