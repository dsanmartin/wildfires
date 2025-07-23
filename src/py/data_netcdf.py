import os
import sys
import numpy as np
from netCDF4 import Dataset

UNITS = {
    'u': 'm s-1',
    'v': 'm s-1',
    'w': 'm s-1',
    'p': 'Pa',
    'T': 'K',
    'Y': '1'
}
LONG_NAMES = {
    'u': 'X velocity',
    'v': 'Y velocity',
    'w': 'Z velocity',
    'p': 'Pressure',
    'T': 'Temperature',
    'Y': 'Mass fraction of species'
}

def read_parameters(dir, bounds=None):
    with open(dir + 'parameters.txt', 'r') as f:
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
    if NT > 1:
        Nt = Nt // NT + 1
    x = np.fromfile(dir + 'x.bin', dtype=np.float64)
    y = np.fromfile(dir + 'y.bin', dtype=np.float64)
    z = np.fromfile(dir + 'z.bin', dtype=np.float64)
    t = np.fromfile(dir + 't.bin', dtype=np.float64)
    return (x, y, z, t), (Nx, Ny, Nz, Nt)

def get_bounds(x, y, z, bounds):
    i_min = np.where(x >= bounds[0])[0][0]
    i_max = np.where(x <= bounds[1])[0][-1]
    j_min = np.where(y >= bounds[2])[0][0]
    j_max = np.where(y <= bounds[3])[0][-1]
    k_min = np.where(z >= bounds[4])[0][0]
    k_max = np.where(z <= bounds[5])[0][-1]
    return (i_min, i_max+1, j_min, j_max+1, k_min, k_max+1)


def save_variables_to_netcdf(data, domain, output_path):
    x, y, z, t = domain
    nnx, nny, nnz, nnt = x.shape[0], y.shape[0], z.shape[0], t.shape[0]
    # Create a new NetCDF file with the variables in data
    with Dataset(output_path, mode='w', format='NETCDF4') as ncfile:
        # Domain dimensions
        ncfile.createDimension('time', nnt)
        ncfile.createDimension('z', nnz)
        ncfile.createDimension('y', nny)
        ncfile.createDimension('x', nnx)

        # Coordinates variables
        times = ncfile.createVariable('time', 'f8', ('time',))
        levels = ncfile.createVariable('z', 'f8', ('z',))
        ys = ncfile.createVariable('y', 'f8', ('y',))
        xs = ncfile.createVariable('x', 'f8', ('x',))

        # Add data to coordinate variables
        times[:] = t
        levels[:] = z
        ys[:] = y
        xs[:] = x

        # Metadata for CF-compliance
        # Time
        times.units = 's'
        times.standard_name = 'time'
        times.long_name = 'Time'
        times.axis = 'T'
        times.calendar = 'gregorian'
        # Vertical levels
        levels.units = 'm'
        levels.standard_name = 'height'
        levels.long_name = 'Vertical height'
        levels.axis = 'Z'
        levels.positive = 'up'
        # Y and X coordinates
        ys.units = 'm'
        ys.standard_name = 'projection_y_coordinate'
        ys.long_name = 'Y Cartesian coordinate'
        ys.axis = 'Y'
        xs.units = 'm'
        xs.standard_name = 'projection_x_coordinate'
        xs.long_name = 'X Cartesian coordinate'
        xs.axis = 'X'

        # NetCDF variables for each variable in data
        for var_name in data:
            variable = data[var_name]
            # Create variable with zlib compression
            var = ncfile.createVariable(var_name, 'f8', ('time', 'z', 'y', 'x'), zlib=True)
            var.units = variable['unit']
            var.long_name = variable['long_name']
            var.standard_name = var_name
            var.coordinates = 'time z y x'
            var[:, :, :, :] = variable['data']

if __name__ == "__main__":
    # Read input path
    if len(sys.argv) < 2:
        print("Usage: python data_netcdf.py input_dir [output_dir]")
        sys.exit(1)
    elif len(sys.argv) < 3:
        input_dir = sys.argv[1]
        output_dir = input_dir
    else:
        input_dir = sys.argv[1]
        output_dir = sys.argv[2]

    netcdf_dir = os.path.join(output_dir, 'netCDF')
    if not os.path.exists(netcdf_dir):
        os.makedirs(netcdf_dir)
        
    # Bound domain
    x_min, x_max = 0.0, 200
    y_min, y_max = 0.0, 200
    z_min, z_max = 0.0, 40
    bounds = (x_min, x_max, y_min, y_max, z_min, z_max)

    # Read domain from parameters file
    (x, y, z, t), (Nx, Ny, Nz, Nt) = read_parameters(input_dir, bounds)
    
    # Get indexes for the domain bounds
    bounds = get_bounds(x, y, z, bounds)

    # Create bounds for x, y, z
    x_bounds = x[bounds[0]:bounds[1]]
    y_bounds = y[bounds[2]:bounds[3]]
    z_bounds = z[bounds[4]:bounds[5]]
    Nx_, Ny_, Nz_ = len(x_bounds), len(y_bounds), len(z_bounds)
    
    variables = ['u', 'v', 'w', 'p', 'T', 'Y']

    # Read each variable for each time step to numpy array
    var_dict = {}
    for var in variables:
        var_dict[var] = {
            'data': np.zeros((Nt, Nz_, Ny_, Nx_)),
            'unit': UNITS[var],
            'long_name': LONG_NAMES[var]
        }
        
    for n in range(Nt):
        for var in variables:
            input_raw_file = os.path.join(input_dir, f"{var}.bin.{n:1d}")
            if os.path.exists(input_raw_file):
                data = np.fromfile(input_raw_file, dtype=np.float64)
                if data.size < Nx * Ny * Nz:
                    Nz_Y = data.size // (Nx * Ny)
                    tmp = np.zeros((Nx, Ny, Nz))
                    tmp[:, :, :Nz_Y] = data.reshape((Nx, Ny, Nz_Y))
                    data = tmp
                data = data.reshape((Nx, Ny, Nz))
                data = data[bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]]
                data = data.transpose(2, 1, 0)
                var_dict[var]['data'][n, :, :, :] = data
    
    
    output_data_file = os.path.join(netcdf_dir, "data.nc")
    save_variables_to_netcdf(var_dict, (x_bounds, y_bounds, z_bounds, t), output_data_file)
    print(f"NetCDF file with variables written in: {output_data_file}")