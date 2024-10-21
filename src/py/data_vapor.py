import os
import sys
import numpy as np

def write_bov_file(raw_data_file, bov_file, variable_name, t, nx, ny, nz):
    """
    Convert raw C binary data to a .bov file.
    
    Parameters:
    - raw_data_file (str): Path to the input raw data file.
    - bov_file (str): Path to the output .bov file.
    - variable_name (str): Name of the variable (e.g., temperature).
    - nx, ny, nz (int): Dimensions of the data grid.
    - data_type (str): Data type of the raw data (e.g., float, double, int).
    - byte_order (str): Endianness of the data (LITTLE or BIG).
    """
    # Determine the format for the data type
    # data_format = {'float': 'FLOAT', 'double': 'DOUBLE', 'int': 'INT'}[data_type]
    
    # Write the .bov header file
    with open(bov_file, 'w') as bov:
        bov.write(f"TIME: {t}\n")
        bov.write(f"DATA_FILE: {raw_data_file}\n")
        bov.write(f"DATA_SIZE: {nx} {ny} {nz}\n")
        bov.write(f"DATA_FORMAT: DOUBLE\n")
        bov.write(f"VARIABLE: {variable_name}\n")
        # bov.write(f"DATA_ENDIAN: {byte_order}\n")
        # bov.write(f"CENTERING: nodal\n")
        bov.write(f"BRICK_ORIGIN: 0.0 0.0 0.0\n")
        bov.write(f"BRICK_SIZE: {nx}.0 {ny}.0 {nz}.0\n")
    
    print(f"{bov_file} written successfully.")
    

def convert_c_raw_to_bov(input_raw_file, output_bov_file, output_dat_file, variable_name, t, nx, ny, nz, bounds):
    """
    Convert a C raw binary file to .bov format by writing a .dat file and a .bov header.
    
    Parameters:
    - input_raw_file (str): Path to the input raw data file.
    - output_bov_file (str): Path to the output .bov header file.
    - output_dat_file (str): Path to the output .dat binary file.
    - nx, ny, nz (int): Dimensions of the grid.
    - data_type (str): Type of data in the raw file (e.g., float, double).
    """
    # Load the raw data from the C binary file
    # dtype = np.dtype(data_type)
    raw_data = np.fromfile(input_raw_file, dtype=np.float64)
    
    # Fuel doesn't store all the values in Nz, so we need to adjust the size of the array
    if raw_data.size < nx*ny*nz:
        Nz_Y = raw_data.size // (nx*ny)
        tmp = np.zeros((nx, ny, nz))
        tmp[:, :, :Nz_Y] = raw_data.reshape((nx, ny, Nz_Y))
        raw_data = tmp.copy()
    
    # Reshape data into grid dimensions
    # raw_data = raw_data.reshape((nx, ny, nz))   
    
    # output = np.zeros(nx*ny*nz)
    # r = 0
    # for k in range(nz):
    #     for j in range(ny):
    #         for i in range(nx):
    #             output[r] = raw_data[i, j, k]
    #             r += 1
    # output.tofile(output_dat_file)
                
    # Reshape data into grid dimensions and reorder it
    raw_data = raw_data.reshape((nx, ny, nz))
    
    # Filter the data
    i_min, i_max, j_min, j_max, k_min, k_max = bounds
    raw_data = raw_data[i_min:i_max+1, j_min:j_max+1, k_min:k_max+1]
    nnx, nny, nnz = raw_data.shape
    
    # Transpose the data to match the VAPOR format
    raw_data = raw_data.transpose(2, 1, 0).ravel()
    
    # Write the raw data to a .dat file
    raw_data.tofile(output_dat_file)
    print(f"{output_dat_file} written successfully.")
    
    # Get only the filename
    output_dat_file = os.path.basename(output_dat_file)
    
    # Write the corresponding .bov file
    write_bov_file(output_dat_file, output_bov_file, variable_name, t, nnx, nny, nnz)
    
def read_parameters(dir):
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
    
# Read input path
if len(sys.argv) < 2:
    print("Usage: python data_paraview.py input_dir [output_dir]")
    sys.exit(1)
elif len(sys.argv) < 3:
    input_dir = sys.argv[1]
    output_dir = input_dir
else:
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
# Create the output directory if it does not exist
vapor_dir = os.path.join(output_dir, 'vapor')
if not os.path.exists(vapor_dir):
    os.makedirs(vapor_dir)
    
x_min, x_max = 0.0, 200
y_min, y_max = 0.0, 200
z_min, z_max = 0.0, 20
# List all files in the input directory
#files = os.listdir(input_dir)
(x, y, z, t), (Nx, Ny, Nz, Nt) = read_parameters(input_dir)
# Get indexes for the boundaries
i_min = np.where(x >= x_min)[0][0]
i_max = np.where(x <= x_max)[0][-1]
j_min = np.where(y >= y_min)[0][0]
j_max = np.where(y <= y_max)[0][-1]
k_min = np.where(z >= z_min)[0][0]
k_max = np.where(z <= z_max)[0][-1]
bounds = (i_min, i_max, j_min, j_max, k_min, k_max)
variables = ['u', 'v', 'w', 'p', 'T', 'Y']
for n in range(Nt):
    for var in variables:
        input_raw_file = os.path.join(input_dir, f"{var}.bin.{n:1d}")
        output_bov_file = os.path.join(vapor_dir, f"{var}_{n:04d}.bov")
        output_dat_file = os.path.join(vapor_dir, f"{var}_{n:04d}.dat")
        convert_c_raw_to_bov(input_raw_file, output_bov_file, output_dat_file, var, t[n], Nx, Ny, Nz, bounds)
