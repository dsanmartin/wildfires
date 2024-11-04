import os
import sys
import numpy as np
import vtk
from vtk.util import numpy_support

def load_c(c_dir):
    # Load domain
    x = np.fromfile(c_dir + 'x.bin', dtype=np.float64)
    y = np.fromfile(c_dir + 'y.bin', dtype=np.float64)
    z = np.fromfile(c_dir + 'z.bin', dtype=np.float64)
    t = np.fromfile(c_dir + 't.bin', dtype=np.float64)
    # Domain shapes
    Nx, Ny, Nz, Nt = x.shape[0], y.shape[0], z.shape[0], t.shape[0]
    # Load data
    u = np.zeros((Nt, Nx, Ny, Nz))
    v = np.zeros((Nt, Nx, Ny, Nz))
    w = np.zeros((Nt, Nx, Ny, Nz))
    T = np.zeros((Nt, Nx, Ny, Nz))
    Y = np.zeros((Nt, Nx, Ny, Nz))
    p = np.zeros((Nt, Nx, Ny, Nz))
    # Read the data
    for n in range(Nt):
        data_u = np.fromfile(c_dir + 'u.bin.{}'.format(n), dtype=np.float64)
        data_v = np.fromfile(c_dir + 'v.bin.{}'.format(n), dtype=np.float64)
        data_w = np.fromfile(c_dir + 'w.bin.{}'.format(n), dtype=np.float64)
        data_T = np.fromfile(c_dir + 'T.bin.{}'.format(n), dtype=np.float64)
        data_Y = np.fromfile(c_dir + 'Y.bin.{}'.format(n), dtype=np.float64)
        data_p = np.fromfile(c_dir + 'p.bin.{}'.format(n), dtype=np.float64)
        if n == 0: Nz_Y = data_Y.shape[0] // (Ny * Nx)
        u[n] = data_u.reshape((Nx, Ny, Nz), order='C')
        v[n] = data_v.reshape((Nx, Ny, Nz), order='C')
        w[n] = data_w.reshape((Nx, Ny, Nz), order='C')
        T[n] = data_T.reshape((Nx, Ny, Nz), order='C')
        p[n] = data_p.reshape((Nx, Ny, Nz), order='C')
        Y[n, :, :, :Nz_Y] = data_Y.reshape((Nx, Ny, Nz_Y), order='C')
    variables = [u, v, w, T, Y, p]
    domain = [x, y, z, t]
    return variables, domain

def write_vtk_v1(filename, variables, domain):
    u, v, w, T, Y, p = variables
    x, y, z, t = domain
    
    # Create a VTK temporal structured grid
    temporal_grid = vtk.vtkMultiBlockDataSet()
    # temporal_grid = vtk.vtkTemporalStructuredGrid()

    # Create a VTK structured grid
    structured_grid = vtk.vtkStructuredGrid()    

    # Set the dimensions of the grid
    Nx, Ny, Nz, Nt = x.shape[0], y.shape[0], z.shape[0], t.shape[0]
    structured_grid.SetDimensions(Nx, Ny, Nz)

    # Create points
    points = vtk.vtkPoints()
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                points.InsertNextPoint(x[i], y[j], z[k])
    structured_grid.SetPoints(points)
    
    # Function to add a data array to the grid
    def add_array_to_grid(array, name):
        flat_array = array.flatten(order='C')
        vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
        vtk_array.SetName(name)
        structured_grid.GetPointData().AddArray(vtk_array)

    # Add data arrays to the grid for each timestep
    for n in range(Nt):
        add_array_to_grid(u[n], 'u')
        add_array_to_grid(v[n], 'v')
        add_array_to_grid(w[n], 'w')
        add_array_to_grid(T[n], 'T')
        add_array_to_grid(Y[n], 'Y')
        add_array_to_grid(p[n], 'p')
        
        # Set the time for this timestep
        # temporal_grid.SetTimeValue(t[n], structured_grid)
        temporal_grid.SetBlock(n, structured_grid)

    # Write the grid to a .vtk file
    # writer = vtk.vtkStructuredGridWriter()
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(temporal_grid)
    writer.Write()
    
def write_vtk(filename, variables, domain):
    u, v, w, T, Y, p = variables
    x, y, z, t = domain
    
    # Dimensions of the grid
    Nx, Ny, Nz, Nt = x.shape[0], y.shape[0], z.shape[0], t.shape[0]

    # Create a VTK multi-block dataset to hold multiple structured grids
    multi_block = vtk.vtkMultiBlockDataSet()

    # Iterate over time steps
    for time_index in range(len(t)):
        # Create a VTK structured grid for the current time step
        structured_grid = vtk.vtkStructuredGrid()

        # Set the dimensions of the grid
        structured_grid.SetDimensions(Nx, Ny, Nz)

        # Create points
        points = vtk.vtkPoints()
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    points.InsertNextPoint(x[i], y[j], z[k])
        structured_grid.SetPoints(points)

        # Function to add a data array to the grid
        def add_array_to_grid(array, name):
            flat_array = array.flatten(order='C')
            vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
            vtk_array.SetName(name)
            structured_grid.GetPointData().AddArray(vtk_array)

        # Add data arrays to the grid for the current time step
        add_array_to_grid(u[time_index], 'u')
        add_array_to_grid(v[time_index], 'v')
        add_array_to_grid(w[time_index], 'w')
        add_array_to_grid(T[time_index], 'T')
        add_array_to_grid(Y[time_index], 'Y')
        add_array_to_grid(p[time_index], 'p')

        # Add the structured grid to the multi-block dataset
        multi_block.SetBlock(time_index, structured_grid)

    # Write the multi-block dataset to a .vtm file
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(multi_block)
    writer.Write()
    
def create_vtk_spatio_temporal_data_v1(output_file, variables, domain):
    """
    Create a VTK file for spatio-temporal data using numpy arrays and the vtk library.

    Parameters:
    - u, v, w, T, Y, p: 4D numpy arrays of shape (Nt, Nx, Ny, Nz), representing time-dependent 3D fields.
    - output_file: The path to the output VTK file.
    """
    u, v, w, T, Y, p = variables
    x, y, z, t = domain
    
    Nt, Nx, Ny, Nz = t.shape[0], x.shape[0], y.shape[0], z.shape[0]

    # Create a VTK MultiBlock dataset to hold each timestep
    multi_block = vtk.vtkMultiBlockDataSet()
    
    # Loop over timesteps
    for n in range(Nt):
        # Create a VTK points object to define the grid points
        points = vtk.vtkPoints()
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    points.InsertNextPoint(x[i], y[j], z[k])

        # Create a vtkStructuredGrid to store the grid and data
        structured_grid = vtk.vtkStructuredGrid()
        structured_grid.SetDimensions(Nx, Ny, Nz)
        structured_grid.SetPoints(points)

        # Function to add a scalar field to the structured grid
        def add_scalar_field(name, data):
            data_array = vtk.vtkFloatArray()
            data_array.SetName(name)
            data_array.SetNumberOfComponents(1)
            data_array.SetNumberOfTuples(Nx * Ny * Nz)

            for i in range(Nx):
                for j in range(Ny):
                    for k in range(Nz):
                        idx = i + j * Nx + k * Nx * Ny
                        data_array.SetValue(idx, data[n, i, j, k])

            structured_grid.GetPointData().AddArray(data_array)

        # Add each data field for this timestep to the structured grid
        add_scalar_field("u", u)
        add_scalar_field("v", v)
        add_scalar_field("w", w)
        add_scalar_field("T", T)
        add_scalar_field("Y", Y)
        add_scalar_field("p", p)

        # Add this timestep's structured grid to the multi-block dataset
        multi_block.SetBlock(n, structured_grid)

    # Write the multi-block dataset to a VTK file using vtkXMLMultiBlockDataWriter
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(multi_block)
    writer.Write()
    
def create_vtk_spatio_temporal_data(output_directory, variables, domain):
    """
    Create a VTK file for spatio-temporal data using numpy arrays and the vtk library.

    Parameters:
    - u, v, w, T, Y, p: 4D numpy arrays of shape (Nt, Nx, Ny, Nz), representing time-dependent 3D fields.
    - output_directory: The directory where the VTK files will be saved.
    """
    u, v, w, T, Y, p = variables
    x, y, z, t = domain
    
    Nt, Nx, Ny, Nz = t.shape[0], x.shape[0], y.shape[0], z.shape[0]

    # Create a VTK MultiBlock dataset to hold each timestep
    multi_block = vtk.vtkMultiBlockDataSet()
    
    # Loop over timesteps
    for n in range(Nt):
        # Create a VTK points object to define the grid points
        points = vtk.vtkPoints()
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    points.InsertNextPoint(i, j, k)#x[i], y[j], z[k])

        # Create a vtkStructuredGrid to store the grid and data
        structured_grid = vtk.vtkStructuredGrid()
        structured_grid.SetDimensions(Nx, Ny, Nz)
        structured_grid.SetPoints(points)

        # Function to add a scalar field to the structured grid
        def add_scalar_field(name, data):
            data_array = vtk.vtkFloatArray()
            data_array.SetName(name)
            data_array.SetNumberOfComponents(1)
            data_array.SetNumberOfTuples(Nx * Ny * Nz)

            for i in range(Nx):
                for j in range(Ny):
                    for k in range(Nz):
                        idx = i + j * Nx + k * Nx * Ny
                        data_array.SetValue(idx, data[n, i, j, k])

            structured_grid.GetPointData().AddArray(data_array)

        # Add each data field for this timestep to the structured grid
        add_scalar_field("u", u)
        add_scalar_field("v", v)
        add_scalar_field("w", w)
        add_scalar_field("T", T)
        add_scalar_field("Y", Y)
        add_scalar_field("p", p)

        # Write each structured grid to a separate VTK file
        grid_file = f"{output_directory}/structured_grid_timestep_{n}.vts"
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetFileName(grid_file)
        writer.SetInputData(structured_grid)
        writer.Write()

        # Add this timestep's structured grid to the multi-block dataset
        reader = vtk.vtkXMLStructuredGridReader()
        reader.SetFileName(grid_file)
        reader.Update()
        multi_block.SetBlock(n, reader.GetOutput())

    # Write the multi-block dataset to a VTK file using vtkXMLMultiBlockDataWriter
    multi_block_file = f"{output_directory}/spatio_temporal_data.vtm"
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(multi_block_file)
    writer.SetInputData(multi_block)
    writer.Write()
    
def write_array(x, y, z, data, variable_name, filename):
    # Create a VTK structured grid
    structured_grid = vtk.vtkStructuredGrid()
    # Set the dimensions of the grid
    Nx, Ny, Nz = x.shape[0], y.shape[0], z.shape[0]
    structured_grid.SetDimensions(Nx, Ny, Nz)
    # Create points
    points = vtk.vtkPoints()
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                points.InsertNextPoint(x[i], y[j], z[k])
    structured_grid.SetPoints(points)
    # Add array to the grid
    flat_array = data.flatten(order='C')
    vtk_array = numpy_support.numpy_to_vtk(flat_array, deep=True)
    vtk_array.SetName(variable_name)
    structured_grid.GetPointData().AddArray(vtk_array)
    # Write the grid to a .vtk file
    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(structured_grid)
    writer.Write()
    
def write_each_time(domain, variables, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    x, y, z, t = domain
    u, v, w, T, Y, p = variables
    for n in range(t.shape[0]):
        write_array(x, y, z, u[n], 'u', f"{output_directory}/u_{n}.vtk")
        write_array(x, y, z, v[n], 'v', f"{output_directory}/v_{n}.vtk")
        write_array(x, y, z, w[n], 'w', f"{output_directory}/w_{n}.vtk")
        write_array(x, y, z, T[n], 'T', f"{output_directory}/T_{n}.vtk")
        write_array(x, y, z, Y[n], 'Y', f"{output_directory}/Y_{n}.vtk")
        write_array(x, y, z, p[n], 'p', f"{output_directory}/p_{n}.vtk")
    

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

output_path = output_dir + 'vtk/'
    
variables, domain = load_c(input_dir)  # Replace with actual function to read data
write_each_time(domain, variables, output_path)
