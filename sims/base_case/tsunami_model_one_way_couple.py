from thetis import *
import numpy as np
import boundary_forcing
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import global_params_tsunami


#timestepping options
dt = global_params_tsunami.time_step # reduce if solver does not converge
t_export = global_params_tsunami.output_time
t_end = global_params_tsunami.end_time
t_start = global_params_tsunami.start_time
output_dir = global_params_tsunami.output_dir
utm_zone = global_params_tsunami.utm_zone
utm_band=global_params_tsunami.utm_band
cent_lat = global_params_tsunami.cent_lat
cent_lon = global_params_tsunami.cent_lon

# read bathymetry code
chk = CheckpointFile('bathymetry.h5', 'r')
mesh2d = chk.load_mesh()
bathymetry2d = chk.load_function(mesh2d,'bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = CheckpointFile('viscosity.h5', 'r')
h_viscosity = chk.load_function(mesh2d,'viscosity')
chk.close()
chk = CheckpointFile('manning.h5', 'r')
manning = chk.load_function(mesh2d, 'manning')
chk.close()

# function to set up the Coriolis force
# Depends on a "central" lat/lon point in
# your mesh
def coriolis(mesh, lat, lon):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = global_params_tsunami.from_latlon(lat, lon)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))

    return coriolis_2d

#account for Coriolis code - mesh, centre lat, centre lon
coriolis_2d = coriolis(mesh2d, cent_lat, cent_lon)

# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.manning_drag_coefficient = manning #the manning function we created in initialisation & loaded above
options.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
options.coriolis_frequency = coriolis_2d
options.timestep = dt
options.use_automatic_wetting_and_drying_alpha = True
options.wetting_and_drying_alpha_min = global_params_tsunami.alpha_min
options.wetting_and_drying_alpha_max = global_params_tsunami.alpha_max
options.use_wetting_and_drying = True
options.element_family = "dg-dg"
options.swe_timestepper_type = 'DIRK22'

# boundary conditions
tsunami_elev = Function(FunctionSpace(mesh2d, "CG", 1), name='tsunami_elev')
solverObj.bnd_functions['shallow_water'] = {
    global_params_tsunami.forcing_boundary: {'elev': tsunami_elev},
    1000: {'un': 0.0},
}

# Set up usual SWE terms
solverObj.create_equations()

def update_forcings(t):
    t += t_start
    boundary_forcing.set_tsunami_field(tsunami_elev, t)

update_forcings(0.0)
solverObj.assign_initial_conditions(uv=Constant(("1e-7","0.0")), elev=Constant(0.))

# Run model
solverObj.iterate(update_forcings=update_forcings)


