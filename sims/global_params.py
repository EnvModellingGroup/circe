# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 11:54:47 2025

@author: jenny
"""

# path relative to the root dir of this template. Leave as mesh/blah.msh in most cases
mesh_file = 'mesh/my_mesh.msh'
forcing_boundary = 666
utm_zone = 56
utm_band="K"
cent_lat = -24.15
cent_lon = 151.8
spin_up = 432000 # 5 days
end_time = 3456000 # 40 days
output_dir = "output"
output_time = 900
constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4']
# year, month, day, hour, min, sec
start_datetime = datetime.datetime(2000,1,1,0,0,0) 
time_diff = 0

#get anything with a hard coded option of some kind
#from model.py
#timestepping
dt = 180 # reduce if solver does not converge
# read bathymetry code
chk = CheckpointFile('bathymetry.h5', 'r') #this type of thing can have file specified somewhere?
bathymetry2d = chk.load_function(mesh2d,'bathymetry')
#read viscosity / manning boundaries code
chk = CheckpointFile('viscosity.h5', 'r')
h_viscosity = chk.load_function(mesh2d,'viscosity')
chk = CheckpointFile('manning.h5', 'r')
manning = chk.load_function(mesh2d, 'manning')
##from coriolis func
def coriolis(mesh, lat, lon):
    R = 6371e3
    Omega = 7.292e-5
    
#from solver section
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d'] #probably can leave
options.wetting_and_drying_alpha_min = Constant(0.1)
options.wetting_and_drying_alpha_max = Constant(75.0)
options.element_family = "dg-dg"
options.swe_timestepper_type = 'DIRK22'

alpha_min = Constant(0.1)
alpha_max = Constant(75.0)

options.swe_timestepper_options.solver_parameters = {
      'snes_type': 'newtonls',
      'snes_rtol': 1e-3,
      'ksp_rtol': 1e-3,
      'ksp_type': 'gmres',
      'pc_type': 'fieldsplit',
  }

# set boundary/initial conditions code
tidal_elev = Function(bathymetry2d.function_space())
solverObj.bnd_functions['shallow_water'] = {
    #set open boundaries to tidal_elev function
    params.forcing_boundary: {'elev': tidal_elev},
    #set closed boundaries to zero velocity
    1000: {'un': 0.0},
}


solverObj.assign_initial_conditions(uv=Constant((1.0e-12,1.0e-12)))

