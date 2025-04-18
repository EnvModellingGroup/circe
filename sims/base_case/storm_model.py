from thetis import *
# this imports our tidal forcing. If you want fes, comment out and uncomment the FES line
import tidal_forcing_tpxo as tidal_forcing
#import tidal_forcing_fes as tidal_forcing
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import global_params_storm
import utm
import atmospheric_forcing

#timestepping options
dt = global_params_storm.time_step # reduce if solver does not converge
t_export = global_params_storm.output_time
t_end = global_params_storm.end_time
output_dir = global_params_storm.output_dir
utm_zone = global_params_storm.utm_zone
utm_band=global_params_storm.utm_band
cent_lat = global_params_storm.cent_lat
cent_lon = global_params_storm.cent_lon

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
    x_0, y_0, utm_zone, zone_letter = global_params_storm.from_latlon(lat, lon)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))

    return coriolis_2d

#account for Coriolis code - mesh, centre lat, centre lon
coriolis_2d = coriolis(mesh2d, cent_lat, cent_lon)

# set up storm data
CG_2d = FunctionSpace(mesh2d, 'CG', 1)
CG_2d_vec = VectorFunctionSpace(mesh2d, 'CG',1)

sim_tz = timezone.pytz.utc
tau = Function(CG_2d_vec, name="tau")
pressure = Function(CG_2d, name="pressure")
tau_file = File(outputdir + '/tau.pvd')
pressure_file = File(outputdir + '/pressure.pvd')

coord_system = coordsys.UTMCoordinateSystem(utm_zone=30, south=True)
start_datetime = datetime.datetime(2020,2,14,0,0,0,tzinfo=sim_tz)
era5_file = global_params_storm.era5_file
forcing = atmospheric_forcing.ERA5Interpolator(CG_2d,tau,pressure,coord_system, era5_file, start_datetime)
forcing.set_fields(t_start)

# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.atmospheric_pressure = pressure
options.wind_stress = tau
#options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export = []
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.manning_drag_coefficient = manning #the manning function we created in initialisation & loaded above
options.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
options.coriolis_frequency = coriolis_2d
options.timestep = dt
options.use_automatic_wetting_and_drying_alpha = True
options.wetting_and_drying_alpha_min = global_params_storm.alpha_min
options.wetting_and_drying_alpha_max = global_params_storm.alpha_max
options.use_wetting_and_drying = True
options.element_family = "dg-dg"
options.swe_timestepper_type = 'DIRK22'


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
    global_params_storm.forcing_boundary: {'elev': tidal_elev},
    #set closed boundaries to zero velocity
    1000: {'un': 0.0},
}


solverObj.assign_initial_conditions(uv=Constant((1.0e-10,0.0)))

#work out our coords in lat/lon to save doing this every timestep.
mesh2d = tidal_elev.function_space().mesh()
xvector = mesh2d.coordinates.dat.data
llvector = []
for i,xy in enumerate(xvector):
    ll = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
    llvector.append(ll)

def update_forcings(t):

  # Save extra functions to separate files every export time
  if np.mod(t,t_export) == 0:
        tau_file.write(tau)
        pressure_file.write(pressure)

  with timed_stage('update forcings'):
    print_output("Updating tidal & atmos fields at t={}".format(t))
    tidal_forcing.set_tidal_field(tidal_elev, t, llvector)
    forcing.set_fields(t_start + t)
    print_output("Done updating tidal and atmospheric fields")

solverObj.iterate(update_forcings=update_forcings)
