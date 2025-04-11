##import required scripts
from thetis import *
import numpy as np
import sys
import slide_movement
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params
# for storm model import atmospheric_forcing only
import atmospheric_forcing
#timestepping options
DT = 180 # reduce if solver does not converge
T_EXPORT = params.output_time
T_END = params.end_time
T_START = params.start_datetime
OUTPUT_DIR = params.output_dir
UTM_ZONE = params.utm_zone
UTM_BAND=params.utm_band
CENT_LAT = params.cent_lat
CENT_LON = params.cent_lon
coord_system = coordsys.UTMCoordinateSystem(utm_zone=UTM_ZONE)
STORM_FILE = params.storm_file
model_type=params.model_type

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

if model_type is tsunami:
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)

# function to set up the Coriolis force
# Depends on a "central" lat/lon point in
# your mesh
def coriolis(mesh2d, lat, lon):
    r = 6371e3
    omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * omega * sin(lat_r)
    beta = (1 / r) * 2 * omega * cos(lat_r)
    x = SpatialCoordinate(mesh2d)
    x_0, y_0, UTM_ZONE, zone_letter = params.from_latlon(lat, lon)
    coriolis_2d = Function(FunctionSpace(mesh2d, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))

    return coriolis_2d

# Set up extra things to export for tsunami only
if model_type is tsunami:
    slide_height = Function(P1_2d, name="slide_height")         # Slide movement
    slide_height_file = File(OUTPUT_DIR + '/slide_height.pvd')

#account for Coriolis code - mesh, centre lat, centre lon
coriolis_2d = coriolis(mesh2d, CENT_LAT, CENT_LON)

# set up storm data #
if model_type is storm or storm_tides:
    CG_2d = FunctionSpace(mesh2d, 'CG', 1)
    CG_2d_vec = VectorFunctionSpace(mesh2d, 'CG',1)
    sim_tz = timezone.pytz.utc
    tau = Function(CG_2d_vec, name="tau")
    pressure = Function(CG_2d, name="pressure")
    tau_file = File(OUTPUT_DIR + '/tau.pvd')
    pressure_file = File(OUTPUT_DIR + '/pressure.pvd')
    era5_file = STORM_FILE
    forcing = atmospheric_forcing.ERA5Interpolator(
	CG_2d,tau,pressure,coord_system, era5_file, T_START
	)
    forcing.set_fields(T_START)
#set up tsunami data
if model_type is tsunami:
    def create_hs(t):
        hs = Function(P1_2d)
        mesh2d = hs.ufl_domain()
        xy_vector = mesh2d.coordinates.dat.data
        hs_vector = hs.dat.data
        assert xy_vector.shape[0] == hs_vector.shape[0]
        for i, xy in enumerate(xy_vector):
            hs_vector[i] = set_slide_height(xy,t,form=params.FORM,profile=params.PROFILE)
        return hs

## end of storm/tsunami data setup ##
# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = T_EXPORT
options.simulation_end_time = T_END
options.output_directory = OUTPUT_DIR
options.check_volume_conservation_2d = True
##options for storm model [need to turn on/off]
if model_type is storm or storm_tides:
    options.atmospheric_pressure = pressure
    options.wind_stress = tau
#options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export = []
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
#the manning function we created in initialisation & loaded above
options.manning_drag_coefficient = manning
#the viscosity 'cushion' we created in initialisation & loaded above
options.horizontal_viscosity = h_viscosity
options.coriolis_frequency = coriolis_2d
options.timestep = DT
if model_type is tsunami:
    options.volume_source_2d = Function(P1_2d)
options.use_automatic_wetting_and_drying_alpha = True
options.wetting_and_drying_alpha_min = Constant(0.1)
options.wetting_and_drying_alpha_max = Constant(75.0)
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
    params.forcing_boundary: {'elev': tidal_elev},
    #set closed boundaries to zero velocity
    1000: {'un': 0.0},
}
##TXPO
data_dir = os.path.join(os.environ.get("DATA", "./data"), "tpxo")
if not os.path.exists(data_dir):
    raise IOError(f"Data directory {data_dir} does not exist")
    tbnd = forcing.TPXOTidalBoundaryForcing(
        tidal_elev,
        T_START,
        coord_system,
        data_dir=data_dir,
        constituents=params.constituents,
        boundary_ids=params.forcing_boundary,
    )
##set initial conditions
if model_type is tsunami:
    solverObj.assign_initial_conditions(uv=as_vector((1e-10, 0.0)), elev = Constant(0.0))
    solverObj.create_equations()
else:
    solverObj.assign_initial_conditions(uv=Constant((1.0e-12,1.0e-12)))


#work out our coords in lat/lon to save doing this every timestep.
if not model_type is tsunami:
    mesh2d = tidal_elev.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    llvector = []
    for i,xy in enumerate(xvector):
        ll = params.utm.to_latlon(xy[0], xy[1], UTM_ZONE, UTM_BAND)
        llvector.append(ll)

def update_forcings(t):
  # Save extra functions to separate files every export time for storms
    if model_type is storm or storm_tides :
        if np.mod(t,T_EXPORT) == 0:
            tau_file.write(tau)
            pressure_file.write(pressure)
#for tsunami models only
    if model_type is tsunami:
        landslide_source = solverObj.timestepper.fields.get('volume_source')
        hs = (create_hs(t + options.timestep) - create_hs(t))/options.timestep
        landslide_source.project(hs)
        uv, elev = solverObj.timestepper.solution.split()
        eta = Function(P1_2d).project(elev)
    # Save extra functions to separate files every export time for tsunami models only
        if np.mod(t,T_EXPORT) == 0:
            slide_height_file.write(slide_height.project(create_hs(t)))

    with timed_stage('update forcings'):
        if model_type is tides:
            print_output("Updating tidal field at t={}".format(t))
        if model_type is storm or storm_tides:
            print_output("Updating tidal & atmos fields at t={}".format(t))
    if not model_type is tsunami:
        tidal_forcing.set_tidal_field(tidal_elev, t, llvector)
    if model_type is storm or storm_tides:
        forcing.set_fields(T_START + t)
    if model_type is tides:
        print_output("Done updating tidal field")
    if model_type is storm or storm_tides:
        print_output("Done updating tidal and atmospheric fields")
solverObj.iterate(update_forcings=update_forcings)
