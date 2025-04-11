import datetime
import utm
from utm import *

# bathymetry parameters
# path relative to the root dir of this template. Leave as mesh/blah.msh in most cases
mesh_file = 'mesh/my_mesh.msh'

# first raster in list is REQUIRED (whole bathymetry), additonal rasters follow
rasters = ["../../data/gbr_400_utm56S.tif", "../../data/gbr_100_utm56S_cropped.tif", "../../data/oti_bathy_utm56S_filled_cropped.tif"]
distances = [500.0, 10.0] #only for additional rasters, if only 1 raster, comment out
forcing_boundary = 666
utm_zone = 56
utm_band="K"
cent_lat = -24.15
cent_lon = 151.8

# model start and end parameters
spin_up = 432000 # 5 days
end_time = 3456000 # 40 days
# year, month, day, hour, min, sec
start_datetime = datetime.datetime(2000,1,1,0,0,0) 
time_diff = 0

# model output parameters
output_dir = "output"
output_time = 900

# physical parameters
constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4']
viscosity = 1.0 # viscosity, obvs. 1.0 is a decent value. 10 is high, 0.001 is very low. Lower means more
manning_drag = 0.025
blend_dist = 50000 # what distance should be used for the boundary blending (m)


