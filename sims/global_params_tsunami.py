from utm import *

# path relative to the root dir of this sim. Leave as mesh/blah.msh in most cases
mesh_file = 'mesh/my_mesh.msh'
forcing_boundary = 666
utm_zone = 30
utm_band = "U"
cent_lat = 55.696
cent_lon = -1.812
# If nothing happens in that for e.g. 3 hours, then alter this
start_time = 23400 # linked to the forcing data. 
end_time = 54000 # 15 hours
output_dir = "output"
output_time = 60 # 1 minute

#from solid block model
time_step = 2
alpha_min = Constant(0.1)
alpha_max = Constant(75.0)

#from post-processing.py
manning_drag_coefficient = Constant(1.0)
horizontal_viscosity = Constant(1.0)

#from slide_movement.py
slide_breadth = 5000.0
slide_width = 3430.
slide_max_thickness = 500.
#not sure if I should change these??
# start x-location
#slide_start_x = 725000+(b/2.)
# start y-location
#slide_start_y = 0
terminal_velocity = 21.09
