# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 11:54:47 2025

@author: jenny
"""

#only tidal stuff
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


#from tidal model.py
time_step = 180 # reduce if solver does not converge
alpha_min = Constant(0.1)
alpha_max = Constant(75.0)



