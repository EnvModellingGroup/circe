In here are your main sim directories. 
For each sim, simply copy the base directory.

stag\_water is a stagnant water test to check your mesh
works well with wetting and drying. If you don't have
wetting and drying enabled, you can delete this.

global\_params.py contains the main parameters to alter, 
as well as allows you to enable processes. Edit this
as you need to.

Inside each simulation folder should be a params.py file.
In this file, you can overwrite any of the variables in
global\_params.py which will be used for that simulation only,
e.g. using a different mesh, bathymetry files, forcing file, etc

Finally, each simulation can be set up with the command line, but
this is best left for overwriting one or two parameters only!

The simulaiton will use parameters in the following order:
 - command line set
 - params.py
 - global\_params.py
