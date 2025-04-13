# circe

A modelling template for the [EMRG](https://envmodellinggroup.github.io/)
[thetis](https://thetisproject.org/) models.

[Circe is a Greek goddess, daughter of sun god Helios and Perse, 
one of the three thousand Oceanid nymphs](https://en.wikipedia.org/wiki/Circe).
She lived on an island, Aeaea, and was an enchantress. Around her home prowl 
strangely docile lions and wolves. She lures any who land on 
the island to her home with her lovely singing while weaving 
on an enormous loom, but later drugs them so that they change shape.

However, this model is set up to run tidal, storm or tsunami
simulations. To get started you'll need some data and a mesh.

 - The input folder contains all input files; bathymetry, mesh, etc
 - The data folder contains any other data needed for your
model, including tidal gauges, for example
 - The sims folder contain all your model simulations.
 - The scripts folder contains a bunch of generic processing
and analysis scripts.

To use this, fork this repository on GitHub (or the command line and set your
upstream correctly). You can then set options in the `sim/global_options.py`
file, such as your bathymetry data, mesh file, time step etc, as well as which 
processes you want to enable.

If you want to create a suite of runs, then copy the `base_case`
directory and in the params.py file there, edit any parameters from the `global_params.py`
you wish to change.
Alternatively, use the command line option `--param` to set a name-value pair. 

We have exposed the parameters most users may wish to change so the main `cerci.py`
script should not need to be altered unless you know what you are doing.
