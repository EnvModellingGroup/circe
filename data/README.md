Place any data for post-processing in here, such as tidal gauge data.

For tidal gauge data, the format is a CSV file that at minimum must contain
a Name and X, Y columns, e.g.

Name, X, Y
G1, 10, 10
G2, 10, 20

It can also contain phases and amplitudes of the tides, e.g.

Name, X, Y, M2 amp, M2 phase
G1, 10, 10, 1.2, 345
G2, 10, 20, 2.3, 12

The default tide gauge data file is tide\_gauges.csv
If this file exists, the scripts don't need telling about it

