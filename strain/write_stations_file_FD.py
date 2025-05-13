# Writes stations in SPECFEM STATIONS format for FD calc.
# rewrites the station file for FD strain calculations:
import numpy as np
network = 'YY'
input_stns = [['KS1', 38.5,   140.5],
              ['KS2', 37.75,  139.5],
              ['KS3', 37.0,   140.5],
              ['KS4', 45.0,   135.3],
              ['KS5', 45.0,   134.0],
              ['KS6', 43.7,   134.0],
             ]

dlat = 0.0001
dlon = 0.0001

text_file = open("EXAMPLE_FD_STATIONS", "w")

for istn in input_stns:
    name = istn[0]
    lat  = istn[1]
    lon  = istn[2]

    for x in range(1,6):
        for y in range(1,6):
            text_file.write(f'{name}_{y}{x}  YY  {np.around(lat+(y-3)*dlat,5 )}   {np.around(lon+(x-3)*dlon,5)}    0     0\n')

text_file.close()
