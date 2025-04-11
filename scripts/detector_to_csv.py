#!/usr/bin/env python
from helper_functions import *
import numpy as np
import sys
import os
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params
import argparse
import h5py

def main():

    parser = argparse.ArgumentParser(
         prog="detector_to_csv",
         description="""Convert a detector file to CSV that matches the extract_gauges output"""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )
    parser.add_argument(
        'detector_file',
        help='The model run detector file with the gauges used'
    )

    args = parser.parse_args()
    verbose = args.verbose    
    detector_file = args.detector_file

    start_date = params.start_datetime

    # output dir is where the detector file is
    output_dir, filename = os.path.split(detector_file)

    thetis_gauges = h5py.File(detector_file, 'r+')
    thetis_times = thetis_gauges["time"]
    dt = thetis_times[1] - thetis_times[0]
    spin_up = int(params.spin_up / dt) -1
    model_elev_data = pd.DataFrame({"Time":thetis_times[spin_up:,0]})
    model_speed_data = pd.DataFrame({"Time":thetis_times[spin_up:,0]})


    for name, data in thetis_gauges.items():
        if name == "time":
            continue
        thetis_elev = data[spin_up:,0]
        thetis_u =  data[spin_up:,1]
        thetis_v =  data[spin_up:,2]
        thetis_speed = np.sqrt(thetis_u*thetis_u + thetis_v*thetis_v)
        model_elev_data = pd.concat([model_elev_data, pd.DataFrame({name:thetis_elev})], axis=1)
        model_speed_data = pd.concat([model_speed_data, pd.DataFrame({name:thetis_speed})], axis=1)

    filename = output_dir + "/../model_gauges_elev.csv"
    model_elev_data.to_csv(filename, index=False, header=True)
    filename = output_dir + "/../model_gauges_speed.csv"
    model_speed_data.to_csv(filename, index=False, header=True)


    


if __name__ == "__main__":
    main()
