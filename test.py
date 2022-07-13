import os
import numpy as np
import pandas as pd
from pytec import stations as st
from pytec import tec
import sys

st.root_dir = "./example/"

# Route of RINEX file to be computed
f_rinex = "./example/ade30700.22o"
#Route of bias file of satellites, should not be necessary in some future version
f_bias = "./example/P1P22203.DCB"

#Output where will be stored the STEC/VTEC computed
f_feather = "./example/jicamarca_349.feather"

if os.path.exists(f_rinex) and not os.path.exists(f_feather):
    # Tells pytec where is the rinex file to process
    tec_station = tec.tec(f_rinex)
    # Compute STEC of pseudo range and code phase
    if not tec_station.rinex_to_stec(60): sys.exit()
    # Compute position of satellite
    tec_station.add_satellite_pos()
    # Compute baseline from STEC data
    tec_station.add_baseline(f_bias)
    # Add receiver biais to STEC
    tec_station.add_receiver_bias()
    # Convert data to feather
    tec_station.to_feather(f_feather)

print (pd.read_feather(f_feather))
