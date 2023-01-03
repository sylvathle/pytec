#import pytec
import os
import numpy as np

#from pytec import stations as st
from pytec import tec
from pytec import graph
import sys


f_base = "/home/sylvain/Documents/spaceweather/tec/testingfolder/"

f_obs = [f_base + "antc3480.20o", f_base + "antc3490.20o", f_base + "antc3500.20o"]
f_nav = [f_base + "amco3491.20n", f_base + "babj3491.20n", f_base + "amte3481.20n", f_base + "naus3481.20n", f_base + "naus3501.20n"]
#f_feather = "/home/sylvain/Documents/spaceweather/tec/CSN_TEC/2020/349/AAtest.feather"
f_bias = f_base + "P1P22203.DCB"

# Se deberia crear f"antc3490tec.feather"
station_tec = tec.tec(f_obs, f_nav, f_bias)
station_tec.compute_vtec()
print (station_tec.df_obs)

graph.plot_station(station_tec.df_obs,f_base+"antc349")
graph.plot_station(station_tec.df_obs,f_base+"antc349",mozaic=True)


