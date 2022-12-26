#import pytec
import os
import numpy as np

from pytec import stations as st
from pytec import tec
import sys


f_base = "/home/sylvain/Documents/spaceweather/tec/testingfolder/"

f_obs = f_base + "antc3490.20o"
f_nav = f_base + "amco3491.20n"
#f_feather = "/home/sylvain/Documents/spaceweather/tec/CSN_TEC/2020/349/AAtest.feather"
f_bias = f_base + "P1P22203.DCB"

# Se deberia crear f"antc3490tec.feather"
tec.rinex_to_feather(f_obs=f_obs,f_nav=f_nav,f_bias=f_bias)


