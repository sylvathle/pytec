#import pytec
import os
import numpy as np

from pytec import stations as st
from pytec import tec
import sys


f_base = "/home/sylvain/Documents/spaceweather/tec/testingfolder/"

f_rinex = f_base + "antc3490.20o"
#f_feather = "/home/sylvain/Documents/spaceweather/tec/CSN_TEC/2020/349/AAtest.feather"
f_bias = f_base + "P1P22203.DCB"

tec.rinex_to_feather(f_rinex=f_rinex,f_bias=f_bias)

sys.exit()

def compute_tec(f_rinex,f_feather):
    # Access observation RINEX file
    print (f_rinex)
    tec_station = tec.tec(f_rinex)
    if not tec_station.rinex_to_stec(60): return
    print ("HOLA!!!!")
    #tec_station.sv = ["G22"]
    tec_station.add_satellite_pos()
    print ("HOLAHOLA!!!!")

    tec_station.add_baseline(f_bias="./example/P1P22203.DCB")
    tec_station.add_receiver_bias()
    tec_station.to_feather(f_feather)

i=0
list_stations = []
if len(sys.argv)>=2: list_stations.append(sys.argv[1])
else:
    list_stations = st.get_closest_stations([1491869.7698,-4687403.9181,-4046803.3685])
    print (list_stations)
    #np.random.shuffle(list_stations)


#d=int(sys.argv[2])
#f_rinex = st.root_dir+str(year)+"/"+str(d)+"/" + sys.argv[1] + str(d)+'0.'+str(year)[2:]+'o'
doy_list = [349]

for station in list_stations:
    for doy in doy_list:
        # The source of the rinex data to be converted into TEC data
        f_rinex = st.root_dir+str(year)+"/"+str(doy)+"/" + station + str(doy)+'0.'+str(year)[2:]+'o'
        print ("doy=",doy)
        print (f_rinex)
        # Filename of TEC data in feather format
        f_feather = st.root_dir + str(year) + "/"\
                                + str(doy) + "/" + station + ".feather"
        # Case feather file already exists, go for next doy
        if os.path.exists(f_feather): continue
        if os.path.exists(f_rinex):
            compute_tec(f_rinex,f_feather)




sys.exit()
##if not os.path.exists(f_rinex): continue
compute_tec(f_rinex)
sys.exit()

#doy = int(sys.argv[2])
#list_stations = st.get_closest_stations([1491869.7698,-4687403.9181,-4046803.3685])
print (list_stations)
#list_stations = ["each","gva1","cano","lios","gana"]
for i,station in enumerate(list_stations):#[10:50]:
#if True:
    #station = "savo"
    if os.path.exists(st.root_dir+station+".feather"):
        print (st.root_dir+station+".feather  exists")
        continue
    print (f"Computing station {station}, {i/len(list_stations)*100}%")
