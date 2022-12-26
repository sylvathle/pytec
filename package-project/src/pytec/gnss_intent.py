#############################################################################3
#
# MIT License
#
# Copyright (c) 2021 Sylvain Blunier
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#
#############################################################################3

import os,sys,shutil

#import sys
import numpy as np
import georinex as gr
import datetime
import math
import numpy as np
import urllib.request
#from unlzw import unlzw

import pandas as pd

import random
import time

import matplotlib.pyplot as plt
from . import stations as st

from os import listdir
from os.path import isfile, join

#Ionosphere altitude
Alt_i = 300000
# Earth Radius
R_E = 6370000
omega_e = 7.292115e-5
GM = 3.986005e14
C20 = -1082.63e-6
a_major_ellipsoid_axis = 6378136
WG84_to_PZ90=np.matrix([[1,-0.33/3600.0,0],[0.33/3600.0,1,0],[0,0,1]])


#root_dir = "/media/sblunier/CSN_TEC/"
#root_dir = os.getenv("ROOT_RINEX")
#root_dir = "/home/sblunier/Work/spaceweather/tec/"
#dir_rinex_nav = "/home/sblunier/Work/spaceweather/data/GNSS/"
#dir_rinex_nav = root_dir+"2020/349/"
#dir_sat_bias = "/home/sblunier/Work/spaceweather/data/GNSS/CODE/"
#dir_sat_bias = dir_rinex_nav

def F1(pos,v):
	return np.array(v)

def F2(pos,v,a):
	r=np.linalg.norm(pos)
	xpp = -GM*pos[0]/r**3 + 1.5*C20*GM*a_major_ellipsoid_axis**2 * pos[0] * abs(1-5*pos[2]**2/r**2) / r**5  + omega_e**2*pos[0] + 2*omega_e*v[1]+a[0]
	ypp = -GM*pos[1]/r**3 + 1.5*C20*GM*a_major_ellipsoid_axis**2 * pos[1] * abs(1-5*pos[2]**2/r**2) / r**5 + omega_e**2*pos[1] - 2*omega_e*v[0]+a[1]
	zpp = -GM*pos[2]/r**3 + 1.5*C20*GM*a_major_ellipsoid_axis**2 * pos[2] * abs(3-5*pos[2]**2/r**2) / r**5
	return np.array([xpp,ypp,zpp])

def getDateFromSat(weekNumber,seconds):
	date_ref = datetime.datetime(1980,1,6,0,0,0)
	days = weekNumber*7
	sec = int(seconds) #- 604800
	while sec > 86399:
		days = days + 1
		sec -= 86400
	date = date_ref + datetime.timedelta(days=days) + datetime.timedelta(seconds=sec)
	return date

def gps_nav_to_XYZ(data,date):
    Toe = data['Toe']
    TGD = data['TGD']
    IDOT = data['IDOT']
    IODC = data['IODC']
    GPSWeek = data['GPSWeek']
    TransTime = data['TransTime']
    #clock correction biais
    clockBiais = data['SVclockBias']
    #clock correction drift
    clockDrift = data['SVclockDrift']
	#clock correction drift rate
    clockDriftRate = data['SVclockDriftRate']
	#Square root of the semi-major axis m^{1/2}
    sqrtA =  data['sqrtA']
	#Eccentricity
    Eccentricity =  data['Eccentricity']
	#Inclination angle at reference time
    Io =  data['Io']
	#Longitude of ascending node at reference time
    Omega0 =  data['Omega0']
	#Argument of perigee (semicircles)
    omega = data['omega']
	#Mean anomaly at reference time (semicircles)
    M0 = data['M0']
	#Mean motion difference from computed value
    DeltaN = data['DeltaN']
	#Rate of change of right ascension
    OmegaDot = data['OmegaDot']
	#Rate of change of inclination
    IDOT = data['IDOT']
	#Amplitude of the sine harmonic correction term to the argument of latitude (rad)
    Cus = data['Cus']
	#Amplitude of the cosine harmonic correction term to the argument of latitude (rad)
    Cuc = data['Cuc']
	#Amplitude of the sine harmonic correction term to the angle of inclination (rad)
    Cis = data['Cis']
	#Amplitude of the cosine harmonic correction term to the angle of inclination (rad)
    Cic = data['Cic']
	#Amplitude of the sine harmonic correction term to the orbit radius (m)
    Crs = data['Crs']
	#Amplitude of the cosine harmonic correction term to the orbit radius (m)
    Crc = data['Crc']

    d0 = getDateFromSat(GPSWeek,Toe)

    tk = (date-d0).seconds

    n0 = math.sqrt(GM)/sqrtA**3
    n = n0 + DeltaN

    Mk = M0 + n*tk
    Ek = Mk
    for i in range(3): Ek = Mk + Eccentricity*math.sin(Ek)

    sin_vk = math.sqrt(1-Eccentricity**2)*math.sin(Ek) / (1.0-Eccentricity*math.cos(Ek))
    cos_vk = (math.cos(Ek)-Eccentricity)/(1.0-Eccentricity*math.cos(Ek))
    vk=0
    if sin_vk>0: vk = math.acos(cos_vk)
    else: vk = -math.acos(cos_vk)

    Phik = vk + omega
    delta_uk = Cuc * math.cos(2*Phik) + Cus * math.sin(2*Phik)
    delta_rk = Crc * math.cos(2*Phik) + Crs * math.sin(2*Phik)
    delta_ik = Cic * math.cos(2*Phik) + Cis * math.sin(2*Phik)

    uk = Phik + delta_uk
    rk = sqrtA**2 * (1 - Eccentricity * math.cos(Ek)) + delta_rk
    ik = Io + IDOT*tk + delta_ik

    Xk_p = rk * math.cos(uk)
    Yk_p = rk * math.sin(uk)

    Omegak = Omega0 + (OmegaDot-omega_e)*tk - omega_e * Toe

    Xk = Xk_p * math.cos(Omegak) - Yk_p * math.sin(Omegak)* math.cos(ik)
    Yk = Xk_p * math.sin(Omegak) + Yk_p * math.cos(Omegak)* math.cos(ik)
    Zk = Yk_p * math.sin(ik)
    return [Xk,Yk,Zk]

dir_sats = "sats/"


class gnss:
    '''
        gnss Class to perform eficiently the position and elevation of the
        satellites
        for now it only works on GPS constellation, future version should deal
        with GLONASS
    '''
    def __init__(self,f_nav=""):

        # Create buffer directory if does not exist to store satellite
        # positions
        #if not os.path.exists(dir_sats): os.makedirs(dir_sats)

        self.f_rinex_nav = f_nav
        self.nav = None
        try:
            self.nav = gr.load(self.f_rinex_nav)
        except ValueError:
            print ("Warning: Problem with navigation file",f_nav)

        self.nav_sats = {}
        self.xyz_sats = {}
        for i in range(1,33):
            #Name of the satellite
            sat = ""
            if i<10: sat = "G0"+str(i)
            else: sat = "G"+str(i)

            self.nav_sats[sat]=pd.DataFrame()

            d_sat = {"time":[],"X":[],"Y":[],"Z":[],"elevation":[]}
            self.xyz_sats[sat]=pd.DataFrame(d_sat).set_index("time")

    def getPos(self,sat,date):
        doy = str(date.timetuple().tm_yday)
        year = str(date.year)
        #gps_dir = root_dir + year+ "/"+ doy
        if len(self.nav_sats[sat])==0 or date<self.nav_sats[sat].index[0]:

        
            data = self.nav.sel(method='pad')
            #print (data)

            #sat_date = date.strftime("_%Y%m%d")
            #print (sat)
            d_csv = {"time":[],"Toe":[],"TGD":[],"IDOT":[],"IODC":[],
					"GPSWeek":[],"TransTime":[],"SVclockBias":[],"SVclockDrift":[],
					"SVclockDriftRate":[],"sqrtA":[],"Eccentricity":[],
					"Io":[],"Omega0":[],"omega":[],"M0":[],"DeltaN":[],
					"OmegaDot":[],"IDOT":[],"Cus":[],"Cuc":[],"Cis":[],
					"Crs":[],"Crc":[],"Cic":[],"sat":[]}

            nav_time_coords = self.nav.sel(sv=sat).to_dataframe().dropna().index.tolist()
            #print (nav_time_coords)

            for t in nav_time_coords:
                data = self.nav.sel(sv=sat,time=t,method='pad')
                if (math.isnan(data['Eccentricity'].values)): continue
                for k in d_csv.keys(): 
                    if k!="sat": d_csv[k].append(data[k].values)
                d_csv["sat"].append(sat)
            if len(d_csv["time"])<2:
                print ("Warning: No heatly navigation information for satellite",sat, "in file",self.f_rinex_nav)
            
            #print (pd.DataFrame(d_csv))
            
            #sys.exit()
          	  #print (d_csv)        
        
        
            #csv_nav_dir = st.root_dir + year+"/"+ doy +"/" + "GPS/" + sat + "/"
            
            #print (csv_nav_dir)
            #if not os.path.exists(csv_nav_dir):  
                #print ("make csv 1", date)
            #    makeCSV(date)
            #nav_csv_files = [f for f in listdir(csv_nav_dir) if isfile(join(csv_nav_dir, f))]
            #while (len(nav_csv_files)==0):
            #    doy=str(int(doy)-1)
            #    csv_nav_dir = st.root_dir + year+"/"+ doy +"/" + "GPS/" + sat + "/"
            #    if not os.path.exists(csv_nav_dir):  
            #        print ("make csv 2", date)
            #        makeCSV(date)
            #    nav_csv_files = [f for f in listdir(csv_nav_dir) if isfile(join(csv_nav_dir, f))]
            if len(d_csv["time"])==0: return [float("nan"),float("nan"),float("nan")]
            #csv_nav = csv_nav_dir+nav_csv_files[0]
            self.nav_sats[sat] = pd.DataFrame(d_csv) #pd.read_csv(csv_nav)
            self.nav_sats[sat]["time"]=pd.to_datetime(self.nav_sats[sat]["time"])
            self.nav_sats[sat].set_index("time",inplace=True)
        dt=datetime.timedelta(days=2)
        closest_t = date
        for t in self.nav_sats[sat].index:
            # Test for dates before time seeked
            #if d<t:
            if (t-date)<dt:
                closest_t = t
                dt = t-date
        data = self.nav_sats[sat].loc[t]
        return gps_nav_to_XYZ(data,date)

    # return elevation of satellite given the ID os the satellite, the position of the antenna.
    def getElevation(self,sat,pos_antena,date,mode="RAD"):
        #return 1.0, [20000000,2000000,-1000000]
        pos_sat=self.getPos(sat,date)
        n = np.array([pos_antena[0],pos_antena[1],pos_antena[2]])
        v = np.array([pos_sat[0]-pos_antena[0],pos_sat[1]-pos_antena[1],pos_sat[2]-pos_antena[2]])
        nv = np.dot(n,v)
        sin_phi = nv/(np.linalg.norm(n)*np.linalg.norm(v))
        return math.asin(sin_phi),pos_sat


def resetCSV(date):
    #gps_dir = dir_rinex_nav+"GPS/"

    gps_dir = st.root_dir + str(date.year)+"/"+ str(date.timetuple().tm_yday) +"/GPS/"
    print ("In resetcsv",gps_dir)
    try: os.mkdir(gps_dir)
    except OSError as e:
        if e.errno!=17: print ("FAIL creation of directory "+gps_dir, e )
    else: print ("Successfully created the directory "+gps_dir)
    for n_sat in range(1,33):
        if n_sat<10: sat = "G0" + str(n_sat)
        else: sat = "G" + str(n_sat)
        sat_dir = gps_dir+sat
        try: os.mkdir(sat_dir)
        except OSError as e:
            if e.errno!=17: print ("FAIL creation of directory "+gps_dir, err)
        else: print ("Successfully created the directory "+sat_dir)
        for filename in os.listdir(sat_dir):
            file_path = os.path.join(sat_dir, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path): os.unlink(file_path)
                elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e: print('Failed to delete %s. Reason: %s' % (file_path, e))

def get_pos_antenna(station,year,doy):
    df_antennas = pd.read_csv(st.root_dir+str(year)+"/"+str(doy)+"/antennas.csv")
    df_antennas.set_index("station",inplace=True)
    pos_antena = []
    pos_antena.append(df_antennas.loc[station,"x"])
    pos_antena.append(df_antennas.loc[station,"y"])
    pos_antena.append(df_antennas.loc[station,"z"])
    return pos_antena

def getIonosphereIntersec(pos_antena,pos_sat,h=400000):
    R_I = R_E + h
    #### Calculation of the position of the intersection of Satellite-Antena line with ionosphere which is a second degree equation
    Xas = pos_sat[0]-pos_antena[0]
    # return elevation of satellite given the ID os the satellite, the position of the antenna.
    Yas = pos_sat[1]-pos_antena[1]
    Zas = pos_sat[2]-pos_antena[2]
    AS_squared = Xas**2 + Yas**2 + Zas**2
    #Parameters of the polynomial to solve
    a_poly = AS_squared
    b_poly = 2*Xas*(pos_antena[0]*Xas + pos_antena[1]*Yas + pos_antena[2]*Zas)
    c_poly = Xas**2 * ( pos_antena[0]**2 + pos_antena[1]**2+ pos_antena[2]**2 - R_I**2 )
    Delta = b_poly**2 - 4*a_poly * c_poly
    x_ion, y_ion, z_ion = float('nan'),float('nan'),float('nan')
    pos_ion = [float('nan'),float('nan'),float('nan')]
    if Delta > 0:
        x_ion_1 = (-b_poly+np.sqrt(Delta))/(2*a_poly) + pos_antena[0]
        y_ion_1 = Yas/Xas*(x_ion_1-pos_antena[0])+pos_antena[1]
        z_ion_1 = Zas/Xas*(x_ion_1-pos_antena[0])+pos_antena[2]

        x_ion_2 = (-b_poly-np.sqrt(Delta))/(2*a_poly) + pos_antena[0]
        y_ion_2 = Yas/Xas*(x_ion_2-pos_antena[0])+pos_antena[1]
        z_ion_2 = Zas/Xas*(x_ion_2-pos_antena[0])+pos_antena[2]

        # point must be between satellite and antena, este must be the other solution
        X1xposAntena=(x_ion_1-pos_antena[0])*Xas + (y_ion_1-pos_antena[1])*Yas + (z_ion_1-pos_antena[2])*Zas
        X2xposAntena=(x_ion_2-pos_antena[0])*Xas + (y_ion_2-pos_antena[1])*Yas + (z_ion_2-pos_antena[2])*Zas
        if X1xposAntena<0 or (X1xposAntena>0 and X2xposAntena>0 and X2xposAntena<X1xposAntena): pos_ion = [x_ion_2,y_ion_2,z_ion_2]
        else: pos_ion = [x_ion_1,y_ion_1,z_ion_1]
    elif Delta==0:
        x_ion = -b_poly/(2*a_poly)
        pos_ion = [x_ion,Yas/Xas*(x_ion-pos_antena[0])+pos_antena[1],Zas/Xas*(x_ion-pos_antena[0])+pos_antena[2]]
    return pos_ion

def makeCSV(date):
    resetCSV(date)
    gps_dir = st.root_dir + str(date.year)+"/"+ str(date.timetuple().tm_yday)+"/GPS"

    print ("in makecsv", gps_dir)

    doy = str(date.timetuple().tm_yday)
    while len(doy)<3: doy = "0" + doy
    year = str(date.year)
	#i_c = 0
    print ("in make csv",doy,year)

    print ("in make csv", os.listdir(gps_dir))
	#print (dir_rinex_nav)
    for file_nav in os.listdir(gps_dir):
        print (file_nav)
        if file_nav[-1]!="n":
            continue
        pos_file_nav = gps_dir+"/"+file_nav
        print (pos_file_nav)

        station = file_nav[:4]
        try:
            nav = gr.load(pos_file_nav)
        except ValueError:
            continue
        #print (nav)
        data = nav.sel(method='pad')
        #print (data)

        #sat_date = date.strftime("_%Y%m%d")
        for sat in nav.coords["sv"].values:
            #print (sat)
            #csv_file = gps_dir + "/GPS/" + sat + "/"+station+".csv"
            csv_file = gps_dir + "/" + sat + "/"+station+".csv"
            d_csv = {"time":[],"Toe":[],"TGD":[],"IDOT":[],"IODC":[],
						"GPSWeek":[],"TransTime":[],"SVclockBias":[],"SVclockDrift":[],
						"SVclockDriftRate":[],"sqrtA":[],"Eccentricity":[],
						"Io":[],"Omega0":[],"omega":[],"M0":[],"DeltaN":[],
						"OmegaDot":[],"IDOT":[],"Cus":[],"Cuc":[],"Cis":[],
						"Crs":[],"Crc":[],"Cic":[],"sat":[]}

            nav_time_coords = nav.sel(sv=sat).to_dataframe().dropna().index.tolist()
            #print (nav_time_coords)

            for t in nav_time_coords:
                data = nav.sel(sv=sat,time=t,method='pad')
                if (math.isnan(data['Eccentricity'].values)): continue
                for k in d_csv.keys(): d_csv[k].append(data[k].values)
                d_csv["sat"].append(sat)
            if len(d_csv["time"])<2: continue
            #print (d_csv)
            pd.DataFrame(d_csv).to_csv(csv_file,index=False)


def getGPSSatPosition(sat,d,pos_antena=None):
	#csv_nav_dir = dir_rinex_nav + "GPS/" + sat + "/"
	csv_nav_dir = st.root_dir + str(d.year)+"/"+ str(d.timetuple().tm_yday) +"/" + "GPS/" + sat + "/"
	if not os.path.isdir(csv_nav_dir):
		print (f"Calculating RINEX navigation to csv for year:{d.year}, doy:{d.timetuple().tm_yday}")
		makeCSV(d)
		print ("...Done")
	nav_csv_files = [f for f in listdir(csv_nav_dir) if isfile(join(csv_nav_dir, f))]
	if len(nav_csv_files)==0:
		print (f"Calculating RINEX navigation to csv for year:{d.year}, doy:{d.timetuple().tm_yday}")
		makeCSV(d)
		print ("...Done")
		nav_csv_files = [f for f in listdir(csv_nav_dir) if isfile(join(csv_nav_dir, f))]
	csv_nav = csv_nav_dir+nav_csv_files[0]

	df_data = pd.read_csv(csv_nav)
	df_data["time"]=pd.to_datetime(df_data["time"])
	df_data.set_index("time",inplace=True)
	dt=datetime.timedelta(days=2)
	closest_t = d
	for t in df_data.index:
		# Test for dates before time seeked
		#if d<t:
		if (t-d)<dt:
			closest_t = t
			dt = t-d

	data = df_data.loc[t]
	Toe = data['Toe']
	TGD = data['TGD']
	IDOT = data['IDOT']
	IODC = data['IODC']
	GPSWeek = data['GPSWeek']
	TransTime = data['TransTime']
	#clock correction biais
	clockBiais = data['SVclockBias']
	#clock correction drift
	clockDrift = data['SVclockDrift']
	#clock correction drift rate
	clockDriftRate = data['SVclockDriftRate']
	#Square root of the semi-major axis m^{1/2}
	sqrtA =  data['sqrtA']
	#Eccentricity
	Eccentricity =  data['Eccentricity']
	#Inclination angle at reference time
	Io =  data['Io']
	#Longitude of ascending node at reference time
	Omega0 =  data['Omega0']
	#Argument of perigee (semicircles)
	omega = data['omega']
	#Mean anomaly at reference time (semicircles)
	M0 = data['M0']
	#Mean motion difference from computed value
	DeltaN = data['DeltaN']
	#Rate of change of right ascension
	OmegaDot = data['OmegaDot']
	#Rate of change of inclination
	IDOT = data['IDOT']
	#Amplitude of the sine harmonic correction term to the argument of latitude (rad)
	Cus = data['Cus']
	#Amplitude of the cosine harmonic correction term to the argument of latitude (rad)
	Cuc = data['Cuc']
	#Amplitude of the sine harmonic correction term to the angle of inclination (rad)
	Cis = data['Cis']
	#Amplitude of the cosine harmonic correction term to the angle of inclination (rad)
	Cic = data['Cic']
	#Amplitude of the sine harmonic correction term to the orbit radius (m)
	Crs = data['Crs']
	#Amplitude of the cosine harmonic correction term to the orbit radius (m)
	Crc = data['Crc']

	d0 = getDateFromSat(GPSWeek,Toe)

	tk = (d-d0).seconds

	n0 = math.sqrt(GM)/sqrtA**3
	n = n0 + DeltaN

	Mk = M0 + n*tk
	Ek = Mk
	for i in range(3): Ek = Mk + Eccentricity*math.sin(Ek)

	sin_vk = math.sqrt(1-Eccentricity**2)*math.sin(Ek) / (1.0-Eccentricity*math.cos(Ek))
	cos_vk = (math.cos(Ek)-Eccentricity)/(1.0-Eccentricity*math.cos(Ek))
	vk=0
	if sin_vk>0: vk = math.acos(cos_vk)
	else: vk = -math.acos(cos_vk)

	Phik = vk + omega
	delta_uk = Cuc * math.cos(2*Phik) + Cus * math.sin(2*Phik)
	delta_rk = Crc * math.cos(2*Phik) + Crs * math.sin(2*Phik)
	delta_ik = Cic * math.cos(2*Phik) + Cis * math.sin(2*Phik)

	uk = Phik + delta_uk
	rk = sqrtA**2 * (1 - Eccentricity * math.cos(Ek)) + delta_rk
	ik = Io + IDOT*tk + delta_ik

	Xk_p = rk * math.cos(uk)
	Yk_p = rk * math.sin(uk)

	Omegak = Omega0 + (OmegaDot-omega_e)*tk - omega_e * Toe

	Xk = Xk_p * math.cos(Omegak) - Yk_p * math.sin(Omegak)* math.cos(ik)
	Yk = Xk_p * math.sin(Omegak) + Yk_p * math.cos(Omegak)* math.cos(ik)
	Zk = Yk_p * math.sin(ik)
	return [Xk,Yk,Zk]



def getElevation(sat,pos_antena,date,mode="RAD"):
	#return 1.0, [20000000,2000000,-1000000]
	pos_sat=getGPSSatPosition(sat,date,pos_antena)
	n = np.array([pos_antena[0],pos_antena[1],pos_antena[2]])
	v = np.array([pos_sat[0]-pos_antena[0],pos_sat[1]-pos_antena[1],pos_sat[2]-pos_antena[2]])
	nv = np.dot(n,v)
	sin_phi = nv/(np.linalg.norm(n)*np.linalg.norm(v))
	if mode=="DEG": return math.asin(sin_phi)*180.0/math.pi, pos_sat
	return math.asin(sin_phi),pos_sat

def getBias(sat,date,mode="P1P2"):
	if mode=="P1P2": filename = mode+date.strftime("%y%m")+".DCB"
	elif mode=="P1C1": filename = mode+date.strftime("%y%m")+".DCB"
	path_local_file = st.root_dir + str(date.year)+"/"+ str(date.timetuple().tm_yday) +"/" + filename
	if not os.path.exists(path_local_file):
		compressed_file = filename + ".Z"
		url_compressed_file = "http://ftp.aiub.unibe.ch/CODE/" + str(date.year) + "/" + compressed_file
		path_compressed_file = dir_sat_bias + compressed_file
		#print (path_local_file)
		try:
			print ("Download " + compressed_file)
			urllib.request.urlretrieve(url_compressed_file,path_compressed_file)
		except urllib.error.URLError as e:
			print ("Problem in retrieving: " + url_compressed_file +"\n" + e.reason)
			return 0
		with open(path_compressed_file,'rb') as fh:
			f=open(path_local_file,"w")
			f.write(unlzw(fh.read()).decode("utf-8"))
			f.close()
		os.remove(path_compressed_file)
	f = open(path_local_file)
	for lin in f:
		splt_lin = lin.split()
		if len(splt_lin)==0: continue
		if splt_lin[0]!=sat: continue
		else: return float(splt_lin[1])
	print ("No biais found")
	return 0.0

def getBias_fromfile(sat,f_bias):
    path_local_file = f_bias
    f = open(path_local_file)
    for lin in f:
         splt_lin = lin.split()
         if len(splt_lin)==0: continue
         if splt_lin[0]!=sat: continue
         else: return float(splt_lin[1])
    print ("No biais found")
    return 0.0


'''
	return list of dictionnaries, each correspond to information of the arcs in he time series
	each arc is described with keys:
		- "start" : First point of satellite arc
		- "end" : Last point of satellite arc
		- "tmax": time of maximum elevation
		- "full": True if begining and end of arc is contained is the time serie, if truncated, False
	Input: pandas Serie of elevations with time index.
		t_begin and t_end are to indicate time borders of the time serie.
'''
def get_arcs(elevations,t_begin=None,t_end=None):
    elevations = elevations.dropna()
    values = elevations.values
    index = elevations.index

    if t_begin==None: t_begin = index[0]
    if t_end==None: t_end = index[-1]

    # List of dictionaries containing first and last date of arc and if arc is complete in the time serie
    list_arcs=[]

    # Declare first dict, first arc.
    arc={}

    # Variable to track max elevation during the arc.
    max_elevation = 0
    iin=0

    # If first date has data, the first arc is already running
    # First date of arc is first non nan value, and arc not fully available in time serie
    if len(index)==0 or len(values)==0: return list_arcs
    ioldold,voldold = index[0],values[0]
    iold,vold = index[0],values[0]
    i,v = index[0],values[0]
    #isnan1 = math.isnan(v)
    is_in_arc = False
    if v>0:
        if len(index)>1: inew,vnew=index[1],values[1]
        else:  inew,vnew=index[0],values[0]
        arc["start"]=i
        iin=0
        max_elevation=v
        arc["end"]=None
        if len(values)>1:
            if values[1]<values[0]:
                arc["max_ele"]=values[0]
                arc["tmax"]=index[0]
                arc["imax"]=0
            else:
                arc["max_ele"]=values[1]
                arc["tmax"]=index[1]
                arc["imax"]=0

        if i-t_begin<datetime.timedelta(minutes=5) and v<15.0*math.pi/180: arc["full"]=False
        else: arc["full"]=True
        is_in_arc = True



	# Get elevation sense, rising (positive) or falling (negative)
    last_last_rising=False
    last_rising=False
    rising=False

    if len(elevations)==1:
        arc["tmax"]=index[0]
        arc["imax"]=0
        arc["max_ele"]=values[0]
    for k in range(len(elevations)-1):
		# Get new value of elevation
        inew,vnew = index[k+1],values[k+1]
        rising=vnew>v
		#print (inew,vnew,last_rising,rising,is_in_arc)
		#isnan2 = math.isnan(vnew)
        if not is_in_arc and vnew<0:
            ioldold,voldol=iold,vold
            iold,vold=i,v
            i,v=inew,vnew
            last_rising = rising
            continue

        if not is_in_arc and vnew>0:
            print ("not in arc and vnew>0")
            is_in_arc=True
            arc["start"]=inew
            iin=k
            arc["end"]=None
            arc["tmax"]=i
            arc["imax"]=0
            arc["full"]=True
            max_elevation = v
            arc["max_ele"]=max_elevation

        # if elevation is changing from up to down, max reached
        if k>0 and last_rising and not rising and max_elevation<vnew:
            max_elevation=vnew
            arc["tmax"]=i
            arc["imax"]=k-iin-1
            arc["max_ele"]=max_elevation

		# Case arc is finished when satellite goes down skyline
        if not rising and (vnew<0):
            #print ("going down, skyline passed")
            arc["end"]=i
            list_arcs.append(arc)
            is_in_arc=False
            arc={} # empty dictionnary waiting for next arc

		# Case time serie between falling and next rising is truncated
        if k>1 and  not last_last_rising and rising and (i-iold>datetime.timedelta(hours=1)):
            arc["end"]=iold#inew
            list_arcs.append(arc)
            arc={}
            arc["start"]=i
            iin=k
            arc["end"]=None
            arc["tmax"]=inew
            arc["imax"]=0
            arc["full"]=True
            max_elevation = vnew
            arc["max_ele"]=max_elevation
            is_in_arc=True
			#print ("GAP")
			#print (i,iold,i-iold)
			#print (list_arcs)

        if k>0 and not last_rising and rising and (v<0):
            arc["end"]=i
            list_arcs.append(arc)
            arc={}
            arc["start"]=i
            iin=k
            arc["end"]=None
            arc["tmax"]=i
            arc["imax"]=0
            arc["full"]=True
            max_elevation = v
            arc["max_ele"]=max_elevation
            is_in_arc=True
			#print ("CHANGING ARC v<0")
			#print (i,iold,i-iold)
			#print (list_arcs)

        iold,vold=i,v
        i,v=inew,vnew

        last_last_rising = last_rising
        last_rising = rising


	# Special treatment for last data of time series, if an arc was running end it and set no "full" arc
    if is_in_arc:
        arc["end"]=inew
        if rising:
            arc["tmax"]=inew
            arc["imax"]=k-iin
        else: arc["full"]=True
        arc["max_ele"]=elevations[arc["tmax"]]
        list_arcs.append(arc)

	#for arc in list_arcs: print (arc)
    return list_arcs
