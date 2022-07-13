#from __future__ import print_function
import os,sys
from spawe.dst import download_last_dst,get_dst
from spawe.kp import download_kp,get_kp
from spawe.omni import download_last_omni, get_omni_data
from spawe.swpc import download_swpc_data,get_swpc_data
from spawe.sunspots import download_sunspots,get_sunspots
from spawe.sdo import download_images,make_video

import datetime
import matplotlib.pyplot as plt

import pandas as pd

if not os.environ.get("SW_DATA_PATH"):
    print ("First of all you need to define the root directory where all data will be stored, can be an external disk (consider at least 10GB)")
    answer_for_main_dir = input("Is \""+os.environ.get("HOME")+"/spawe-data\" good for you? (will create it) (y/n)")
    spawe_dir = ""
    if answer_for_main_dir=="y":
        spawe_dir = os.environ.get("HOME")+"/spawe-data/"
    else:
        spawe_dir = input("Where do you want to store your data? (write full route, finish with /. If folder is missing it will be created)")

    answer_for_bashrc = input("Will define permanently a global variable in .bashrc (export SW_DATA_PATH="+spawe_dir+") are you OK with this? (y/n)")

    if (answer_for_bashrc=="y"):
        bashrc = open(os.environ.get("HOME")+'/.bashrc', 'a')
        bashrc.write("\nexport SW_DATA_PATH="+spawe_dir)
        bashrc.close()
    else:
        print ("You will need to run this command anytime you need spawe, or define it manually in you .bashrc")

    print ("Running export SW_DATA_PATH="+spawe_dir)
    os.environ["SW_DATA_PATH"] = spawe_dir
else:
    spawe_dir = os.environ["SW_DATA_PATH"]

if not os.path.exists(spawe_dir):
    try: os.mkdir(spawe_dir)
    except ValueError:
        print ("Problem with "+spawe_dir+" creation, try again with valid path")


### DST testing ####
if False:
    print ("testing dst module")
    start_date = datetime.datetime(2000,12,31,0,0,0,0)
    end_date = datetime.datetime(2001,1,1,0,0,0,0)
    download_last_dst()
    df_dst = get_dst(start_date,end_date)
    print (df_dst)
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #    print (df_dst)
    sys.exit()
    plt.plot(df_dst)
    #plt.show()


### Kp testing ###
if False:
    start_date = datetime.datetime(2000,12,31,0,0,0,0)
    end_date = datetime.datetime.now()
    download_kp()
    df_kp = get_kp(start_date,end_date,numberized=True,hourly=False)
    print (df_kp)
    df_kp = df_kp.resample('1W').mean()
    plt.plot(df_kp)
    plt.show()

### swpc testing ####
if False:
    print ("testing swpc module")
    start_date = datetime.datetime(2021,11,14,0,0,0,0)
    end_date = datetime.datetime.now()
    start_date = end_date - datetime.timedelta(days=1)
    download_swpc_data()
    df_swpc = get_swpc_data(start_date,end_date,["Bz"])
    print (df_swpc)
    #plt.plot(df)
    #plt.show()


### omni testing ###
if False:
    print ("testing omni module")
    start_date = datetime.datetime(2020,11,10,0,0,0,0)
    end_date = datetime.datetime(2020,12,10,0,0,0,0)
    download_last_omni()
    df_omni = get_omni_data(start_date,end_date,"Bz")
    #plt.plot(df_omni)
    #plt.show()

### sunspots testing ##
if False:
    print ("testing sunspots module")
    start_date = datetime.datetime(1818,1,1,0,0,0,0)
    end_date = datetime.datetime.now()
    download_sunspots()
    df_sunspots = get_sunspots(start_date,end_date)
    print (df_sunspots)
    #plt.plot(df_sunspots)
    #plt.show()

### sdo testing ##
if False:
    print ("testing sdo module")
    start_date = datetime.datetime(2020,11,11,0,0,0,0)
    end_date = datetime.datetime(2020,11,11,5,0,0,0)
    download_images(start_date,end_date,res=512)
    #make_video(start_date,end_date,channel="0171",res=512)



