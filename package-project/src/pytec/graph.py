import matplotlib.pyplot as plt
import pandas as pd
import sys


def plot_station(df_station,png_file_name):

    list_sats = df_station["sv"].unique()

    fig, axs = plt.subplots(4,2,figsize=(14,10),sharex=True,sharey=True)

    i,j=0,0
    n_sat = 0
    n=0
    N=0
    for sat in list_sats:
        i = int(n/2)
        j = n%2
        df_sat = df_station[df_station["sv"]==sat]
        axs[i,j].plot(df_sat["VTEC"],'b.')
        axsx = axs[i,j].twinx()
        axsx.plot(df_sat["elevation"]*180/3.1415926535,'r.',markersize=1,alpha=0.5)
        axsx.set_ylim([0,90])
        axs[i,j].set_title(sat)
        axs[i,j].grid(True)
        axs[i,j].set_ylabel("VTEC")
        axsx.set_ylabel("elevation")
        axsx.spines['right'].set_color('red')
        axsx.tick_params(axis='y', colors='red')
        axsx.yaxis.label.set_color('red')

        n = n+1
        n_sat = n_sat+1

        if n==8:
            n=0
            plt.savefig(png_file_name+"-"+str(N)+".png",bbox_inches='tight')
            plt.close()
            fig, axs = plt.subplots(4,2,figsize=(14,10),sharex=True,sharey=True)
            N = N+1
    if n!=0:
        plt.savefig(png_file_name+"-"+str(N)+".png",bbox_inches='tight')
        plt.close()
        
