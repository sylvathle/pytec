# pytec
Python library for Total Electron Content computing from Rinex files

## Requirements

Tested on python 3.9.2

## Installation

Go where you want pytec to be stored:
```bash
cd /where/you/want/pytec/to/be/stored
```
Clone the git repostory on your computer:
```
git clone https://github.com/sylvathle/pytec.git
```

Install pytec following instructions. You will be asked to set up a repository and an environment variable "PYTEC_PATH" so pytec can store intermediate the results. You can skip it but VTEC might take longer to compute:

```bash
cd pytec
./install-devmode-package.sh
```

### Known Errors

**Incompatible xarray module version with georinex**

Georinex has some specific dependencies that we have observed could be broken.
In particular latest versions of the "xarray" module are not compatible with georinex and might make the georinex "load" method fail.

A known solution is to force a working version of "xarray":

```bash
pip3 install xarray=0.20.1
```

**Other problems**

Please report problems you find to install pytec to *sylvain.blunier@gmail.com*

## Simple example



## Code

### Loading module

To load the pytec module do:

`from pytec import tec`

### Defining input data

Compute VTEC requires to have the following elements:

1. A list Rinex observation file
This is the file list where the code and phase data will be extracted. The files will be concatenated in the process they should come in continous and chronological order.
* Format: List of string
* Example: `f_obs=["/home/me/pytec/amtc1200.21o","/home/me/pytec/amtc1210.21o","/home/me/pytec/amtc1220.21o"]`

2. List of Rinex navigation file
This is the file list where the position of the satellites will be extracted, usually with a few well chosen files, positions of all satellites can be calculated. Of course, they should correspond to the periods of observation files.
* Format: List of string
* Example: `f_nav=["/home/me/pytec/amtc1201.21n","/home/me/pytec/amtc1211.21n","/home/me/pytec/amtc1221.21n"]`

3. File containing bias of satellites in DCB
File containing the bias of satellites corresponding to the month of observation files
* Format: string pointing to DCB file
* Example: `f_bias="/home/me/pytec/P1P22104.DCB"`

### Computing VTEC

Once these elements are defined, an instance of pytec can be created:

`station_tec = tec.tec(f_obs,f_nav_f_bias)`

And VTEC can be computed with default parameters:

`station_tec.compute_vtec()`

This will create one file in *feather* format in the same path as the second item of `f_obs`, or the first one in case `f_obs` is a one item list, so in the case of the variables defined upper the files should be found at `/home/me/pytec/amtc121-tec.feater`.
It will contain the following information:

* **time** (string): *YYYY-MM-DD HH:mm:ss* format
* **sv** (string): satellite from which is calculated VTEC
* **lat** (float): latitude of calculted piercing point
* **lon** (float):  longitude of calculated piercing point
* **elevation** (float): elevation of satellite in Radians
* **STEC_slp** (float): slant TEC in TEC Units calculated with pseudo range (P1 and C1 in rinex observation files)
* **STEC_sl** (float): slant TEC in TEC Units calculated with code phases (L1 and L2 in rinex observation files)
* **VTEC** (float): Vertical TEC in TEC Units, corrected with baseline, satellite bias, receiver bias.

The pandas DataFrame of this data can be directly accessed in the instanciated object through: `station_tec.df_obs`.

### Plots

Fast PNGs of data can be directly produced using the `pytec.graph` submodule:

```
from pytec import tec
```

All satellites in one graph
```
graph.plot_station(station_tec.df_obs,destination)
```

or 8 satellited per graph in a mozaic fashion
```
graph.plot_station(station_tec.df_obs,destination,mozaic=True)
```

Of course the `destination` argument should contain the path where you wish to keep the produced graphs, do not add any extension, the method does it alone. 


