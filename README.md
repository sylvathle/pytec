# pytec
Python library for Total Electron Content computing from Rinex files

## Installation

Go where you want pytec to be stored:
```
cd where/you/want/pytec/to/be/stored
```
Download the git repostory on your computer:
```
wget https://github.com/sylvathle/pytec.git
```

Clone it to get the folder structure:
```
git clone pytec.git
```

Remove the .git file
``` 
rm pytec.git
```

Install pytec following instructions. You will be asked to set up a repository so pytec can store intermediate the results. You can skip it but processes might take longer to compute:

```
cd pytec
./install-devmode-package.sh
```

## Paths
Source and target files must be hosted in a folder which path is contained in a environment variable called RINEX\_PATH.
It is recommended to save the RINEX\_PATH in the .bashrc file so that there will be no need to redefine at each new terminal session.

This folder is structured according to year and day of year, for example:
RINEX\_PATH:
 - 2020
  --> 349
  --> 350
 - 2021
  --> 37
In each day folder must be hosted the RINEX observation files and will be written feather files containing the computed TEC (same name, but extension changed from statDOY0.YYo to stat.feather where stat stands for the 4 letter designing the GNSS station). 
In order to work, this folder must also contain the navigation rinex files to get satellite ephemeris information.

## Code

Here is a simple code:
```bash
f_rinex="position of the rinex file"
# Instanciate a tec object from a f_rinex observation file
tec_station = tec.tec(f_rinex)
# Compute STEC from the rinex file
tec_station.rinex_to_stec(60)
# Compute the satellite position from navigation file
tec_station.add_satellite_pos()
# Compute TEC baseline
tec_station.add_baseline(f_bias="./example/P1P22203.DCB")
# Compute station biais with method of maximazing correlation between TEC from each satellites from the station
tec_station.add_receiver_bias()
# Store TEC data in feather format
tec_station.to_feather(f_feather)
```



