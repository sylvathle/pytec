# pytec
Python library for Total Electron Content computing from Rinex files

## Installation

Go where you want pytec to be stored:
```bash
cd where/you/want/pytec/to/be/stored
```
Download the git repostory on your computer:
```
wget https://github.com/sylvathle/pytec.git
```

Clone it to get the folder structure:
```bash
git clone pytec.git
```

Remove the .git file
``` bash
rm pytec.git
```

Install pytec following instructions. You will be asked to set up a repository and an environment variable "PYTEC_PATH" so pytec can store intermediate the results. You can skip it but processes might take longer to compute:

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



