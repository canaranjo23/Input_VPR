# Input VPR
This repository contains the module to create the input file for the VPR model of the INGV, Roma.

There is a Parent Class which is called: ***ingestion***. Also, there are 2 Childs Classes called: ***MODIS, SLSTR***.
On the other hand, for SEVIRI there is an unique Class called: ***SEVIRI***

## Structure 
The module *invpr.py* contains the classes and functions to build the Input file for VPR model. In this moment, the module allows to create the input file for 3 satellites data: **MODIS, SLSTR** and **SEVIRI**. 

The detailed structure of module is:

- Libraries used

- Classes:
	- ingestion:
		- MODIS()
		- SLSTR()
	- SEVIRI()

- Functions:
	- Used for the MODIS class 
	- Used for the SLSTR class 
	- Used for the SEVIRI class 

## Quick start guide

Firstly, the module needs a directory which must include:
- Satellite images (According to each satellite)
- Mask Files 
- paramFile.txt (With information about height and temperature of the volcanic ash cloud)

The *paramFile.txt* file must have the follow structure (only can change the value for each parameter):
```
name;value
t_ash(C);-51.326544
t_so2(C);-51.326544
h_ash(Km);10.54
h_so2(Km);10.54 
```

Then, to create the Input VPR file, the steps are like follow:

1) Creating an object using some of the 3 classes available: 
```
in_vpr = MODIS()
```

2) Set the **Path** of the root directory
```
path = 'E:/Test VPR/MODIS/'
```

3) Creating the Input VPR file using the function *generator()*. The output file will be storage in the same path. 
```
in_vpr.generator(path)
```
