# Input VPR
This repository contains the module for create the input file for VPR model from INGV, Roma.

There is a Parent Class which is called: ***ingestion***. Also, there are 2 Childs Classes called: ***MODIS, SLSTR***.
In other hand, for SEVIRI there is an unique Class called: ***SEVIRI***

## Structure 
The module *invpr.py* contains the classes and functions that build the Input file for VPR model. In this moment, the module allows to create the input file for 3 satellites data: **MODIS, SLSTR** and **SEVIRI**. 

The detailed structure of module is:

- Libraries used

- Classes:
	- ingestion:
		- MODIS()
		- SLSTR()
	- SEVIRI()

- Functions:
	- Used for the class MODIS
	- Used for the class SLSTR
	- Used for the class SEVIRI 

## Quick start guide

Firstly, the module needs one directory which must include:
- Satellite mages data (According to each satellite)
- Mask Files 
- paramFile.txt (With the information about the height and temperature of the volcanic ash cloud)

The *paramFile.txt* file must have the follow structure (only can change the value for each parameter):
```
name;value
t_ash(C);-51.326544
t_so2(C);-51.326544
h_ash(Km);10.54
h_so2(Km);10.54 
```

Then, for create the Input VPR file, the steps are like follow:

1) Create a object using some of the 3 classes: 
```
in_vpr = MODIS()
```

2) Set the **Path** of Directory
```
path = 'E:/Test VPR/MODIS/'
```

3) Create the Input VPR file using the function *generator()*. The output file will be storage in the same path. 
```
in_vpr.generator(path)
```