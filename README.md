# Input_VPR
This repository contains the module for create the input file for VPR model from INGV, Roma.

There is a Parent Class which is called: _*Ingestion*_. Also, there are 2 Childs Classes called: _*MODIS, SLSTR*_.
In other hand, for SEVIRI there is an unique Class called: _*SEVIRI*_

## Structure 
The module *invpr.py* contains the classes and functions that build the Input file for VPR model. In this moment, the module allows to create the input file for 3 satellites data: **MODIS, SLSTR and SEVIRI**. 

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

