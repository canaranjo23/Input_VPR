# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 12:16:58 2022
INGV, Rome - Italy

This module contains the classes and functions that build the Input file for 
VPR model (INGV). 

The Parent Class is called: Ingestion 
Also, there are 2 Childs Classes called: MODIS, SLSTR. 
In other hand, for SEVIRI there is an unique Class called: SEVIRI

@author: Camilo Naranjo
"""

# %% - Libraries
import numpy as np
from osgeo import gdal
import os
from pyhdf.SD import SD, SDC
import datetime
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.interpolate import griddata
import time
import shutil 

# %% - Classes

# ============================================================================
    # Parent Class
# ============================================================================
class ingestion: 
    """Primary class for the Input VPR file generation"""
    
    def __init__(self):
        """It creates the all essential variables."""
        self.name = ''
        self.path = ''
        self.d = {}             # Dictionary for save temperature and height data
        self.classes = {'No-Plume':0, 'SO2':1, 'Vol':2, 'And':3, 'Obs':4, 'Ice':5, 'H2O':6, 'SiO':7} 
        self.b85um = []         # 8.5um Band
        self.b11um = []         # 11um Band
        self.b12um = []         # 12um Band
        self.lat = []           # Latitude Data
        self.lon = []           # Longitude Data
        self.vza = []           # View Zenith Angle (Degrees)
        self.pmask = []         # Plume mask
        self.temp_ash = []      # Plume mean temperature particiles (°C)
        self.temp_so2 = []      # Plume mean temperature S02 (°C)
        self.h_ash = []         # Plume mean height particles (Km)
        self.h_so2 = []         # Plume mean height SO2 (Km)
        self.pixel_area = []    # Pixel area (Km2)
        self.lsm = []           # Land-Sea Mask
        self.cm = []            # Cloud Mask
        
    def read_pmask(self):
        """This function read the Plume Mask. 
        The plume mask in this point must be ready with the classification (0/7)
        Formats: .dat / .img"""
        
        # Let's create the default Plume Mask -> 'No plume = 0'
        self.pmask = np.zeros((self.b11um).shape)
        
        # Search in the path automatically
        os.chdir(self.path)
        ls = os.listdir(self.path)
        chararray = np.array(ls)
        
        # Conditions to build the Plume mask
        #   -SO2- [1]
        if query(chararray, 'S02') > 0:
            self.get_mask(chararray, 'S02')
        #   -Vol- [2]
        if query(chararray, 'Vol') > 0:
            self.get_mask(chararray, 'Vol')
        #   -And- [3]
        if query(chararray, 'And') > 0:
            self.get_mask(chararray, 'And')
        #   -Obs- [4]
        if query(chararray, 'Obs') > 0:
            self.get_mask(chararray, 'Obs')
        #   -Ice- [5]
        if query(chararray, 'Ice') > 0:
            self.get_mask(chararray, 'Ice')
        #   -H2O- [6]
        if query(chararray, 'H2O') > 0:
            self.get_mask(chararray, 'H2O')
        #   -SiO- [7]
        if query(chararray, 'SiO') > 0:
            self.get_mask(chararray, 'SiO')
        
    def get_mask(self, m, cat):
        """This function query if there is a file mask in the directory, then 
        get the filename and finally save the data"""
    
        # Get the filename
        index_array = np.char.find(m, cat)
        bool_buffer = np.ones(m.shape, dtype=bool)
        bool_buffer[index_array == -1] = False
        bl = np.logical_and(bool_buffer, np.char.endswith(m,'.hdr'))
        filename = m[bl][0][:-4]
        
        # Read the mask
        datapm = gdal.Open(self.path+filename, gdal.GA_ReadOnly) 
        pmask = datapm.GetRasterBand(1) 
        pmask = pmask.ReadAsArray()
        
        self.pmask[pmask > 0] = self.classes[cat]
        
    def read_temperature(self):
        """This function read the plume Temperature in °C, previously 
        calculated. 
        Format: Float (xx.xxxx)"""
        
        t_ash = np.float(self.d['t_ash(C)'])
        t_so2 = np.float(self.d['t_so2(C)'])

        self.temp_ash = (self.pmask > 0)*t_ash
        self.temp_so2 = (self.pmask > 0)*t_so2
        
    def read_height(self):
        """This function read the plume Height in KM, previously calculated. 
        Format: Float (xx.xxxx)"""
        
        h_ash = np.float(self.d['h_ash(Km)'])
        h_so2 = np.float(self.d['h_so2(Km)'])
        
        self.h_ash = (self.pmask > 0)*h_ash
        self.h_so2 = (self.pmask > 0)*h_so2
        
    def stack(self):    
        """This function stack all variables in an numpy array. This array is 
        neccesary for create the Input VPR file."""
        stacked_variables = np.stack((self.b85um, 
                          self.b11um,
                          self.b12um,
                          self.lat,
                          self.lon,
                          self.vza,
                          self.pmask,
                          self.temp_ash,
                          self.temp_so2,
                          self.h_ash,
                          self.h_so2,
                          self.pixel_area,
                          self.lsm,
                          self.cm), axis=0)     
        return stacked_variables
    
    def save(self):
        """
        This routine allow save a stack of array to ENVI file format.    
            > Input: 
                array = Stacked array with the information for each band.          
            < Output:
                file = ENVI file and .hdr file will are create and storage  
        """
        
        # Call to stack function.
        array = self.stack()
        
        band_names = ['IR_8.5 (W/(m2*sr*micron)',
                      'TIR_11.0 (W/(m2*sr*micron)',
                      'TIR_12.0 (W/(m2*sr*micron)',
                      'Latitude (deg)',
                      'Longitude (deg)',
                      'View Zenith Angle (deg)',
                      'Plume Mask',
                      'Plume Mean Temperature Particles (C)',
                      'Plume Mean Temperature SO2 (C)',
                      'Plume Mean Height Particles (km)',
                      'Plume Mean Height SO2 (km)',
                      'Pixel Area (km2)',
                      'Land-Sea Mask',
                      'Cloud Mask']
        
        # Let us define the ENVI format
        driver = gdal.GetDriverByName('ENVI')
        
        # Let us create one OBJECT where storage all data for the file.
        outRaster = driver.Create(self.path + self.name, 
                                  array.shape[2], array.shape[1], array.shape[0], 
                                  gdal.GDT_Float32)
        
        # Let us storage each array data into its respective band
        for i, image in enumerate(array[:,:], 1):
            b = outRaster.GetRasterBand(i)
            b.SetDescription(band_names[i-1])
            b.WriteArray( image )
        outRaster.FlushCache()
        
        print('Input VPR file created and storaged sucessfully!')
        
    def plot(self, vble, c="gray"):
        """This function allow to plot fastly some attribute from Input VPR file"""
        fig = plt.figure(dpi=900)
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.imshow(vble, cmap=c)

# ============================================================================
    # MODIS Class
# ============================================================================
class MODIS(ingestion):
    """Derived class for MODIS data.
    The MODIS data is requeried in .hdf format. Web-page information:
        https://ladsweb.modaps.eosdis.nasa.gov/search/
        https://lpdaac.usgs.gov/data/get-started-data/collection-overview/missions/modis-overview/
        """
    
    def generator(self, path):
        """This function allows reading all data from the .hdf file. 
        From here is extracted: 8.5um Band, 11um Band, 12um Band, Latitude, 
        Longitude, View Zenith Angle(VZA), Land-Sea Mask. Also is calculated 
        the Pixel Area (using a previous function) and the Cloud Mask."""
        
        self.path = path
        
        # Search in the path automatically
        os.chdir(path)
        ls = os.listdir(path)
        chararray = np.array(ls)
        
        # Get the filename
        bl = np.logical_and(np.char.startswith(chararray,'MOD021KM.'), np.char.endswith(chararray,'.hdf'))
        filename = chararray[bl][0]
        
        # Get the geo filename 
        bl = np.logical_and(np.char.startswith(chararray,'MOD03.'), np.char.endswith(chararray,'.hdf'))
        geo_filename = chararray[bl][0]
        
        # Let's manage the date and day of image.
        info = filename.split('.')      # Split the all information from filename
        jdate = info[1][1:]             # Get the julian date
        date = jd2stdd(jdate)           # Convert the julian date to standar date format
        
        # Identify the MODIS Platform 
        if info[0][:3] == 'MOD':
            sat = 'Terr'
        elif info[0][:3] == 'MYD':
            sat = 'Aqua'
        
        # Define the name of Input file VPR
        self.name = sat+'_'+(date.strftime('%Y%m%d'))+'_'+info[2]+'00'+'_VPRinput'
        pp(1, 'The name of file is: ' + self.name)
        
        # --------------------------------------------------------------------
        # Open file data
        hdf = SD( (path+filename) , SDC.READ)
        pp(0, 'Files .hdf opened successfully.')
        
            # Read datasets of my interest
        datafield_name = 'EV_1KM_Emissive'
        TIR_bandsOBJ = hdf.select(datafield_name)        # Contain the attributes and information
        pp(0, 'Data extracted successfully.')
        TIR_bands = calibration_TIR_MODIS(TIR_bandsOBJ)  # Calibration of Data
        
        # Extract important bands
        self.b85um = TIR_bands[:,:,8]
        self.b11um = TIR_bands[:,:,10]
        self.b12um = TIR_bands[:,:,11]
        
        # --------------------------------------------------------------------
        # Open geolocation data
        hdf_geo = SD( (path+geo_filename) , SDC.READ)
        
            # Read geolocation information from MOD03 product.
        latp = hdf_geo.select('Latitude')       # Get latitude information
        self.lat = latp[:,:]
        lonp = hdf_geo.select('Longitude')      # Get longitude information
        self.lon = lonp[:,:]
        
            # Get the information about View Zenith Angle (VZA)
        vza = hdf_geo.select('SensorZenith')
        scale_factor_g = 0.01
        self.vza = vza[:,:]*scale_factor_g
        
            # Get the Land-Sea mask
        lsmask = hdf_geo.select('Land/SeaMask')
        lsmp = lsmask[:,:]
        self.lsm = ((lsmp > 0) & (lsmp < 4))*1
        
            # Create the Cloud-Mask
        # MODIS Cloud Mask User Guide http://cimss.ssec.wisc.edu/modis/CMUSERSGUIDE.PDF
        self.cm = np.zeros(self.lsm.shape)
        
        # --------------------------------------------------------------------
        # Pixel area calculation 
        self.pixel_area = pixel_area_calculator()
        
        # --------------------------------------------------------------------
        # Read the Plume Mask
        self.read_pmask()
        
        # --------------------------------------------------------------------
        # Read data of Temperature and Height
        with open(self.path+"paramFile.txt") as f:
            next(f)
            for line in f:
                (key, val) = (line.rstrip('\n')).split(';')
                self.d[key] = val        
        
        self.read_temperature()
        self.read_height()
        
        # --------------------------------------------------------------------
        # Save input VPR file
        self.save()
    
# ============================================================================
    # SLSTR Class
# ============================================================================
class SLSTR(ingestion):
    """Derived class for SLSTR data
    https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-3-slstr/level-1/tir-brightness-temperature
    """
    
    def generator(self, path):
        
        self.path = path
        
        # Search in the path automatically
        os.chdir(path)
        ls = os.listdir(path)
        chararray = np.array(ls)
        
        # Get the filename
        bl = np.logical_and(np.char.startswith(chararray,'S3'), np.char.endswith(chararray,'.SEN3'))
        dirf = chararray[bl][0] + '/'
       
        # Let's read the bands and data
        b11um_ob = Dataset(path+dirf+'S8_BT_in.nc', 'r')
        b11um = np.array(b11um_ob['S8_BT_in'][:]) 
        
        b12um_ob = Dataset(path+dirf+'S9_BT_in.nc', 'r')
        b12um = np.array(b12um_ob['S9_BT_in'][:]) 
        
        # Manage the no data pixels
        nan_value = -32768.0                 # Value = NAN 
        b11um[b11um == nan_value] = np.NAN
        b12um[b12um == nan_value] = np.NAN
        b85um = np.ones(b11um.shape)         # Create the 8.5um band (SLSTR = 0.0)
        
        # Convert the BT bands to Radiance
        self.b11um = rad(11, b11um)
        self.b12um = rad(12, b12um)
        self.b85um = b85um*0.1
        
        # -------------------------------------------------------------------
        # Read the Name and Date of image. 
        with open(path+dirf+'xfdumanifest.xml') as f:
            read_data = f.read()
            
            start = read_data.find('<sentinel-safe:startTime>')
            end = read_data.find('</sentinel-safe:startTime>')
            date = read_data[start+25:end]
            
            start_id = read_data.find('<sentinel3:productName>')
            sat_id = read_data[start_id+23:start_id+23+3]
        
        # Define the name of Input file VPR
   
            # Identify the MODIS Platform 
        if sat_id == 'S3A':
            sat = 'SS3A'
        elif sat_id == 'S3B':
            sat = 'SS3B'
        
            # Identify date and time
        dt = ''.join(date[:10].split('-'))
        time = date[10+1:].split(':')
        last = str(round(float(time[-1][:5])))
        time[-1] = last
        t = ''.join(time)
        self.name = sat+'_'+dt+'_'+t+'_VPRinput'
        
        # Get Latitude and Longitude data
        geoloc_ob = Dataset(path+dirf+'geodetic_in.nc', 'r')
        self.lat = np.array(geoloc_ob['latitude_in'][:])
        self.lon = np.array(geoloc_ob['longitude_in'][:])
        
        # Get the information about View Zenith Angle (VZA)
        geodetic_ob = Dataset(path+dirf+'geodetic_tx.nc', 'r')
        latx = np.array(geodetic_ob['latitude_tx'][:])
        lonx = np.array(geodetic_ob['longitude_tx'][:])
        
        latxp = latx.flatten()
        lonxp = lonx.flatten()
        
        geometric_ob = Dataset(path+dirf+'geometry_tn.nc', 'r')
        vza_t_ob = np.array(geometric_ob['sat_zenith_tn'][:])
        vza_t_ob = vza_t_ob.flatten()
        self.vza = griddata((latxp,lonxp), vza_t_ob, (self.lat, self.lon), method='linear')
        
        # Get the Land-Sea mask
        # https://sentinel.esa.int/documents/247904/4598082/Sentinel-3-SLSTR-Land-Handbook.pdf/
        lsm_ob = Dataset(path+dirf+'flags_in.nc', 'r')
        lsmask = np.array(lsm_ob['confidence_in'][:]) 
        
        lsmaskp = vec_bin_array(lsmask)             # Convert data to Binary
        self.lsm = lsmaskp[:,:,11]
        
        # Create the Cloud-Mask
        self.cm = np.zeros(self.lsm.shape)
        
        # Create the Pixel-Area 
        self.pixel_area = np.ones(self.lsm.shape)
        
        # --------------------------------------------------------------------
        # Read the Plume Mask
        self.read_pmask()
        
        # --------------------------------------------------------------------
        # Read data of Temperature and Height
        with open(self.path+"paramFile.txt") as f:
            next(f)
            for line in f:
                (key, val) = (line.rstrip('\n')).split(';')
                self.d[key] = val        
        
        self.read_temperature()
        self.read_height()
        
        # --------------------------------------------------------------------
        # Save input VPR file
        self.save()



# ============================================================================
    # SEVIRI Class
# ============================================================================
class SEVIRI: 
    """Unique class for the Input VPR file SEVIRI data. 
    Due to the structure and the format (NetCDF) of data from SEVIRI, it is 
    necessary change an little the generation of Input VPR file."""
    
    def __init__(self):
        """It creates the all essential variables."""
        self.name = ''
        self.path = ''
        self.d = {}             # Dictionary for save temperature and height data
        self.classes = {'No-Plume':0, 'SO2':1, 'Vol':2, 'And':3, 'Obs':4, 'Ice':5, 'H2O':6, 'SiO':7} 
        self.vza = []           # View Zenith angles in degrees
        self.vaa = []           # View Azimuth angles in degrees
        self.pmask = []         # Plume mask
        self.temp_ash = []      # Plume mean temperature particiles (°C)
        self.temp_so2 = []      # Plume mean temperature S02 (°C)
        self.h_ash = []         # Plume mean height particles (Km)
        self.h_so2 = []         # Plume mean height SO2 (Km)
        self.cm = []            # Cloud Mask
        self.pixel_area = []    # Pixel area (Km2)
        
    def generator(self, path):
        """Esta función permite leer el archivo que nos entrega Dario, lo 
        abrimos en modo lectura, extraemos las bandas y atributos necesarios
        para el calculo de VAA y VZA"""
        
        self.path = path
        
        # Search in the path automatically
        os.chdir(path)
        ls = os.listdir(path)
        chararray = np.array(ls)
        
        # Get the filename
        bl = np.char.endswith(chararray,'.nc')
        filename = chararray[bl][0]

        # Let's read the bands and data 
        data_ob = Dataset(path+filename, 'r')
        
        self.name = filename[:-3]
        
        # Get Latitude and Longitude data
        Lon = np.array(data_ob['longitude'][:])
        Lat = np.array(data_ob['latitude'][:])
        
        # Get attributes: Sat_Altitude, Sat_Lat, Sat_Long
        sat_height = getattr(data_ob['IR_108'], 'satellite_altitude')
        sat_lat = getattr(data_ob['IR_108'], 'satellite_latitude')
        sat_long = getattr(data_ob['IR_108'], 'satellite_longitude')
        resolution = getattr(data_ob['IR_108'], 'resolution')
        
        # VZA and VAA calculation
        self.vza, self.vaa = MSG_view_angles(Lat, Lon, sat_lat, sat_long, sat_height)
        
        # Create the Cloud-Mask
        self.cm = np.zeros(self.vza.shape)
        
        # Create the Pixel-Area 
        self.pixel_area = np.ones(self.vza.shape)*(resolution**2)/1000000
        
        # --------------------------------------------------------------------
        # Read the Plume Mask
        self.read_pmask()
        
        # --------------------------------------------------------------------
        # Read data of Temperature and Height
        with open(self.path+"paramFile.txt") as f:
            next(f)
            for line in f:
                (key, val) = (line.rstrip('\n')).split(';')
                self.d[key] = val        
        
        self.read_temperature()
        self.read_height()
        
        # --------------------------------------------------------------------
        # Save input VPR file
        self.save()
        
    def read_pmask(self):
        """This function read the Plume Mask. 
        The plume mask in this point must be ready with the classification (0/7)
        Formats: .dat / .img"""
        
        # Let's create the default Plume Mask -> 'No plume = 0'
        self.pmask = np.zeros(self.vza.shape)
        
        # Search in the path automatically
        os.chdir(self.path)
        ls = os.listdir(self.path)
        chararray = np.array(ls)
        
        # Conditions to build the Plume mask
        #   -SO2- [1]
        if query(chararray, 'S02') > 0:
            self.get_mask(chararray, 'S02')
        #   -Vol- [2]
        if query(chararray, 'Vol') > 0:
            self.get_mask(chararray, 'Vol')
        #   -And- [3]
        if query(chararray, 'And') > 0:
            self.get_mask(chararray, 'And')
        #   -Obs- [4]
        if query(chararray, 'Obs') > 0:
            self.get_mask(chararray, 'Obs')
        #   -Ice- [5]
        if query(chararray, 'Ice') > 0:
            self.get_mask(chararray, 'Ice')
        #   -H2O- [6]
        if query(chararray, 'H2O') > 0:
            self.get_mask(chararray, 'H2O')
        #   -SiO- [7]
        if query(chararray, 'SiO') > 0:
            self.get_mask(chararray, 'SiO')
        
    def get_mask(self, m, cat):
        """This function query if there is a file mask in the directory, then 
        get the filename and finally save the data"""
    
        # Get the filename
        index_array = np.char.find(m, cat)
        bool_buffer = np.ones(m.shape, dtype=bool)
        bool_buffer[index_array == -1] = False
        bl = np.logical_and(bool_buffer, np.char.endswith(m,'.hdr'))
        filename = m[bl][0][:-4]
        
        # Read the mask
        datapm = gdal.Open(self.path+filename, gdal.GA_ReadOnly) 
        pmask = datapm.GetRasterBand(1) 
        pmask = pmask.ReadAsArray()
        
        self.pmask[pmask > 0] = self.classes[cat]
        
    def read_temperature(self):
        """This function read the plume Temperature in °C, previously 
        calculated. 
        Format: Float (xx.xxxx)"""
        
        t_ash = np.float(self.d['t_ash(C)'])
        t_so2 = np.float(self.d['t_so2(C)'])

        self.temp_ash = (self.pmask > 0)*t_ash
        self.temp_so2 = (self.pmask > 0)*t_so2
        
    def read_height(self):
        """This function read the plume Height in KM, previously calculated. 
        Format: Float (xx.xxxx)"""
        
        h_ash = np.float(self.d['h_ash(Km)'])
        h_so2 = np.float(self.d['h_so2(Km)'])
        
        self.h_ash = (self.pmask > 0)*h_ash
        self.h_so2 = (self.pmask > 0)*h_so2
        
    def save(self):
        """
        """
        
        # Copy file and rename
        shutil.copyfile(self.path+self.name+'.nc', self.path+self.name+'_VPR.nc')
        
        # Let's read the data (File readed in mode Read and Write)
        # Opening a file, creating a new Dataset
        try: data.close()  # just to be safe, make sure dataset is not already open.
        except: pass
        data = Dataset(self.path+self.name+'_VPR.nc', 'r+', format='NETCDF4')
        
        # # Opening a file, creating a new Dataset
        # try: data.close()  # just to be safe, make sure dataset is not already open.
        # except: pass
        # data = Dataset('nevvvv.nc', 'w', format='NETCDF4')
        
        # # Creating dimensions
        # x = data.createDimension('x', 672)
        # y = data.createDimension('y', 464)

        # Creating nc attributes
        data.title = 'Input VPR data'
        data.subtitle = 'Complete file .nc input VPR data SEVIRI'
        data.history = "Created on " + time.ctime(time.time())
        
        # # Creating variables
        tvza = data.createVariable('VZA', np.dtype('float64'), ('y','x'))
        tvza.long_name = 'View Zenith Angle'
        tvza.units = 'Deg'
        
        tvaa = data.createVariable('VAA', np.float64, ('y','x'))
        tpmask = data.createVariable('plume_mask', np.int8, ('y','x'))
        tpmask.description = '0=no_plume; 1=SO2; 2=SO2+Vol; 3=SO2+And; 4=SO2+Obs; 5=SO2+Ice; 6=SO2+H2O; 7=SO2+SiO'
        
        tplume_t = data.createVariable('plume_t', np.float64, ('y','x'))
        tplume_t.setncattr('Profile used', 'NCEP_2021102312_37.5N_15.0E.txt')
        
        tplume_t_so2 = data.createVariable('plume_t_so2', np.float64, ('y','x'))
        tplume_h = data.createVariable('plume_h', np.float64, ('y','x'))
        tplume_h_so2 = data.createVariable('plume_h_so2', np.float64, ('y','x'))
        tcloud_mask = data.createVariable('cloud_mask', np.int8, ('y','x'))
        tpixel_area = data.createVariable('pixel_area', np.float64, ('y','x'))
        
        # Writing data 
        tvza[:,:] = self.vza
        tvaa[:,:] = self.vaa
        tpmask[:,:] = self.pmask
        tplume_t[:,:] = self.temp_ash
        tplume_t_so2[:,:] = self.temp_so2
        tplume_h[:,:] = self.h_ash
        tplume_h_so2[:,:] = self.h_so2
        tcloud_mask[:,:] = self.cm
        tpixel_area[:,:] = self.pixel_area
        
        data.close()
        
    def plot(self, vble, c="gray"):
        """This function allow to plot fastly some attribute from Input VPR file"""
        fig = plt.figure(dpi=900)
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.imshow(vble, cmap=c)

# %% - Functions

# ============================================================================
    # General
# ============================================================================    
def pp(i, msg):
    """This function aims to give format to print messages. 
    0 = Message without superior space 
    1 = Message with superior space """
    if i == 1:
        print(' ')
    print('+ ' + msg)
    
def query(m, cat):
    """This function makes the query for know if in the directory there is 
    the file mask according the dictionary. This step is neccesary for built 
    the total Plume Mask.   
    
    Return: 0 = Empty
            >1 = There is a file
    """
    
    index_array = np.char.find(m, cat)
    bool_buffer = np.ones(m.shape, dtype=bool)
    bool_buffer[index_array == -1] = False
    
    return np.sum(bool_buffer)

# ============================================================================
    # For MODIS 
# ============================================================================
def jd2stdd(jdate):
    """Convert Julian date to Standar date"""
    fmt = '%Y%j'
    datestd = datetime.datetime.strptime(jdate, fmt).date()
    return(datestd)

def calibration_TIR_MODIS(TIR_bands2D):
    """This function calibrate the MODIS images"""
        # Retrieve attributes.
    attrs = TIR_bands2D.attributes(full=1)
    
    lna=attrs["long_name"]
    long_name = lna[0]
    
    aoa=attrs["radiance_offsets"]
    add_offset = aoa[0]
    add_offset = np.array(add_offset)
    
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    
    sfa=attrs["radiance_scales"]
    scale_factor = sfa[0]   
    scale_factor = np.array(scale_factor)
         
    vra=attrs["valid_range"]
    valid_min = vra[0][0]        
    valid_max = vra[0][1]  
          
    ua=attrs["radiance_units"]
    units = ua[0]
    
    # Bands data
    TIR_bands = TIR_bands2D[:,:].astype(np.double)  # Contain data
    
    # Calibration of data
        # Correction of Invalid values
    invalid = np.logical_or(TIR_bands > valid_max, TIR_bands < valid_min)
    invalid = np.logical_or(invalid, TIR_bands == _FillValue)
    TIR_bands[invalid] = np.nan
    
        # Calibration
    TIR_bands = np.rot90(np.rot90(TIR_bands,k=1,axes=(2,1)),k=1,axes=(0,2))
    TIR_bands = (TIR_bands - add_offset) * scale_factor 
    TIR_bands = np.ma.masked_array(TIR_bands, np.isnan(TIR_bands))
    
    return TIR_bands

def pixel_area_calculator():
    """
    This function allow to compute the Pixel Area for MODIS 1KM data.
    
    Size data: 2030 x 1354 
    
    Procedure taken from: 
    https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=2018
    """
    # Parameters for MODIS 1KM
    hp =  676.5 # is 1/2 the total number of pixels (zero-based)
    H = 700     # is the sensor altitude divided by the pixel size (700km/1km)
    I = np.arange(1354)
    
    # Compute the scan angle, S (in radians), given pixel number:
    S = (I-hp)/H
    
    # Compute the zenith angle:
    Z = np.arcsin(1.111*np.sin(S)) 
    
    # Compute the Along-track pixel size:
    Pn = 1    # is the nadir pixel size (1km)
    Pt = Pn*9*np.sin(Z-S)/np.sin(S)
    
    # Compute the Along-scan pixel size:
    Ps = Pt/np.cos(Z)
    
    # Thus, the pixel area is ~ 
    pixel_area = Pt * Ps
    pixel_area = np.tile(pixel_area, (2030, 1))
    
    return pixel_area

# ============================================================================
    # For SLSTR
# ============================================================================
def vec_bin_array(arr, m = 15):
    """
    This function converts decimal array to binary vector array. 
    
    Arguments: 
    arr: Numpy array of positive integers
    m: Number of bits of each integer to retain

    Returns a copy of arr with every element replaced with a bit vector.
    Bits encoded as int8's.
    """
    to_str_func = np.vectorize(lambda x: np.binary_repr(x).zfill(m))
    strs = to_str_func(arr)
    ret = np.zeros(list(arr.shape) + [m], dtype=np.int8)
    for bit_ix in range(0, m):
        fetch_bit_func = np.vectorize(lambda x: x[bit_ix] == '1')
        ret[...,bit_ix] = fetch_bit_func(strs).astype("int8")

    return ret 

def rad(b, tb):
    """This function allows to calculate the radiance by a wavelength and BT given. 
    The Planck Function."""
    # Constants
    c1 = 1.191042e8         # (W m-2 sr-1 um4)
    c2 = 1.4387752e4        # (um K)
    
    # Seleccionamos la longitud de onda correspondiente a la banda. 
    # https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-3-slstr/level-1/observation-mode-desc
    if (b == 11):
        lmda = 10.85
    elif (b == 12):
        lmda = 12.0
    
    radiance = c1/((lmda**5)*(np.exp(c2/(lmda*tb))-1))
    return radiance


# ============================================================================
    # For SEVIRI
# ============================================================================
def EvalGeoDLat(Latitude): 
    '''Converts geographical to geodetic latitude

    INPUT
        Latitude: Latitude (deg)

    OUTPUT
        GeoLat: geodetic latitude (deg)

    MODIFICATION HISTORY:
        Written by:     Marco Clerici, 27.02.04 - Original MATLAB
        Rewritten by:   Camilo Naranjo, 29.03.2022 - PYTHON
    '''
    
    EARTH_RADIUS_M    =     6378.214
    RPO               =     6356.829
    DOUBLE_EPSILON    =     1.0E-10

    RadEq = (EARTH_RADIUS_M - RPO )/EARTH_RADIUS_M
    RadEq = (1 - RadEq)*(1 - RadEq)
    
    LatRad = Latitude * (np.pi/180)
    
    Tmpout = np.arctan(np.tan(LatRad)*RadEq)
    
    Diff = abs(2*np.pi - abs(LatRad))
    
    # Original Version - If conditional
    # if (Diff < DOUBLE_EPSILON):
    #     Tmpout = LatRad
    # else:
    #     RadEq = (1 - RadEq)*(1 - RadEq)
    #     Tmpout = np.arctan(np.tan(LatRad)*RadEq)
    #     GeoLat = Tmpout / (np.pi/180)
        
    GeoLat = np.zeros(Latitude.shape)
    
    GeoLat[(Diff < DOUBLE_EPSILON)] = LatRad[(Diff < DOUBLE_EPSILON)] / (np.pi/180)
    GeoLat[(Diff >= DOUBLE_EPSILON)] = Tmpout[(Diff >= DOUBLE_EPSILON)] / (np.pi/180)
    
    return GeoLat


def PixRadius(Latitude):
    '''PURPOSE:
        Calculate the radius of the Earth at a given latitude

     INPUT
         Latitude: Latitude (deg)

     OUTPUT
         Radius: radius of the Eath at the current latitude (km)

     MODIFICATION HISTORY:
         Written by:     Marco Clerici, 27.02.04 - Original MATLAB
         Rewritten by:   Camilo Naranjo, 29.03.2022 - PYTHON
    '''
         
    EARTH_RADIUS_M    =     6378.214;
    RPO               =     6356.829;

    Exc     = np.sqrt((EARTH_RADIUS_M * EARTH_RADIUS_M) - (RPO * RPO))/(EARTH_RADIUS_M)

    Exc2    = Exc*Exc
    Exc4    = Exc2*Exc2

    LatRad  = Latitude * np.pi/180
    Sp1     = np.sin(LatRad)
    Sp12    = Sp1*Sp1
    Sp14    = Sp12*Sp12

    Radius = EARTH_RADIUS_M * (1. - (Exc2*Sp12)/2. + (Exc4*Sp12)/2. - (5.*Sp14*Exc4)/8.)
    return Radius


def MSG_view_angles(Lat, Lon, Sat_Lat=0, Sat_Lon=0, Sat_Height=42185):
    '''PURPOSE:
        Compute the satellite viewing angles with respect to the normal to a given
        geographical location

    INPUTS:

        Lat:        array of geographical latitude from -90 (South) to 90 (North)
        Lon:        array of geographical longitude from -90 (West) to 90 (East),
                    when REF_LON=0

    KEYWORD PARAMETERS:

        HEIGHT:     height of the location (km)
        SAT_LAT:    latitude of the SSPoint ( default is 0.)
        SAT_LON:    longitude of the SSPoint ( default is 0.)
        SAT_HEIGHT: height of the satellite  ( default is 42185. Km)
        VERBOSE:    print verbose output
        SILENT:     force silent behaviour

    OUTPUTS:

        VZA:    array of View Zenith angles in degrees  ( Zenith is 0. degrees )
        VAA:    array of View Azimuth angles in degrees ( reference is North,
                                                         clockwise )

    MODIFICATION HISTORY:
    Written by: Marco Clerici, 17.02.04 - Original MATLAB

    SWb modif --> 13/01/2005
    To be rigourous, one must calculate exactly the height of the
    satellite as a function of its position.
    SWe modif --> 13/01/2005
    
    Check Lorenzo 15/02/2022
    Rewritten by:   Camilo Naranjo, 29/03/2022 - PYTHON
        Comments: I have replaced all if conditionals by booleans operators and 
        indixes numpy array. 
    '''      
    
    # Variables definition 
    DOUBLE_EPSILON = 1.0E-10
    Height = 0
    VZA = np.zeros(Lat.shape)
    VAA = np.zeros(Lat.shape)
    
    # Handle default value Sat_Height and correction
    Sat_Height = Sat_Height/1000    # Convert Sat_Height to Km
    if (Sat_Height != 42185):
        SAT_EIGHT_M = Sat_Height + 6371
    else: 
        SAT_EIGHT_M = Sat_Height
    
    # Validation of size Lat and Lon
    if (Lat.shape != Lon.shape):
        raise Exception('ERROR! Lat and Lon array must have the same dimension')
        
    # Transform satellite latitude from geographical to geocentric
    GeoSatLonRad   = Sat_Lon * np.pi/180
    GeoDSatLat=EvalGeoDLat(Sat_Lat)
    GeoDSatLatRad = GeoDSatLat * np.pi/180
    
    # Terrestrial coordinates (meters) of satellite
    XSat = SAT_EIGHT_M * np.cos(GeoDSatLatRad)*np.cos(GeoSatLonRad)
    YSat = SAT_EIGHT_M * np.cos(GeoDSatLatRad)*np.sin(GeoSatLonRad)
    ZSat = SAT_EIGHT_M * np.sin(GeoDSatLatRad)
        
    # Convert from degrees to radians
    SatLatRad = Lat * np.pi/180
    PixLonRad = Lon * np.pi/180
    PixLatRad = Lat * np.pi/180
    
    # Calculate geocentric latitude of observation point
    GeoDLat = EvalGeoDLat(Lat)
    GeoDLatRad = GeoDLat * np.pi/180
    
    # Earth radius at geographical latitude of observation point
    Radius = PixRadius(Lat)
    PixHeight = Radius + Height
    
    # Terrestrial coordinates (meters) of observation point
    XPix = PixHeight * np.cos(GeoDLatRad)*np.cos(PixLonRad)
    YPix = PixHeight * np.cos(GeoDLatRad)*np.sin(PixLonRad)
    ZPix = PixHeight * np.sin(GeoDLatRad)
    
    # Components of a vector orthogonal to the surface of the earth geoid
    Dro = np.ones(Lat.shape)

    Dro[abs(Lon) == 90] = 0
    Dro[(Lon == 0)] = np.sqrt(XPix[(Lon == 0)]*XPix[(Lon == 0)] + YPix[(Lon == 0)]*YPix[(Lon == 0)])
    Dro[np.logical_and((abs(Lon) > 0), (abs(Lon) < 90))] = ZPix[np.logical_and((abs(Lon) > 0),(abs(Lon) < 90))]/np.tan(PixLatRad[np.logical_and((abs(Lon) > 0),(abs(Lon) < 90))])

    XGeo = Dro * np.cos(PixLonRad)
    YGeo = Dro * np.sin(PixLonRad)
    ZGeo = ZPix
    
    # Vector norms
    TGeo = np.sqrt((XGeo*XGeo) + (YGeo*YGeo) + (ZGeo*ZGeo))
    TSat = np.sqrt((XSat*XSat) + (YSat*YSat) + (ZSat*ZSat))
    TSatGeo = np.sqrt((XSat-XPix)*(XSat-XPix) + (YSat-YPix)*(YSat-YPix) + (ZSat-ZPix)*(ZSat-ZPix))
    
    # Geocentric angle between directions of satellite and observation point
    GeoSatDAngle = (XSat*XPix + YSat*YPix + ZSat*ZPix)/(TSat*PixHeight)
    
    GeoSatDAngle[np.logical_and((abs(1 - abs(GeoSatDAngle)) < DOUBLE_EPSILON), (GeoSatDAngle > 0.))] = 1
    GeoSatDAngle[np.logical_and((abs(1 - abs(GeoSatDAngle)) < DOUBLE_EPSILON), (GeoSatDAngle < 0.))] = -1
    
    GeoSatDAngle = np.arccos(GeoSatDAngle)/(np.pi/180)
    
    # Check if geocentric directions of satellite and observer are coincident
    iLabSat = np.zeros(GeoSatDAngle.shape)
    iLabSat[np.logical_or((GeoSatDAngle > 179.9), (GeoSatDAngle <= 0.))] = 1
    
    # ----------------------------------------------------------------------------------------------------------------------------
    # Satellite Zenith Angle (VZA)
    Arg2 = (XPix*(XSat-XPix) + YPix*(YSat-YPix) + ZPix*(ZSat-ZPix))/(PixHeight*TSatGeo)
    
    Arg2[np.logical_and((abs(1. - abs(Arg2)) < DOUBLE_EPSILON), (Arg2 > 0.))] = 1
    Arg2[np.logical_and((abs(1. - abs(Arg2)) < DOUBLE_EPSILON), (Arg2 < 0.))] = -1
    
    VZA = np.arccos(Arg2) / (np.pi/180);
    VZA = np.real(VZA)
    
    # ----------------------------------------------------------------------------------------------------------------------------
    # Satellite Azimuth Angle (VAA)
    
    # Define vector v0 in the meridional plane, perpendicular
    # to vertical of obervation point and pointing to north pole
    # [see that vertical direction is not necessarily geocentric direction]
    V90X = np.zeros(Lat.shape)
    V90Y = np.zeros(Lat.shape)
    V90Z = np.zeros(Lat.shape)
    
        # Conditional 1
    mbool = np.logical_and((abs(90. - abs(GeoDLat)) < DOUBLE_EPSILON), (GeoDLat < 0.))
    V90X[mbool] = (1)*PixHeight[mbool]*np.cos(PixLonRad[mbool])
    mbool = np.logical_and((abs(90. - abs(GeoDLat)) < DOUBLE_EPSILON), (GeoDLat < 0.))
    V90Y[mbool] = (1)*PixHeight[mbool]*np.sin(PixLonRad[mbool])
    
        # Conditional 2
    mbool = np.logical_and((abs(90. - abs(GeoDLat)) < DOUBLE_EPSILON), (GeoDLat >= 0.))
    V90X[mbool] = (-1)*PixHeight[mbool]*np.cos(PixLonRad[mbool])
    mbool = np.logical_and((abs(90. - abs(GeoDLat)) < DOUBLE_EPSILON), (GeoDLat >= 0.))
    V90Y[mbool] = (-1)*PixHeight[mbool]*np.sin(PixLonRad[mbool])
    
    V90Z[(abs(90. - abs(GeoDLat)) < DOUBLE_EPSILON)] = 0
    
        # Conditional 3
    mbool = (abs(90. - abs(GeoDLat)) >= DOUBLE_EPSILON)
    V90X[mbool] = -XGeo[mbool]*ZGeo[mbool]
    V90Y[mbool] = -YGeo[mbool]*ZGeo[mbool]
    V90Z[mbool] = (XGeo[mbool]*XGeo[mbool]) + (YGeo[mbool]*YGeo[mbool])
    
    V90 = np.sqrt((V90X*V90X) + (V90Y*V90Y) + (V90Z*V90Z))
    
    # Define vector Y0 orthogonal to the vectors V90 and Geo
    # and forming with them a 3-axis right system (Y0 points eastwards)
    Y0X = V90Y*ZGeo - YGeo*V90Z
    Y0Y = V90Z*XGeo - V90X*ZGeo
    Y0Z = V90X*YGeo - XGeo*V90Y
    Y0 = np.sqrt((Y0X*Y0X) + (Y0Y*Y0Y) + (Y0Z*Y0Z))
    
    # Define the vector SS pointing to the satellite at observation point
    SSX = XSat - XPix
    SSY = YSat - YPix
    SSZ = ZSat - ZPix
    SS = np.sqrt((SSX*SSX) + (SSY*SSY) + (SSZ*SSZ))
    
    # Calculate vector on the plane tangent to the observation point
    # and pointing to the satellite = Geo x SS x Geo
    
        # Conditional 1 - (iLabSat == 0) and (Prod > 0.) 
    TmpAr1 = YGeo*SSZ - ZGeo*SSY
    TmpAr2 = ZGeo*SSX - XGeo*SSZ
    TmpAr3 = XGeo*SSY - YGeo*SSX

    ProdX = TmpAr2*ZGeo - TmpAr3*YGeo
    ProdY = TmpAr3*XGeo - TmpAr1*ZGeo
    ProdZ = TmpAr1*YGeo - TmpAr2*XGeo

    Prod = np.sqrt((ProdX*ProdX) + (ProdY*ProdY) + (ProdZ*ProdZ))   
    
    # Calculate angles wrt vectors Y0 and V90
    # If angle wrt to y0 exceeds 90 deg => satellite azimuth > 180 deg
    C1Sat = (V90X*ProdX + V90Y*ProdY + V90Z*ProdZ)/(V90*Prod)
    C2Sat = (Y0X*ProdX + Y0Y*ProdY + Y0Z*ProdZ)/(Y0*Prod)
    
        # Conditions
            # C1SAT
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod > 0.)), np.logical_and((abs(1. - abs(C1Sat)) < DOUBLE_EPSILON), (C1Sat > 0.)))
    C1Sat[mbool] = 1  
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod > 0.)), np.logical_and((abs(1. - abs(C1Sat)) < DOUBLE_EPSILON), (C1Sat < 0.)))
    C1Sat[mbool] = -1
            
            # C2SAT
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod > 0.)), np.logical_and((abs(1. - abs(C2Sat)) < DOUBLE_EPSILON), (C2Sat > 0.)))
    C2Sat[mbool] = 1
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod > 0.)), np.logical_and((abs(1. - abs(C2Sat)) < DOUBLE_EPSILON), (C2Sat < 0.)))
    C2Sat[mbool] = -1
    
        # Conditional 2 - (iLabSat == 0) and (Prod <= 0.) 
        # Special case of exterior product = 0 <> SS and Geo have same
        #  direction It occurs if the satellite is at the zenith (or nadir).
        #  Its azimuth is undetermined. For geosynchronous satellites we adopt the same solution as for Sun
        #  Minor modifications will be needed for non-geosynchronous platforms
    
        # (ZPix >= ZSat) and (VZA(ii) <= 90.)
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod <= 0.)), np.logical_and((ZPix >= ZSat), (VZA <= 90.)))
    C1Sat[mbool] = -1 
    C2Sat[mbool] = 0
    
        # (ZPix >= ZSat) and (VZA(ii) > 90.)
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod <= 0.)), np.logical_and((ZPix >= ZSat), (VZA > 90.)))
    C1Sat[mbool] = +1 
    C2Sat[mbool] = 0
    
        # (ZPix < ZSat) and (VZA(ii) <= 90.)
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod <= 0.)), np.logical_and((ZPix < ZSat), (VZA <= 90.)))
    C1Sat[mbool] = +1 
    C2Sat[mbool] = 0
    
        # (ZPix < ZSat) and (VZA(ii) > 90.)
    mbool = np.logical_and( np.logical_and((iLabSat == 0), (Prod <= 0.)), np.logical_and((ZPix < ZSat), (VZA > 90.)))
    C1Sat[mbool] = -1 
    C2Sat[mbool] = 0
    
    # Calculating VAA
    Saph1 = np.arccos(C1Sat)/(np.pi/180)
    Saph2 = np.arccos(C2Sat)/(np.pi/180)
    
    TmpVal = np.zeros(Lat.shape)
    TmpVal[Saph2 <= 90.] = Saph1[Saph2 <= 90.]
    TmpVal[Saph2 > 90.] = 360. - Saph1[Saph2 > 90.] 
    
        # Conditional 3
    TmpVal[np.logical_and(iLabSat == 1,  np.logical_and((ZPix >= ZSat), (GeoSatDAngle < 90.)))] = 180
    TmpVal[np.logical_and(iLabSat == 1,  np.logical_and((ZPix >= ZSat), (GeoSatDAngle > 90.)))] = 0
    TmpVal[np.logical_and(iLabSat == 1,  np.logical_and((ZPix <  ZSat), (GeoSatDAngle < 90.)))] = 0
    TmpVal[np.logical_and(iLabSat == 1,  np.logical_and((ZPix <  ZSat), (GeoSatDAngle > 90.)))] = 180
    
    VAA = TmpVal
    VAA = np.real(VAA)
    
    return VZA, VAA
