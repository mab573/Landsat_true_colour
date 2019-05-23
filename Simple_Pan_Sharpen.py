#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# SCRIPT NAME		:: Simple_Pan_Sharpen.py 

# PURPOSE 		:: 
#				   

# SYNOPSIS 		::  
#				:: 
#				:: 
#				::  

# MODULES CALLED 	:: 
#					:: 
#					:: 
#					:: 
#
# SUBROUTINES		:: 
#					::  
#					:: 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


from osgeo import gdal
from osgeo import osr
import matplotlib as plt
from pylab import *
import numpy as np
import getopt, sys, math, os, glob, re, subprocess
import time
from pyproj import Proj
import scipy as sci 

# My/local subroutines

# Don't use congrid, use gdal instead
## TODO write code to do resampling using GDAL
#import congrid #copied from the interweb


def Simple_Pan_Sharpen(blue_band, green_band, red_band, pan_band):


    blue_max=np.max(blue_band)
    green_max=np.max(green_band)
    red_max=np.max(red_band)	

    # Letting GDAL increase resolution so there is no need to do congrid here
    new_blue_band=blue_band
    new_green_band=green_band
    new_red_band=red_band	

    #new_blue_band=congrid.congrid(blue_band, pan_band.shape, method='neighbour', centre=True, minusone=False)
    #new_green_band=congrid.congrid(green_band, pan_band.shape, method='neighbour', centre=True, minusone=False)
    #new_red_band=congrid.congrid(red_band, pan_band.shape, method='neighbour', centre=True, minusone=False)

    # Remove any NaN values if congrid does wierd stuff
    whereAreNaNs = np.isnan(new_blue_band);
    new_blue_band[whereAreNaNs] = 0;
	
    whereAreNaNs = np.isnan(new_green_band);
    new_green_band[whereAreNaNs] = 0;
	
    whereAreNaNs = np.isnan(new_red_band);
    new_red_band[whereAreNaNs] = 0;
	
    dims=pan_band.shape
	
    tmp = pan_band[0,0]
    dt = tmp.dtype

    new_blue_band_pan = np.zeros((dims[0],dims[1]),dtype = dt)
    new_green_band_pan = np.zeros((dims[0],dims[1]),dtype = dt)
    new_red_band_pan = np.zeros((dims[0],dims[1]),dtype = dt)
    new_pan_band = np.asarray(pan_band,dtype = dt)
	

    ## Pan sharpening method from (https://github.com/gina-alaska/dans-gdal-scripts/wiki/Gdal_landsat_pansharp)
    # Original
    #ratio_b=0.25
    #ratio_g=0.23
    #ratio_r=0.52

    ratio_b=0.25
    ratio_g=0.25
    ratio_r=0.25

    new_blue_band_pan=new_blue_band*(pan_band/((new_blue_band*ratio_b) + (new_green_band*ratio_g) + (new_red_band*ratio_r)))
    new_green_band_pan=new_green_band*(pan_band/((new_blue_band*ratio_b) + (new_green_band*ratio_g) + (new_red_band*ratio_r)))
    new_red_band_pan=new_red_band*(pan_band/((new_blue_band*ratio_b) + (new_green_band*ratio_g) + (new_red_band*ratio_r)))	

    # Once again remove any NaN values incase there is a divide by zero or something going on here
    whereAreNaNs = np.isnan(new_blue_band_pan);
    new_blue_band_pan[whereAreNaNs] = 0;

    whereAreNaNs = np.isnan(new_green_band_pan);
    new_green_band_pan[whereAreNaNs] = 0;

    whereAreNaNs = np.isnan(new_red_band_pan);
    new_red_band_pan[whereAreNaNs] = 0;

    return {'blue': new_blue_band_pan, 'green': new_green_band_pan, 'red': new_red_band_pan}

