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
import congrid #copied from the interweb


def Simple_Pan_Sharpen(blue_band, green_band, red_band, pan_band):


	blue_max=np.max(blue_band)
	green_max=np.max(green_band)
	red_max=np.max(red_band)	

        # Letting GDAL increase resolution so there is no need to do congrid here
        new_blue_band=blue_band
        new_green_band=green_band
        new_red_band=red_band	


	new_blue_band=congrid.congrid(blue_band, pan_band.shape, method='neighbour', centre=True, minusone=False)
	new_green_band=congrid.congrid(green_band, pan_band.shape, method='neighbour', centre=True, minusone=False)
	new_red_band=congrid.congrid(red_band, pan_band.shape, method='neighbour', centre=True, minusone=False)

	print 'Minimum after congrid ', np.min(new_blue_band), np.min(new_green_band), np.min(new_red_band)
        print 'Maximum after congrid ', np.max(new_blue_band), np.max(new_green_band), np.max(new_red_band)

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

	new_blue_band_pan=new_blue_band*(pan_band/((new_blue_band*0.25) + (new_green_band*0.23) + (new_red_band*0.52)))
	new_green_band_pan=new_green_band*(pan_band/((new_blue_band*0.25) + (new_green_band*0.23) + (new_red_band*0.52)))
	new_red_band_pan=new_red_band*(pan_band/((new_blue_band*0.25) + (new_green_band*0.23) + (new_red_band*0.52)))	

	# Once again remove any NaN values incase there is a divide by zero or something going on here
	whereAreNaNs = np.isnan(new_blue_band_pan);
        new_blue_band_pan[whereAreNaNs] = 0;

        whereAreNaNs = np.isnan(new_green_band_pan);
        new_green_band_pan[whereAreNaNs] = 0;

        whereAreNaNs = np.isnan(new_red_band_pan);
        new_red_band_pan[whereAreNaNs] = 0;
	
	#print pan_band[5000,5000],new_blue_band[5000,5000],new_green_band[5000,5000],new_red_band[5000,5000]

        print 'Minimum after pan values ', np.min(new_blue_band_pan), np.min(new_green_band_pan), np.min(new_red_band_pan)
        print 'Maximum after pan values ', np.max(new_blue_band_pan), np.max(new_green_band_pan), np.max(new_red_band_pan)

	## This scales the images based on the minimum and maximum value of the image when it has been resampled
	## to about 250 m resolution. This should remove skewing of the result to single bright or dark pixels.

	#con_arr=congrid.congrid(new_blue_band_pan,[np.round(dims[0]/4),np.round(dims[1]/17)], method='neighbour')
        #np.clip(new_blue_band_pan,np.min(con_arr), np.max(con_arr), out=new_blue_band_pan)

        #con_arr=congrid.congrid(new_red_band_pan,[np.round(dims[0]/4),np.round(dims[1]/17)], method='neighbour')
        #np.clip(new_red_band_pan,np.min(con_arr), np.max(con_arr), out=new_red_band_pan)

        #con_arr=congrid.congrid(new_green_band_pan,[np.round(dims[0]/4),np.round(dims[1]/17)], method='neighbour')
        #np.clip(new_green_band_pan,np.min(con_arr), np.max(con_arr), out=new_green_band_pan)

	#con_arr=sci.ndimage.filters.median_filter(new_blue_band_pan,size=15)
	#np.clip(new_blue_band_pan,np.min(con_arr), np.max(con_arr), out=new_blue_band_pan)

	#con_arr=sci.ndimage.filters.median_filter(new_green_band_pan,size=15)
	#np.clip(new_green_band_pan,np.min(con_arr), np.max(con_arr), out=new_green_band_pan)

	#con_arr=sci.ndimage.filters.median_filter(new_red_band_pan,size=15)
	#np.clip(new_red_band_pan,np.min(con_arr), np.max(con_arr), out=new_red_band_pan)

	## Scale the data so the numbers do not exceed the limit of the data type
	#blue=np.round(new_blue_band_pan/np.max(new_blue_band_pan)*255.0)
        #green=np.round(new_green_band_pan/np.max(new_green_band_pan)*255.0)
        #red=np.round(new_red_band_pan/np.max(new_red_band_pan)*255.0)

	#print 'Min value after pan sharp and clip', np.min(new_blue_band_pan),np.min(new_green_band_pan),np.min(new_red_band_pan)	
	#print 'Max value after pan sharp and clip', np.max(new_blue_band_pan),np.max(new_green_band_pan),np.max(new_red_band_pan)
	
	#return blue,green,red
	return new_blue_band_pan, new_green_band_pan, new_red_band_pan

