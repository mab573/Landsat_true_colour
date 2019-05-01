#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# SCRIPT NAME		:: Atmosphric_correction.py 

# PURPOSE 		:: This script uses lookup files created from previous processes that will be used to calculate surface reflectance
#			:: values. This will be called from a perl script which will link the python modules together in the correct order.
#			:: The lookup files will be specifically created for each Hyperion file. It will also require a water vapour file
#			:: created from the Hyperion file as this will be used to select which lookup files to use. All other atmospheric parameters 
#			:: (at this stage) will be set to fixed average values. This code is based on (largely stolen - with permission) 
#			:: the work of Dr. Andrew Rodger. It's not a direct copy so there may be differences as it is based on a document
#			:: provided by Dr. Rodger and conversations and e-mails with Dr. Rodger.  

# SYNOPSIS 		:: (fill in later)    

# MODULES CALLED 	:: pyhdf, pyhdf.SD
#			:: numpy
#			:: sys
#			:: datetime, datetime
#			:: getopt, sys
#			:: os, glob
#			:: scipy, interpolate
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#**************** Load modules and such ***************

import matplotlib
import matplotlib.pyplot as plt
import pyhdf
from pyhdf.SD import *
import numpy as np
from scipy import *
from scipy import integrate
from datetime import datetime
import getopt, sys, math, Image, os, glob

#*******************************************************

def get_RTC(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][:,:]
    nc.close()
    return v

def check_outdir(in_file,out_file_base):
    # Check and create output directory if it doesn't already exist.
    in_dir, f_name = os.path.split(in_file)
    in_dir_splt = in_dir.split('/')
    a_time = in_dir_splt[len(in_dir_splt) - 1]
    a_day = in_dir_splt[len(in_dir_splt) - 2]
    a_month = in_dir_splt[len(in_dir_splt) - 3]
    a_year = in_dir_splt[len(in_dir_splt) - 4]
    outfile_year = os.path.join(out_file_base, a_year)
    outfile_month = os.path.join(outfile_year, a_month)
    outfile_day = os.path.join(outfile_month, a_day)
    outfile_time = os.path.join(outfile_day, a_time)


    if os.path.isdir(outfile_year):
        yr = ''
    else:
        os.system('mkdir '+outfile_year)
        yr = a_year
    if os.path.isdir(outfile_month):
        mt = ''
    else:
        os.system('mkdir ' +outfile_month)
        mt = a_month
    if os.path.isdir(outfile_day):
        dy=''
    else:
        os.system('mkdir ' + outfile_day)
        dy = a_day
    if os.path.isdir(outfile_time):
        hr = ''
    else:
        os.system('mkdir ' + outfile_time)
        hr = a_time

    #out_string = 'Created the following directories for year, month, day and hour: year /'+yr+ 'month /'+mt+' day /'+dy+' hour /'+hr

    out_path=outfile_time
    return out_path


def Landsat_ATCOR(LS_rad, rtc_path,band_num):


	########################################### Calculate the surface reflectance estimate per band ################################
	#__________________________________________ ################################################### ________________________________

	#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     


	print 'Performing atmospheric compensation.......'
	print ' '

        # Files interpolated from the LUT/files
 
        fileName = rtc_path + '/Band'+str(band_num)+'_RTC_interpolated_values.nc'
        Lp_0 = get_RTC(fileName, 'Lp_0')
        Eg_0 = get_RTC(fileName, 'Eg_0')
        T_up = get_RTC(fileName, 'T_up')
        S = get_RTC(fileName, 'S')

	# I think a scaling factor 10 will work as Hyperion has the same units and the DN had to be scaled by 40 to get to radiance and I had 400 as the scaling factor			
	A=(math.pi*((LS_rad/10.0)-(band_Lp_0)))/((band_Eg_0)*band_Gamma_up)
	rho=np.around(((A)/(1+(A*band_S)))*10000.0, decimals=0)

			
	print '                                   .......Completed atmospheric compensation'
	print ' '
	
	return rho

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



