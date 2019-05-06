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
#import pyhdf
#from pyhdf.SD import *
import numpy as np
from scipy import *
from scipy import integrate
from datetime import datetime
import getopt, sys, math, os, glob
import netCDF4

def get_RTC(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][:,:]
    nc.close()
    return v

#*******************************************************

def Landsat_ATCOR(LS_rad, band, RTC_dir):
   
 
    fileName=RTC_dir+'Band'+str(band)+'_L8_RTC_interpolated_values.nc'

    band_Lp_0 = get_RTC(fileName, 'Lp_0')
    band_Eg_0 = get_RTC(fileName, 'Eg_0')
    band_Gamma_up = get_RTC(fileName, 'T_up')
    band_S = get_RTC(fileName, 'S')
        

    ## This can be adapted later to use fully interpolated rasters for the atmospheric components
    ## Initially just use the centre coordinate and find the nearest table entry

    print('Performing atmospheric compensation.......\n')

    # I think a scaling factor 10 will work as Hyperion has the same units and the DN had to be scaled by 40 to get to radiance and I had 400 as the scaling factor			
    A=(math.pi*((LS_rad/10.0)-(band_Lp_0)))/((band_Eg_0)*band_Gamma_up)
    rho=np.around(((A)/(1+(A*band_S)))*10000.0, decimals=0)

			
    print ('                                   .......Completed atmospheric compensation\n')

	
    return rho

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#####################################################################################
#>
#>
#>
#>
################## Do other stuff #################################

'''
						
#########################################################################
#>
#>
#>
#>
######################### Write out data to a similar hdf file so that I can read it in Envi ##########################

# 	This creates a new hdf file to store the results of the smiley fixes.
#	It then writes out the new data cube and metadata into the hdf file. This essentially means copying
#	The metadata from the input file.
#	We could copy the header (.hdr) file here but there is no real point as this is essentially an interim product
#	The header file info wont be needed I think until we create the final product

# Get todays date to add to the file metadate for when stuff was processed
thingy=datetime.now()
out_date=thingy.strftime("%A, %d %B %Y %I:%M%p")

out_hdf_file = SD(HDF_filename_split[0]+'_AC_cibr_wv.'+HDF_filename_split[1],SDC.WRITE|SDC.CREATE) 

d1 = out_hdf_file.create(Dataset_name[0], SDC.INT16, (colm_size,band_size,line_size))

#print d1.dim(0),d1.dim(1),d1.dim(2)

d1[:,:,:]=rho[:,:,:]



# This writes out the attributes read from the input file
for attrName in file_attributes.keys():
            setattr(out_hdf_file, attrName, file_attributes[attrName])

# Can manually set attributes like this if you want
setattr(out_hdf_file, 'Post processing', 'Atmospheric compensation (reflectance product) using A/CHEAT, processed by Curtin RSSRG on '+out_date )


# NEED TO CHANGE THE DATA_ATTRIBUTES AS THESE ARE NOW IN REFLECTANCE RATHER THAN RADIANCE

# This writes out the attributes read from the input file
#for attrName in data_attributes.keys():
#            setattr(d1, attrName, data_attributes[attrName])

setattr(d1,'Number of Bands', band_size)
setattr(d1,'Number of Cross Track Pixels', line_size)
setattr(d1,'Number of Along Track Pixels', colm_size)
setattr(d1,'Data units', 'reflectance')
setattr(d1,'Scaling Factor', 'reflectance * 10000')



d1.endaccess()
out_hdf_file.end()
#######################################################################################################################

######################### Write data as an ENVI file as I recon this is how most             ##########################
######################### people would want the file anyway.                                 ##########################

# Copy a header file (created earlier) and change the bits and pieces that make it match the file just been processed


## envi_hdr_file contains a template which has parameter names that will match the required parameters from
## L1R_hdr_file. In most cases, when the parameter names match then copy the line from L1R_hdr_file

bil_file_name=HDF_filename_split[0]+'_At_Comp.bil'
bil_hdr_file_name=HDF_filename_split[0]+'_At_Comp.hdr'

# Open new file to write stuff to
new_file = open(bil_hdr_file_name, "w")

for xx in range(len(out_envi_hdr_file_lines)):
	
	new_file.writelines(out_envi_hdr_file_lines[xx])

new_file.close()


rho.tofile(bil_file_name)


print 'Finished processing Atmospheric correction at ', out_date

'''
