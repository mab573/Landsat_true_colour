#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# SCRIPT NAME		:: Landsat_basic_AC.py 

# PURPOSE 		:: This script will do basic atmospheric correction of Landsat8 band 4,3,2 and 8 (red, green, blue and panchromatic).
#				   

# SYNOPSIS 		:: It performs a partial atmospheric correction (accounting for Rayleigh scattering only) using observational geometric parameters 
#				:: from the Landsat metadata. The lookup table is generated from the radiative transfer program MODTRAN5 and the band equivalent
#				:: values are calculated using spectral response functions obtained from the internet.
#				:: These can now be pan sharpened using Brovey function in ENVI (or similar) 

# MODULES CALLED 	:: pyhdf
#					:: numpy
#					:: gdal
#					:: *Other modules may be called that aren't used but I tend to do this for everything.
#
# SUBROUTINES		:: Landsat8_create_tape5_files - Generates tape5 files to run MODTRAN5 based on Landsat8 metadata
#					:: Landsat8_make_final_lookups - Create a lookup file for use with the AC process 
#					:: Landsat8_atmospheric_correction - Performs the atmospheric correction/compensation and produces Bottom of Atmosphere (BOA) reflectance
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

from osgeo import gdal
from osgeo import osr
#import matplotlib as plt
from pylab import *
import numpy as np
import getopt, sys, math, os, glob, re, subprocess
import time
from pyproj import Proj
from PIL import Image, ImageEnhance, ImageFilter, ImageStat
import netCDF4

# My subroutines

# This is I think a non-permanant way to add module paths so that they can be called from
# a directory that is not on the python path.
#sys.path.append('/home/mbroomhall/src/Landsat_code/True_Colour/')

#import Landsat8_create_tape5_files
#import Landsat8_make_final_lookups
import Landsat8_atmospheric_correction
import Simple_Pan_Sharpen as pan
import Create_L8_interp_RTC as C_RTC
import solar_pos as sol_p
import ContEnh

## Written by Passang (3/5/19)
import calc_sat_solar

# This stuff is used to adjust brightness and contrast
# May help, may not.

#im_gamma_blue=1.0
#im_gamma_green=1.0
#im_gamma_red=1.0
#im_gamma=6.7
#im_gamma=10.
#ir_gamma=0.32

# Maximum ref value to set images to
# The higher this value is then the less cloud becomes saturated in the image.
# Very bright cloud at the edges can exceed this value by quite a bit.
# A more sensible value may be 1.1 or 1.2. If the clip_level is adjusted then the
# gamma values need to be adjusted as well.
clip_level_max=4.0
clip_level_min=0.0


# Adjust the contrast and sharpeness
c_enh=2.3
#This specifies the midpoint of the contrast enhancement
#I have pulled out and hacked the PIL contrast enhancement so I 
#can set this to be independant of the actual image. Previously this
#is set as a mean value for each image.
c_mid=191
s_enh=1.5


#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>********************

def usage ():
	print ('\n\n\n')
	print ('###################################################################################################################')
	print ('###################################################################################################################')
	print ('This program is designed to work with Landsat8 datasets downloaded from the USGS (http://earthexplorer.usgs.gov/).')
	print ('You may be able to get the files in exactly the same format elsewhere, I have never bothered.')
	print ('\n')
	print ('These files are in a tarred and zipped file and contain all 11 bands, a QA band and a metadata file')
	print ('\n')
	print ('This code will do the following with the following options:\n')
	print ('-z [.tar.gz file] -- Extract all bands, write these to a directory and produce an output image.\n')
	print ('-s [UL lat,UL lon,LR lat,LR lon] -- Enter like -21.5,115.3,-22.5,115.8 and the code will subset all bands')
	print ('					   and write them to a directory. It will produce an output image for the')
	print ('					   subset. Southern latitudes are -ve, Western longitudes are -ve.\n')
	print ('-c -- This will perform atmospheric correction on the bands used to produce the colour image')
	print ('	     (bands 4,3,2 and 8 or red, green, blue and panchromatic). If you do not have MODTRAN5 installed')
	print ('	     on your system then just ignore this option. If you do you will need to alter the code so that')
	print ('	     VIR_SRF, tape5_dir, data_dir and MOD_exe point to the correct places. This is the responsibility')
	print ('	     of the user.\n')
	print ('-p -- Set this and the images will be pan sharpened to 15 m resolution pixels, otherwise do at the native resolution\n') 
	print ('-b [a number] -- Enter a number to control how much the image is brightened. This is done for the output')
	print ('			image and will depend on the range of values in the image so a full swath and a subset')
	print ('			will likely scale differently. The impact of the scaling factor will also change depending')
	print ('			on which type of data is scaled (radiance or reflectance).')
	print ('			This scaling is similar to the gamma correction that many software packages can do. Entering')
	print ('			a small number (0.00001) will give an almost linear transform with a gradient of 1 which')
	print ('			does nothing. With reflectance data, numbers greater than about 0.01 and less than 0.5 will')
	print ('			in general make the image as bright as you will ever want. You may need to use numbers')
	print ('			greater than 1 for radiance data.\n')
	print ('-t -- If you image the full swath or the subset is at the edge of the swath, you will get a black boarder')
	print ('	     around your image. This switch will make the image boarder transparent. This should be used with some')
	print ('	     caution as this process makes black pixels transparent. Reflectance products especially will have')
	print ('	     black pixels where there are significant shadows i.e. behinds hills and mountains.\n')
	print ('-k -- This will keep all interim files and directories otherwise you will be left with the image and the .tar.gz')
	print ('	     file. The program will create a directory from the name of the .tar.gz file and a sub directory')
	print ('	     called simple_AC. The top level created directory will contain the bands extracted from the')
	print ('	     archive file as GEOTIF files as well as subsetted versions of these bands if you selected the -s option.')
	print ('	     The simple_AC folder will contain .TIF files for either bands 4,3,2,8 reflectance or radiance depending')
	print ('	     on the -c option and a .TIF file with bands 4,3,2 pan sharpened to 15m either as radiance or reflectance.')
	print ('	     You may want to keep these interim files to do something else with.\n')
	print ('-l -- This will output a logo to the bottom left hand corner of the image. It makes no attempt to resize')
	print ('	     the logo so if the subset is small it may cover most of the image. If you image the entire swath')
	print ('	     the logo will be very small and either copy over the black or tranparent swath boarder. If you')
	print ('	     want to use this option you need to alter the code so that logo_file points to your logo file location.\n')
	print ('-f [format extension] -- What format do you want your image in? This will output png, bmp, eps, jpeg and geotiff.')
	print ('				Enter as png, bmp, jpeg or tiff. Default will output bmp.')
	print ('-r -- This will save 30 m atmosperically corrected data as Rrs (pixel reflectances divided by pi). Default is not')
	print ('      to save Rrs') 
	print ('###################################################################################################################')
	print ('###################################################################################################################')
	print ('\n\n\n')  
	

try:
	opts, args = getopt.getopt(sys.argv[1:], "hz:cps:b:tklf:r", ["help", "zip_file","AC","pan","subset","bright","trans","keep","logo","format","Rrs"])

except getopt.GetoptError as err:
        # print help information and exit:
        print (str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-z", "--zip_file"):
            zip_file = a
        elif o in ("-c", "--AC"):
            AC = a
        elif o in ("-p", "--pan"):
            pan_sharpen = a
        elif o in ("-s", "--subset"):
            subset = a
        elif o in ("-b", "--bright"):
            bright = a
        elif o in ("-t", "--trans"):
            trans = a
        elif o in ("-k", "--keep"):
            keep = a
        elif o in ("-l", "--logo"):
            logo = a
        elif o in ("-f", "--format"):
            format = a
        elif o in ("-r", "--Rrs"):
            Rrs = a
        else:
            assert False, "unhandled option"
    # ...

#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>*******************<D-c>*

######$$$$$$$#########$$$$$$$$$##########$$$$$$$$$$

# These are either are required input files or directories.
# If you are not using MODTRAN or adding a logo to the images then
# you don't need any of them

#VIR_SRF='/home/mbroomhall/src/Landsat_code/True_Colour/template_files/Landsat8_OLI_relative_spectral_responses.csv'
#data_dir='/apps/MODTRAN5/DATA/'
#MOD_exe='/apps/MODTRAN5/Mod90_5.2.1.exe'
#logo_file='/home/mbroomhall/src/Landsat_code/True_Colour/template_files/curtin_logo.png'

######$$$$$$$#########$$$$$$$$$##########$$$$$$$$$$


## This is the maximum allowable value after reflectances have been pan sharpend
## All reflectances have been multiplied by 10000 and converted to integers (saves room)
## It is possible to get reflectance above 10000 in and around clouds which is why this
## is left here just in case someone really wants to look at clouds. Anything above 10000
## is set to 10000 so detail within clouds could be lost.

max_refl=10000


#%%%%%%%%^^^^^^^^^^%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%

# Functions for brightening the image. Feel free to find a better way to do this.

def log_bright(ary,bright):

        whereAreNaNs = np.isnan(ary);
        ary[whereAreNaNs] = 0;

        max_scale=np.max(np.log10(1+bright*ary))
        ary_scaled=np.log10(1+bright*ary)/max_scale*np.max(ary)

        return ary_scaled

def gamma(ary,bright):


    ary_scaled=((ary / 255.0) ** (1.0 /bright))*255.0

    return ary_scaled

def do_gamma(im, gamma):
    """Fast gamma correction with PIL's image.point() method"""
    invert_gamma = 1.0/gamma
    lut = [pow(x/255., invert_gamma) * 255 for x in range(256)]
    lut = lut*3 # need one set of data for each band for RGB
    im = im.point(lut)
    return im


#%%%%%%%%^^^^^^^^^^%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%

def open_geotiff(filename):
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    (X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
    srs_wkt = ds.GetProjection()
    Nx = ds.RasterXSize
    Ny = ds.RasterYSize
    #print srs_wkt
    #print X, deltaX, rotation, Y, rotation, deltaY
    #print Nx,Ny
    #print ' '
    ary = []
    ary = ds.GetRasterBand(1).ReadAsArray()
    return ary

def Write_NetCDF(filename,data):

    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(data[:, 0]))
    nc.createDimension("y_dim", len(data[0, :]))
    v = nc.createVariable("data", float, ("x_dim", "y_dim"))
    v[:, :] = data
    return

def gdal_sub(in_netcdf_file,chan_name,out_netcdf_file,out_dir):

    """!@brief subsets and warps himawari obs data to Australian radar composite
    @param in_nc_file float/int: a netcdf file
    @param out_path string: location that the output files are being written to
    @retval out_nc_file: float/int: a netcdf file (only smaller or bigger and different looking)
    """

    ## Australian RADAR netCDF projection
    radar_proj="+proj=aea +lat_1=-18.0 +lat_2=-36.0 +lon_0=132 +lat_0=0 +a=6378137 +b=6356752"
    radar_te="-2301000.000 -5099000.000 2599000.000 -999000.000"
    radar_tr="2000 2000"

    #radar_proj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    #radar_te="110.0 -45.0 155.0 -10.0"
    #radar_te=lon_1+" "+lat_1+" "+lon_0+" "+lat_0
    #radar_tr=res+" "+res

    temp_map=out_dir+'/map.tif'

    # warp the satellite to the radar grid
    os.system('gdalwarp NETCDF:'+in_netcdf_file+':'+chan_name+' '+temp_map+' -ot Float32 -co \"COMPRESS=DEFLATE\" -t_srs \"'+radar_proj+'\" -tr '+radar_tr+' -te '+radar_te+' -overwrite')


    # convert back to netcdf
    os.system('gdal_translate ' + temp_map + ' ' + out_netcdf_file + ' -of NETCDF -co \"COMPRESS=TRUE\"')

    return

## Used to assess elapsed time. Feel free to delete. Do the same at the other end as well
ts_0 = time.time()

'''
hoy=2
latitude_deg=-15.0
longitude_deg=115.0
year=2019

Sol_el=[]
Sol_az=[]

for jj in range (0,40,5):

    for ii in range(0,366):
        print ii,jj
        result=sol_p.calc_sun_position(latitude_deg-float(jj), longitude_deg, year, hoy+(24*ii))
        Sol_el.append(result[0])
        Sol_az.append(result[1])

plt.plot(Sol_az,'b.')
plt.show()

sys.exit()
'''

# Untar and unzip the Landsat8 files in a directory named after the archive file
# The Landsat file names contain info on target path and row, year and day of year of capture which is a unique identifier

# Split file name from the path

in_path, in_file = os.path.split(zip_file)
print('in_path', in_path)
print('in_file', in_file)
get_basename=in_file.split('.')
basename=get_basename[0]
out_dir_top=in_path+'/'+basename
print('out_dir_top', out_dir_top)
os.system("mkdir "+out_dir_top)
# unpack the archive file to the new directory

LS_meta_file=out_dir_top+'/'+basename+'_MTL.txt'
LS_B1=out_dir_top+'/'+basename+'_B1.TIF'
LS_B2=out_dir_top+'/'+basename+'_B2.TIF'
LS_B3=out_dir_top+'/'+basename+'_B3.TIF'
LS_B4=out_dir_top+'/'+basename+'_B4.TIF'
LS_B5=out_dir_top+'/'+basename+'_B5.TIF'
LS_B6=out_dir_top+'/'+basename+'_B6.TIF'
LS_B7=out_dir_top+'/'+basename+'_B7.TIF'
LS_B8=out_dir_top+'/'+basename+'_B8.TIF'
LS_B9=out_dir_top+'/'+basename+'_B9.TIF'
LS_B10=out_dir_top+'/'+basename+'_B10.TIF'
LS_B11=out_dir_top+'/'+basename+'_B11.TIF'

LS_B1_sub=out_dir_top+'/'+basename+'_sub_B1.TIF'
LS_B2_sub=out_dir_top+'/'+basename+'_sub_B2.TIF'
LS_B3_sub=out_dir_top+'/'+basename+'_sub_B3.TIF'
LS_B4_sub=out_dir_top+'/'+basename+'_sub_B4.TIF'
LS_B5_sub=out_dir_top+'/'+basename+'_sub_B5.TIF'
LS_B6_sub=out_dir_top+'/'+basename+'_sub_B6.TIF'
LS_B7_sub=out_dir_top+'/'+basename+'_sub_B7.TIF'
LS_B8_sub=out_dir_top+'/'+basename+'_sub_B8.TIF'
LS_B9_sub=out_dir_top+'/'+basename+'_sub_B9.TIF'
LS_B10_sub=out_dir_top+'/'+basename+'_sub_B10.TIF'
LS_B11_sub=out_dir_top+'/'+basename+'_sub_B11.TIF'

print ('Extracting landsat raster files from the archive file - if required\n\n\n')
# Untar and unzip L8 file, create an output directory
out_dir=in_path+'/'+basename+'/simple_AC/'
if not os.path.isfile(LS_B1):
    os.system("tar -zxvf "+zip_file+" -C "+out_dir_top)
    os.system("mkdir "+out_dir)

try:
  bright
except NameError:
  print ('\n\n\n Image brightening scaling factor has not been set, automatically set to 0.1\n\n\n')
  bright=np.float(0.1)
  bright=np.float(1.0)

else:
  bright=np.float(bright)

## These can be altered individually if needed
im_gamma_blue=bright
im_gamma_green=bright
im_gamma_red=bright
  

#Lists of metadata tags that we want to find values for

rad_scale_tags=['RADIANCE_MULT_BAND_1','RADIANCE_MULT_BAND_2','RADIANCE_MULT_BAND_3','RADIANCE_MULT_BAND_4','RADIANCE_MULT_BAND_5','RADIANCE_MULT_BAND_6',
'RADIANCE_MULT_BAND_7','RADIANCE_MULT_BAND_8','RADIANCE_MULT_BAND_9','RADIANCE_MULT_BAND_10','RADIANCE_MULT_BAND_11']

rad_offset_tags=['RADIANCE_ADD_BAND_1','RADIANCE_ADD_BAND_2','RADIANCE_ADD_BAND_3','RADIANCE_ADD_BAND_4','RADIANCE_ADD_BAND_5',
'RADIANCE_ADD_BAND_6','RADIANCE_ADD_BAND_7','RADIANCE_ADD_BAND_8','RADIANCE_ADD_BAND_9','RADIANCE_ADD_BAND_10','RADIANCE_ADD_BAND_11']

ref_scale_tags=['REFLECTANCE_MULT_BAND_1','REFLECTANCE_MULT_BAND_2','REFLECTANCE_MULT_BAND_3','REFLECTANCE_MULT_BAND_4','REFLECTANCE_MULT_BAND_5',
'REFLECTANCE_MULT_BAND_6','REFLECTANCE_MULT_BAND_7','REFLECTANCE_MULT_BAND_8','REFLECTANCE_MULT_BAND_9','REFLECTANCE_MULT_BAND_10','REFLECTANCE_MULT_BAND_11']

ref_offset_tags=['REFLECTANCE_ADD_BAND_1','REFLECTANCE_ADD_BAND_2','REFLECTANCE_ADD_BAND_3','REFLECTANCE_ADD_BAND_4','REFLECTANCE_ADD_BAND_5',
'REFLECTANCE_ADD_BAND_6','REFLECTANCE_ADD_BAND_7','REFLECTANCE_ADD_BAND_8','REFLECTANCE_ADD_BAND_9','REFLECTANCE_ADD_BAND_10','REFLECTANCE_ADD_BAND_11']

sun_az_tags='SUN_AZIMUTH'
sun_el_tags='SUN_ELEVATION'
sun_dist_tags='EARTH_SUN_DISTANCE'
scene_cent_time='SCENE_CENTER_TIME'
date='DATE_ACQUIRED'
file_name='LANDSAT_SCENE_ID'
elips='ELLIPSOID'
zone='UTM_ZONE'
map_proj='MAP_PROJECTION'

UL_lat='CORNER_UL_LAT_PRODUCT'
UL_lon='CORNER_UL_LON_PRODUCT'
UR_lat='CORNER_UR_LAT_PRODUCT'
UR_lon='CORNER_UR_LON_PRODUCT'
LL_lat='CORNER_LL_LAT_PRODUCT'
LL_lon='CORNER_LL_LON_PRODUCT'
LR_lat='CORNER_LR_LAT_PRODUCT'
LR_lon='CORNER_LR_LON_PRODUCT'

rad_scale=[]
rad_offset=[]
ref_scale=[]
ref_offset=[]


meta_file_lines=[]
with open(LS_meta_file,'r') as f:

	for line in f:
		
		meta_file_lines.append(line)

f.closed


print('meta', meta_file_lines[0])
print('sun', sun_az_tags)

new_list = [j for j in meta_file_lines if re.search(sun_az_tags, j)]
for item in new_list:
    thingy=item.split('=')
    sun_az=float(thingy[1])
    #print 'Solar azimuth ',sun_az
    
new_list = [j for j in meta_file_lines if re.search(sun_el_tags, j)]
for item in new_list:
    thingy=item.split('=')
    sun_el=float(thingy[1])
    #print 'Solar elevation ', sun_el

new_list = [j for j in meta_file_lines if re.search(sun_dist_tags, j)]
for item in new_list:
    thingy=item.split('=')
    sun_dist=float(thingy[1])
    #print 'Earth-Sun distance ',sun_dist

new_list = [j for j in meta_file_lines if re.search(scene_cent_time, j)]
for item in new_list:
    thingy=item.split('=')
    scene_time=thingy[1]
    scene_time=scene_time[2:-2]
    #print 'tile time ', scene_time

new_list = [j for j in meta_file_lines if re.search(date, j)]
for item in new_list:
    thingy=item.split('=')
    scene_date=thingy[1]
    scene_date=scene_date[2:-2]
    #print 'tile date ', scene_time

new_list = [j for j in meta_file_lines if re.search(file_name, j)]
for item in new_list:
    thingy=item.split('=')
    f_name=thingy[1]

new_list = [j for j in meta_file_lines if re.search(UL_lat, j)]
for item in new_list:
    thingy=item.split('=')
    UL_lat_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(UL_lon, j)]
for item in new_list:
    thingy=item.split('=')
    UL_lon_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(UR_lat, j)]
for item in new_list:
    thingy=item.split('=')
    UR_lat_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(UR_lon, j)]
for item in new_list:
    thingy=item.split('=')
    UR_lon_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(LL_lat, j)]
for item in new_list:
    thingy=item.split('=')
    LL_lat_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(LL_lon, j)]
for item in new_list:
    thingy=item.split('=')
    LL_lon_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(LR_lat, j)]
for item in new_list:
    thingy=item.split('=')
    LR_lat_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(LR_lon, j)]
for item in new_list:
    thingy=item.split('=')
    LR_lon_co=float(thingy[1])

new_list = [j for j in meta_file_lines if re.search(map_proj, j)]
for item in new_list:
    thingy=item.split('=')
    map_projection=thingy[1]
    
new_list = [j for j in meta_file_lines if re.search(elips, j)]
for item in new_list:
    thingy=item.split('=')
    elipsoid=thingy[1]
    
new_list = [j for j in meta_file_lines if re.search(zone, j)]
for item in new_list:
    thingy=item.split('=')
    utm_zone=int(thingy[1])

# Estimated central position of the Landsat swath (not particularly accurate but should be within a minute of actual centre - probably)

centre_lat=((LR_lat_co - UL_lat_co)/2)+UL_lat_co
centre_lon=(((UR_lon_co - LL_lon_co)/2)+LL_lon_co)*-1.0		# This is done as it will be fed into MODTRAN where Eastern longitudes are treated as negative.

scene_centre=[]

target_lat="%.3f" % centre_lat
target_lon="%.3f" % centre_lon

scene_centre.append(target_lat)
scene_centre.append(target_lon)

# Extract scene date as the day of the year from the filename

#day_of_year=f_name[13:16]
day_of_year=basename[13:16]		#easier to pull it out of the file name than getting it from the metafile. This is ok unless they change the naming convention.


ST_split=scene_time.split(':')
sec_split=ST_split[2].split('.')

scene_time_float=float(ST_split[0])+float(ST_split[1])/60.0+float(sec_split[0])/3600.0
scene_time="%.3f" % scene_time_float

#print 'UL Lat ',UL_lat_co
#print 'LR_lat ',LR_lat_co
#print 'LR_lon ',LR_lon_co
#print 'Tile centre coordinates '+str(centre_lat)+' '+str(centre_lon)
#print 'Tile time ',scene_time


#### If the subset switch is used then get everything in the right format so that it can be 
#### subsetted properly.
#### These represent the upper left and lower right coordinates of the grid you want to display.
#### This doesn't do any cleaver error checking to make sure your coordinates are sensable so
#### you will have to do this yourself. The following code converts geographic lat/lon to UTM

try:
  subset
except NameError:
  print ('\n\n\n Subsetting has not been selected, imaging the entire tile with brighness factor ',bright,'\n\n\n')
  do_sub=0
else:
  do_sub=1 	
  subset_vals=subset.split(',')
  ulat = float(subset_vals[0])
  ulon = float(subset_vals[1])
  llat = float(subset_vals[2])
  llon = float(subset_vals[3])
  
  ## strip off the double quotes
  map_projection=map_projection[2:-2]
  elipsoid=elipsoid[2:-2]

  # This uses a python module called pyproj. Map projection must be in lowercase - for some reason
  p = Proj(proj=map_projection.lower(),zone=utm_zone,ellps=elipsoid)
  ul_E,ul_N = p(ulon, ulat)
  lr_E,lr_N = p(llon, llat)


for ii in range(0,len(rad_scale_tags)):

    new_list = [j for j in meta_file_lines if re.search(rad_scale_tags[ii], j)]
    for item in new_list:

        thingy=item.split('=')
        stuff=thingy[0].split('_')

        if stuff[3].rstrip() == '1':

            rad_scale.append(float(thingy[1]))

        elif stuff[3].rstrip() != '10' and stuff[3].rstrip() != '11':

            rad_scale.append(float(thingy[1]))	

    new_list = [j for j in meta_file_lines if re.search(rad_offset_tags[ii], j)]
    for item in new_list:
    		
        thingy=item.split('=')
        stuff=thingy[0].split('_')

        if stuff[3].rstrip() == '1':

            rad_offset.append(float(thingy[1]))

        elif stuff[3].rstrip() != '10' and stuff[3].rstrip() != '11':

            thingy=item.split('=')
            rad_offset.append(float(thingy[1]))

    new_list = [j for j in meta_file_lines if re.search(ref_scale_tags[ii], j)]
    for item in new_list:

        thingy=item.split('=')
        ref_scale.append(float(thingy[1]))

    new_list = [j for j in meta_file_lines if re.search(ref_offset_tags[ii], j)]
    for item in new_list:
    			
        thingy=item.split('=')
        ref_offset.append(float(thingy[1]))

# Write tape5 files based on the information from the metadata file. This will be scene centre position, day of year of the scene, scene centre time of acquisition.
# This will used fixed inputs for water vapour, aerosol, atmospheric pressure and so on.
# This should be reasonably equivalent to CREFL proessing for MODIS.

try:
  AC
except NameError:

  ### No need to call MODTRAN, so dont.
  print ('\n\n\n Atmospheric correction has not been selected. Processing radiance data instead \n\n\n')

  do_AC=0

else:
  
  print ('\n\n\n Atmospherically correcting radiance data to surface reflectance. \n\n\n')

  # None of this is required anymore as I have a LUT
  #tape5_dir=out_dir

  #Landsat8_create_tape5_files.write_tape5(day_of_year,scene_centre, scene_time, tape5_dir, data_dir)

  #subprocess.call([MOD_exe])

  #Landsat8_make_final_lookups.Make_final_lookups(tape5_dir)

  do_AC=1

# Read in the spectral response finctions for Landsat 8

#VIS_SRF_table=np.loadtxt(VIR_SRF, dtype=float, delimiter=',', skiprows = 2)

# Band1 0,1
# Band2 2,3
# Band3 4,5
# Band4 6,7
# etc but the order is band 1 - 9 but the spectral ranges are not consecutive

try:
  pan_sharpen
except NameError:

  do_pan=0
  print ('\n\n\n Not pan sharpening the images. \n\n\n')
else:

  ### No need to call MODTRAN, so dont.
  print ('\n\n\n Pan sharpening the image. \n\n\n')

  do_pan=1


try:
  Rrs
except NameError:

  ### No need to call MODTRAN, so dont.
  print ('\n\n\n Not saving the first 4 bands as Rrs. \n\n\n')

  do_Rrs=0

else:

  do_Rrs=1
  print ('\n\n\n Saving the first 4 bands as Rrs. \n\n\n')


# Temporary hard coded path
path='/g/data/if87/ARD_interoperability/ga-packaged_collection/2019-02-20/LC80930852019051LGN00/SUPPLEMENTARY/'
out_file_base=path+'RTC/'
#VZA=open_geotiff(path+'LC80930852019051LGN00_SATELLITE_VIEW.TIF')
#SZA=open_geotiff(path+'LC80930852019051LGN00_SOLAR_ZENITH.TIF')
#VA=open_geotiff(path+'LC80930852019051LGN00_SATELLITE_AZIMUTH.TIF')
#SA=open_geotiff(path+'LC80930852019051LGN00_SOLAR_AZIMUTH.TIF')

VZA=path+'LC80930852019051LGN00_SATELLITE_VIEW.TIF'
SZA=path+'LC80930852019051LGN00_SOLAR_ZENITH.TIF'
VA=path+'LC80930852019051LGN00_SATELLITE_AZIMUTH.TIF'
SA=path+'LC80930852019051LGN00_SOLAR_AZIMUTH.TIF'

'''
## Temporary hard coded path
path='/short/er8/mab573/Landsat/'
out_file_base=path+'RTC/'
VZA=open_geotiff(path+'SATELLITE-VIEW.tif')
SZA=open_geotiff(path+'SOLAR-ZENITH.tif')
VA=open_geotiff(path+'SATELLITE-AZIMUTH.tif')
SA=open_geotiff(path+'SOLAR-AZIMUTH.tif')
'''

path='/short/er8/mab573/Landsat/'
VZA_sub=path+'SATELLITE-VIEW_sub.TIF'
SZA_sub=path+'SOLAR-ZENITH_sub.TIF'
VA_sub=path+'SATELLITE-AZIMUTH_sub.TIF'
SA_sub=path+'SOLAR-AZIMUTH_sub.TIF'

VZA_sub_B8=path+'SATELLITE-VIEW_sub_B8.TIF'
SZA_sub_B8=path+'SOLAR-ZENITH_sub_B8.TIF'
VA_sub_B8=path+'SATELLITE-AZIMUTH_sub_B8.TIF'
SA_sub_B8=path+'SOLAR-AZIMUTH_sub_B8.TIF'


#### In previous versions, subsetting is done near the end of the processing.
#### This causes problems. Subset all bands first.

if do_sub == 1:

    if do_pan == 1:
        ## Subset at 15 m resolution
        Bands = ['B2', 'B3', 'B4', 'B8']
        VZA,SZA,VAA,SAA,B2,B3,B4,B8=calc_sat_solar.main(zip_file, path, subset, Bands)

    else:
        ## Subset at native resolution and we dont need the pan band
        Bands = ['B2', 'B3', 'B4']
        VZA,SZA,VAA,SAA,B2,B3,B4=calc_sat_solar.main(zip_file, path, subset, Bands)

else:

    LS_B1_sub=LS_B1
    LS_B2_sub=LS_B2
    LS_B3_sub=LS_B3
    LS_B4_sub=LS_B4
    LS_B5_sub=LS_B5
    LS_B6_sub=LS_B6
    LS_B7_sub=LS_B7
    LS_B8_sub=LS_B8
    LS_B9_sub=LS_B9
    LS_B10_sub=LS_B10
    LS_B11_sub=LS_B11

    VZA_sub=VZA
    SZA_sub=SZA
    VA_sub=VA
    SA_sub=SA

    print('we are missing calc_sat_solar here')

print('ba-bye')
sys.exit()
print("you'll never see this")

if do_sub == 1:

    if do_pan ==1:
        print ('Subsetting the scene with upper left geographic coordinates',ulat,ulon,' UTM ',ul_E,ul_N)
        print ('and low right geographic coordinates ',llat,llon,' UTM ',lr_E,lr_N,' with brightening factor ',bright,'\n\n\n')

        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B1,LS_B1_sub))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B2,LS_B2_sub))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B3,LS_B3_sub))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B4,LS_B4_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B5,LS_B5_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B6,LS_B6_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B7,LS_B7_sub))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B8,LS_B8_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B9,LS_B9_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B10,LS_B10_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B11,LS_B11_sub))

        ## ADD subsetting of the VZA, SZA, VA and SA and use this to generate RTC raster of the correct size
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, VZA,VZA_sub_B8))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, SZA,SZA_sub_B8))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, VA,VA_sub_B8))
        os.system('gdal_translate  -tr 15 15 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, SA,SA_sub_B8))

    else:
        print ('Subsetting the scene with upper left geographic coordinates',ulat,ulon,' UTM ',ul_E,ul_N)
        print ('and low right geographic coordinates ',llat,llon,' UTM ',lr_E,lr_N,' with brightening factor ',bright,'\n\n\n')
 
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B1,LS_B1_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B2,LS_B2_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B3,LS_B3_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B4,LS_B4_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B5,LS_B5_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B6,LS_B6_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B7,LS_B7_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B8,LS_B8_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B9,LS_B9_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B10,LS_B10_sub))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, LS_B11,LS_B11_sub))
 
        ## ADD subsetting of the VZA, SZA, VA and SA and use this to generate RTC raster of the correct size
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, VZA,VZA_sub_B8))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, SZA,SZA_sub_B8))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, VA,VA_sub_B8))
        os.system('gdal_translate  -tr 30 30 -projwin %s %s %s %s %s -b 1 %s' %(ul_E, ul_N, lr_E, lr_N, SA,SA_sub_B8))


else:
    LS_B1_sub=LS_B1
    LS_B2_sub=LS_B2
    LS_B3_sub=LS_B3
    LS_B4_sub=LS_B4
    LS_B5_sub=LS_B5
    LS_B6_sub=LS_B6
    LS_B7_sub=LS_B7
    LS_B8_sub=LS_B8
    LS_B9_sub=LS_B9
    LS_B10_sub=LS_B10
    LS_B11_sub=LS_B11

    VZA_sub=VZA
    SZA_sub=SZA
    VA_sub=VA
    SA_sub=SA 
    ## If you want to do this full scene, you would still need to generate VZA etc at 15 m for band 8 RTC generation

############ Generate RTC interpolated Raster files ############

#V_Z_A=open_geotiff(VZA_sub)
#S_Z_A=open_geotiff(SZA_sub)
#V_A=open_geotiff(VA_sub)
#S_A=open_geotiff(SA_sub)

C_RTC.generate_RTC_rasters(VZA_sub,SZA_sub,VA_sub,SA_sub,VZA_sub_B8,SZA_sub_B8,VA_sub_B8,SA_sub_B8)

#PROCESS THE RED BAND (BAND4). OPEN, ATMOSPERIC COMPENSATION AND WRITE OUT RESULTANT BAND
ds = gdal.Open(LS_B4_sub, gdal.GA_ReadOnly)
(X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
srs_wkt = ds.GetProjection()
Nx = ds.RasterXSize
Ny = ds.RasterYSize
ary = []
ary = ds.GetRasterBand(1).ReadAsArray()
rad_arr=(ary[:].astype(float)*rad_scale[3])+rad_offset[3]


# Call the Landsat AC subroutine and return the resultant reflectance array
# Will probably have to un hard code this path
RTC_dir='/short/er8/mab573/Landsat/RTC/'

if do_AC==1:

	print ('Atmosperically compensating for band 4 ....\n\n')
	rho=Landsat8_atmospheric_correction.Landsat_ATCOR(rad_arr, 4, RTC_dir)

	# Remove any negative values
	rho_out_b4 = np.empty((Ny,Nx),dtype=int) # I don't know why but these need to be swapped
	np.clip(rho,0, max_refl, out=rho_out_b4)

	if do_Rrs==1:

                Rrs_out_b4=rho_out_b4/np.pi

else:
	
	rho_out_b4=rad_arr

# Write reflectance array to a Geotiff file

# Set file vars

output_file_B4 = out_dir+basename+'_B4_SAC.TIF'

# Create gtif
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create(output_file_B4, Nx, Ny, 1 , gdal.GDT_UInt16)

# top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

# set the reference info 
srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")
dst_ds.SetProjection( srs_wkt )

# write the band
dst_ds.GetRasterBand(1).WriteArray(rho_out_b4)



#Open and save image of green band
ds = gdal.Open(LS_B3_sub, gdal.GA_ReadOnly)
(X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
srs_wkt = ds.GetProjection()
Nx = ds.RasterXSize
Ny = ds.RasterYSize
ary = []
ary = ds.GetRasterBand(1).ReadAsArray()
rad_arr=(ary[:].astype(float)*rad_scale[2])+rad_offset[2]


if do_AC==1:

	print ('Atmosperically compensating for band 3 ....\n\n')
	rho=Landsat8_atmospheric_correction.Landsat_ATCOR(rad_arr, 3, RTC_dir)
	rho_out_b3 = np.empty((Ny,Nx),dtype=int)
	np.clip(rho,0, max_refl, out=rho_out_b3)

	if do_Rrs==1:
		
		Rrs_out_b3=rho_out_b3/np.pi

else:

	rho_out_b3=rad_arr



# Set file vars
output_file_B3 = out_dir+basename+'_B3_SAC.TIF'

# Create gtif
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create(output_file_B3, Nx, Ny, 1 , gdal.GDT_UInt16)

# top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

# set the reference info 
srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")
dst_ds.SetProjection( srs_wkt )

# write the band
dst_ds.GetRasterBand(1).WriteArray(rho_out_b3)


#Open and save image of blue band
ds = gdal.Open(LS_B2_sub, gdal.GA_ReadOnly)
(X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
srs_wkt = ds.GetProjection()
Nx = ds.RasterXSize
Ny = ds.RasterYSize
ary = []
ary = ds.GetRasterBand(1).ReadAsArray()
rad_arr=(ary[:].astype(float)*rad_scale[1])+rad_offset[1]

#Band 2

if do_AC==1:

	print ('Atmospherically compensating for band 2....\n\n')
	rho=Landsat8_atmospheric_correction.Landsat_ATCOR(rad_arr, 2, RTC_dir)
	rho_out_b2 = np.empty((Ny,Nx),dtype=int)
	np.clip(rho,0, max_refl, out=rho_out_b2)

	if do_Rrs==1:

		Rrs_out_b2=rho_out_b2/np.pi

else:

	rho_out_b2=rad_arr



# Set file vars
output_file_B2 = out_dir+basename+'_B2_SAC.TIF'

# Create gtif
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create(output_file_B2, Nx, Ny, 1 , gdal.GDT_UInt16)

# top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

# set the reference info 
srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")
dst_ds.SetProjection( srs_wkt )

# write the band
dst_ds.GetRasterBand(1).WriteArray(rho_out_b2)


######### Process band 1 as well. This is not used to generate the true colour image.
'''
#PROCESS THE RED BAND (BAND1). OPEN, ATMOSPERIC COMPENSATION AND WRITE OUT RESULTANT BAND
ds = gdal.Open(LS_B1_sub, gdal.GA_ReadOnly)
(X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
srs_wkt = ds.GetProjection()
Nx = ds.RasterXSize
Ny = ds.RasterYSize
ary = []
ary = ds.GetRasterBand(1).ReadAsArray()
rad_arr=(ary[:].astype(float)*rad_scale[0])+rad_offset[0]


# Call the Landsat AC subroutine and return the resultant reflectance array
# I haven't generated RTC for 1. I don't know why I used to do it.

if do_AC==1:

        print 'Atmosperically compensating for band 1 ....\n\n'
        rho=Landsat8_atmospheric_correction.Landsat_ATCOR(rad_arr,filt_res,tape5_dir+'MOD_ATM-cm_WV_0.10_0.50.7sc_lookup.txt', tape5_dir+'lookup_wavelengths.txt')

        # Remove any negative values
        rho_out_b1 = np.empty((Ny,Nx),dtype=int) # I don't know why but these need to be swapped
        np.clip(rho,0, max_refl, out=rho_out_b1)

	if do_Rrs==1:

                Rrs_out_b1=rho_out_b1/np.pi

else:

        rho_out_b1=rad_arr

# Write reflectance array to a Geotiff file

# Set file vars
output_file_B1 = out_dir+basename+'_B1_SAC.TIF'

# Create gtif
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create(output_file_B1, Nx, Ny, 1 , gdal.GDT_UInt16)

# top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

# set the reference info
srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")
dst_ds.SetProjection( srs_wkt )

# write the band
dst_ds.GetRasterBand(1).WriteArray(rho_out_b1)

#########
'''


########################################### Save the three colour bands to a single tiff file for other processing purposes
# Set file vars

output_file_Lo_res = out_dir+basename+'_RGB_SAC.TIF'

# Create gtif
driver = gdal.GetDriverByName("GTiff")
dst_ds = driver.Create(output_file_Lo_res, Nx, Ny, 3 , gdal.GDT_UInt16)

# top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

# set the reference info 
srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")
dst_ds.SetProjection( srs_wkt )

# write the band

g = gdal.Open(output_file_B2) # blue band
band1 = g.ReadAsArray()
dst_ds.GetRasterBand(1).WriteArray(band1)

g = gdal.Open(output_file_B3) # green band
band1 = g.ReadAsArray()
dst_ds.GetRasterBand(2).WriteArray(band1)

g = gdal.Open(output_file_B4) # red band
band1 = g.ReadAsArray()
dst_ds.GetRasterBand(3).WriteArray(band1)


##############################################


########################################### Save the first 4 bands as Rrs to a geotiff file for other processing purposes

if do_AC==1:

    if do_Rrs==1:

        Rrs_file_Lo_res = out_dir+basename+'_Rrs_SAC.TIF'

        # Create gtif
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(Rrs_file_Lo_res, Nx, Ny, 4 , gdal.GDT_UInt16)

        # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

        # set the reference info 
        srs = osr.SpatialReference()
        srs.SetWellKnownGeogCS("WGS84")
        dst_ds.SetProjection( srs_wkt )

        # write the band

        #g = gdal.Open(output_file_B1) # Coastal band
        #band1 = g.ReadAsArray()
        dst_ds.GetRasterBand(1).WriteArray(Rrs_out_b1)

        #g = gdal.Open(output_file_B2) # blue band
        #band1 = g.ReadAsArray()
        dst_ds.GetRasterBand(2).WriteArray(Rrs_out_b2)

        #g = gdal.Open(output_file_B3) # green band
        #band1 = g.ReadAsArray()
        dst_ds.GetRasterBand(3).WriteArray(Rrs_out_b3)

        #g = gdal.Open(output_file_B4) # red band
        #band1 = g.ReadAsArray()
        dst_ds.GetRasterBand(4).WriteArray(Rrs_out_b4)


        print ('\n\n\n Saving Rrs files to multiband geotiff file \n\n\n')




#Open and save image of the panchromatic band
ds = gdal.Open(LS_B8_sub, gdal.GA_ReadOnly)
(X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
srs_wkt = ds.GetProjection()
Nx = ds.RasterXSize
Ny = ds.RasterYSize
ary = []
ary = ds.GetRasterBand(1).ReadAsArray()
rad_arr=(ary[:].astype(float)*rad_scale[7])+rad_offset[7]


#Band 8

if do_AC==1:

    if do_pan==1:
        print ('Atmospherically compensating for band 8....\n\n')

        ## This is a hack but the band8 rad array is always a different shape to the
        ## RTC arrays. Need to pad in one direction and take in the other

        #rad_arr=np.take(rad_arr)
        #rad_arr=np.pad(rad_arr,((0,0),(0,1)),'edge')
        #rad_arr=rad_arr[0:-1,:]

        #print rad_arr.shape

        rho=Landsat8_atmospheric_correction.Landsat_ATCOR(rad_arr, 8, RTC_dir)
        rho_out_b8 = np.empty((rho.shape),dtype=int)
        np.clip(rho,0, max_refl, out=rho_out_b8)

else:

    rho_out_b8=rad_arr

## Set file vars
#output_file_B8 = out_dir+basename+'_B8_SAC.TIF'

## Create gtif
#driver = gdal.GetDriverByName("GTiff")
#dst_ds = driver.Create(output_file_B8, Nx, Ny, 1 , gdal.GDT_UInt16)

## top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
#dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

## set the reference info 
#srs = osr.SpatialReference()
#srs.SetWellKnownGeogCS("WGS84")
#dst_ds.SetProjection( srs_wkt )

## write the band
#dst_ds.GetRasterBand(1).WriteArray(rho_out_b8)


############# Pan sharpen using Brovey (or something like it) transform

if do_pan ==1:
    print ('Performing pan sharpening using simple Brovey transform\n\n')

    pan_sharp_blue, pan_sharp_green, pan_sharp_red=pan.Simple_Pan_Sharpen(rho_out_b2, rho_out_b3, rho_out_b4, rho_out_b8)

else:
    pan_sharp_blue=rho_out_b2
    pan_sharp_green=rho_out_b3
    pan_sharp_red=rho_out_b4

## This saves the pan sharpened rasters to a file. The GEOTIFF format from the PAN band (band 8) is used to create
## the file and the 3 raster bands written into the new file. 

## Set file vars
#output_file_RGB = out_dir+basename+'_RGB_SAC_PAN.TIF'

## Create gtif
#driver = gdal.GetDriverByName("GTiff")
#dst_ds = driver.Create(output_file_RGB, Nx, Ny, 3 , gdal.GDT_UInt16)

## top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
#dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )

## set the reference info 
#srs = osr.SpatialReference()
#srs.SetWellKnownGeogCS("WGS84")
#dst_ds.SetProjection( srs_wkt )

# write the band

#dst_ds.GetRasterBand(1).WriteArray(pan_sharp_red)
#dst_ds.GetRasterBand(2).WriteArray(pan_sharp_green)
#dst_ds.GetRasterBand(3).WriteArray(pan_sharp_blue)

##########################################


if do_sub == 1:
	
	out_image_name=in_path+'/'+basename+'_RGB_AC_PAN_CC_SUB.bmp'
	out_image_name_tif=in_path+'/'+basename+'_RGB_AC_PAN_CC_SUB.tif'
	out_image_name_png=in_path+'/'+basename+'_RGB_AC_PAN_CC_SUB.png'
	out_image_name_jpeg=in_path+'/'+basename+'_RGB_AC_PAN_CC_SUB.jpeg'
else:
		
	out_image_name=in_path+'/'+basename+'_RGB_AC_PAN_CC.bmp'
	out_image_name_tif=in_path+'/'+basename+'_RGB_AC_PAN_CC.tif'
	out_image_name_png=in_path+'/'+basename+'_RGB_AC_PAN_CC_SUB.png'
	out_image_name_jpeg=in_path+'/'+basename+'_RGB_AC_PAN_CC_SUB.jpeg'


## Create gtif
#driver = gdal.GetDriverByName("GTiff")
#dst_ds = driver.Create(out_image_name_tif, Nx, Ny, 3 , gdal.GDT_Byte)
## top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
#dst_ds.SetGeoTransform( [ X, deltaX, rotation, Y, rotation, deltaY ] )
## set the reference info
#srs = osr.SpatialReference()
#srs.SetWellKnownGeogCS("WGS84")
#dst_ds.SetProjection( srs_wkt )

# The radiance values are quite a bit smaller than the reflectance values (once they are converted to integers)
# The scaling does not work particularly well when most of the nubers are under 100.0
# As we only want to produce a picture at this stage, make the numbers bigger by multiplication.

if do_AC==1:

    print ('Now Scaling red band....\n\n\n')
    #ary_scaled=log_bright(pan_sharp_red)
    # Another , probably better, way to brighten images
    red = gamma(np.clip(pan_sharp_red.astype(float)/np.max(pan_sharp_red.astype(float))*255.0, 0,255), im_gamma_red)
        
    #red=ary_scaled/np.max(ary_scaled)*255.0
    jr=Image.fromarray(red.astype(np.uint8),mode='L')
    jr.save(out_dir+'/Landsat_red.bmp')
    #dst_ds.GetRasterBand(1).WriteArray(red.astype(byte))

    print ('Now scaling green band.... \n\n\n')
    #ary_scaled=log_bright(pan_sharp_green)
    #green=ary_scaled/np.max(ary_scaled)*255.0
    # Another , probably better, way to brighten images
    green = gamma(np.clip(pan_sharp_green.astype(float)/np.max(pan_sharp_green.astype(float))*255.0,0,255), im_gamma_green)
    jg=Image.fromarray(green.astype(np.uint8),mode='L')
    jg.save(out_dir+'/Landsat_green.bmp')
    #dst_ds.GetRasterBand(2).WriteArray(green.astype(byte))

    print ('Now scaling blue band.... \n\n\n')
    #ary_scaled=log_bright(pan_sharp_blue)
    #blue=ary_scaled/np.max(ary_scaled)*255.0
    # Another , probably better, way to brighten images
    blue = gamma(np.clip(pan_sharp_blue.astype(float)/np.max(pan_sharp_blue.astype(float))*255.0, 0,255), im_gamma_blue)
    jb=Image.fromarray(blue.astype(np.uint8),mode='L')
    jb.save(out_dir+'/Landsat_blue.bmp')
    #dst_ds.GetRasterBand(3).WriteArray(blue.astype(byte))

else:
	
    rad_scl_fact=10

    print ('Now Scaling red band....\n\n\n')
    #ary_scaled=log_bright(pan_sharp_red*rad_scl_fact)
    #red=ary_scaled/np.max(ary_scaled)*255.0

    red = gamma(np.clip(pan_sharp_red/np.max(pan_sharp_red)*255.0, 0,255), im_gamma_red)        

    jr=Image.fromarray(red.astype(np.uint8))
    jr.save(out_dir+'/Landsat_red.bmp')
    #dst_ds.GetRasterBand(1).WriteArray(red.astype(byte))

    print ('Now scaling green band.... \n\n\n')
    #ary_scaled=log_bright(pan_sharp_green*rad_scl_fact)
    #green=ary_scaled/np.max(ary_scaled)*255.0

    green = gamma(np.clip(pan_sharp_green/np.max(pan_sharp_green)*255.0,0,255), im_gamma_green)

    jg=Image.fromarray(green.astype(np.uint8))
    jg.save(out_dir+'/Landsat_green.bmp')
    #dst_ds.GetRasterBand(2).WriteArray(green.astype(byte))

    print ('Now scaling blue band.... \n\n\n')
    #ary_scaled=log_bright(pan_sharp_blue*rad_scl_fact)
    #blue=ary_scaled/np.max(ary_scaled)*255.0

    blue = gamma(np.clip(pan_sharp_blue/np.max(pan_sharp_blue)*255.0, 0,255), im_gamma_blue)

    jb=Image.fromarray(blue.astype(np.uint8))
    jb.save(out_dir+'/Landsat_blue.bmp')
    #dst_ds.GetRasterBand(3).WriteArray(blue.astype(byte))

# Use PIL to merge the RGB, adjust contrast and save the image
# Or you can use imagemagick (the system calls)


imrgb = Image.merge('RGB', (jr,jg,jb))
contrast=ContEnh.Contrast(imrgb,c_mid)
imrgb_en=contrast.enhce(c_enh)
#sharp=ImageEnhance.Sharpness(imrgb_en)
#imrgb_en_sh=sharp.enhance(s_enh)
#imrgb_en.save(out_image_file,quality=95)
imrgb_en.save(out_image_name)


#print 'Compositing to RGB and sharpening Image\n\n' 
#os.system("composite -compose CopyGreen "+out_dir+"/Landsat_green.bmp "+out_dir+"/Landsat_red.bmp "+out_dir+"/Landsat_red_green.bmp")
#os.system("composite -compose CopyBlue "+out_dir+"/Landsat_blue.bmp "+out_dir+"/Landsat_red_green.bmp  "+out_dir+"/Landsat_RGB.bmp")
#os.system("convert -quality 100 -modulate 105,125 -sharpen 3  "+out_dir+"/Landsat_RGB.bmp -compress none bmp3:"+out_image_name)

## This does a few other things once the image is produced. None of this is needed but its all controlled by switches 

try:
  trans
except NameError:

  print ('\n\n\n Leaving image boarder black \n\n\n')
 
else:
  print ('\n\n\n Making the image boarder transparent \n\n\n') 
  os.system("mogrify -transparent-color black -transparent black "+out_image_name)

try:
    logo
except NameError:

    print ('\n\n\n No logo \n\n\n')

else:

    if format=='tif':

        print ('\n\n\n Logo adding function does not work for geotiff. Sorry! \n\n\n')

    else:

        print ('\n\n\n Logo added to bottom left hand corner of the screen \n\n\n')

        os.system("convert "+out_image_name+" -gravity SouthWest "+logo_file+" -compose Over -composite "+out_image_name)


## Convert output image png into other formats, as required (except geotiff which is done above)

try:
    format
except NameError:

    print ('Outputting final image to '+out_image_name+' as no option was selected.\n\n')
    os.system("rm -fr "+out_image_name_tif)

else:

    if format=='bmp':

        print ('Outputting final image to '+out_image_name+'\n\n')
        os.system("rm -fr "+out_image_name_tif)

    elif format=='png':
	
        print ('Outputting final image to '+out_image_name_png+'\n\n')
        os.system("convert "+out_image_name+" -compress none png:"+out_image_name_png)
        os.system("rm -fr "+out_image_name)
        os.system("rm -fr "+out_image_name_tif)

    elif format=='jpeg':

        print ('Outputting final image to '+out_image_name_jpeg+'\n\n')
        os.system("convert "+out_image_name+" -compress none jpeg:"+out_image_name_jpeg)
        os.system("rm -fr "+out_image_name)
        os.system("rm -fr "+out_image_name_tif)

    elif format=='tif':

        print ('Outputting final image to '+out_image_name_tif+'\n\n')
        os.system("rm -fr "+out_image_name)

    else:

        print ('Outputting final image to '+out_image_name+' as your option was not recognised.\n\n')
        os.system("rm -fr "+out_image_name_tif)

try:
  keep
except NameError:

  print ('\n\n\n Deleting all of the interim files and directories. \n\n\n')
  os.system("rm -fr "+out_dir_top)

else:

  print ('\n\n\n Interim processing files have not been deleted. \n\n\n')


ts_1 = time.time()

time_int=(ts_1-ts_0)/60.0

elapsed_time="%.3f" % time_int

sys.exit('Finished True Colour generation from Landsat8 file in '+elapsed_time+' mins')

