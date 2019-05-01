#### Import Many python modules

import matplotlib
matplotlib.use("agg")
import getopt, sys, math, os, glob, fnmatch
import os.path
import struct
from pylab import *
from scipy import *
from PIL import Image, ImageEnhance, ImageFilter, ImageStat
from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.patches as patches
import netCDF4

# Local module
import ContEnh

# Apply a gamma transform to each colour band individually
# This should make the darker blue bit preferentially brighter
# and brighten the whole image
im_gamma_blue=6.5
im_gamma_green=6.5
im_gamma_red=6.5
im_gamma=6.5
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
c_mid=170
s_enh=1.5


## Use with a user defined area to subset. This matches the IMERG products
def gdal_sub(in_netcdf_file,chan_name,out_netcdf_file,radar_tr,lats,lons):

    """!@brief chops out an area of himawari8 obs or RADAR rainfall data
    @param in_nc_file float/int: a netcdf file
    @param out_path string: location that the output files are being written to
    @retval out_nc_file: float/int: a netcdf file (only smaller or bigger and different looking)
    """

    ## Australian RADAR netCDF projection
    radar_proj="+proj=aea +lat_1=-18.0 +lat_2=-36.0 +lon_0=132 +lat_0=0 +a=6378137 +b=6356752"

    ## Convert lats and longs to eastings and northings
    radproj=pyproj.Proj(radar_proj)
    eastings, northings = radproj(lons, lats)

    radar_te = str(eastings[0]) + " " + str(northings[0]) + " " + str(eastings[1]) + " " + str(northings[1])

    out_dir='/short/er8/mab573/AHI/netCDF/test/'
    temp_map=out_dir+'/map.tif'

    # warp the satellite to the radar grid
    os.system('gdalwarp NETCDF:'+in_netcdf_file+':'+chan_name+' '+temp_map+' -ot Float32 -co \"COMPRESS=DEFLATE\" -t_srs \"'+radar_proj+'\" -tr '+radar_tr+' -te '+radar_te+' -overwrite >/dev/null' )

    # Convert back to netcdf
    os.system('gdal_translate ' + temp_map + ' ' + out_netcdf_file + ' -of NETCDF >/dev/null')

def getData(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][0,:,:]
    nc.close()
    return v

def gamma(ary,bright):


    ary_scaled=((ary / 255.0) ** (1.0 /bright))*255.0

    return ary_scaled

# Subset the shit out of everything using a lat long box (bottom left and top right)
#MT Sinobung
#lats=[2.0,6.0]
#lons=[96.0,100.0]

#Gulf of Carp. Flooding
#lats=[-19.7,-12.5]
#lons=[135.0,144.0]

#Melbourne
#lats=[-40.,-36.0]
#lons=[142.9,146.9]

#20171122
#lats=[-40.,-36.0]
#lons=[146.,152.]

#20171206
#lats=[-34.,-31.]
#lons=[150.,154.]

#20180307
#lats=[-39.,-35.]
#lons=[135.5,142.]

lats=[-47.,-10.]
lons=[110.,160.]


tr='500 500'

input_files=[]
#in_dir='/short/er8/mab573/AHI/netCDF/test/'
#in_dir='/short/er8/vov548/mark/helen/20171122/'
#in_dir='/short/er8/vov548/mark/helen/20171206'
in_dir='/short/er8/mab573/AHI/netCDF/2019/02'
for file_data in glob.glob(in_dir + "/*-P1S-ABOM_CREFL_B01-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc"):
    input_files.append(file_data)

input_files.sort()

in_dir_t='/short/er8/mab573/AHI/netCDF/2019/02'

for aa in range(0,len(input_files)):

    in_dir, f_name = os.path.split(input_files[aa])
    splt_f_name=f_name.split('-')

    # Check to see if image file already exists and step to next file if it does
    out_image_file=in_dir_t+"/images/"+splt_f_name[0]+"_500-HIMAWARI8-AHI_Ray_BOA_RGB.png"
    if os.path.isfile(out_image_file): continue

    in_netcdf_file=in_dir+'/'+splt_f_name[0]+'-P1S-ABOM_CREFL_B03-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc'
    if not os.path.isfile(in_netcdf_file): continue
    chan_name='channel_0003_corrected_reflectance'
    out_netcdf_file_B3_1000='/short/er8/mab573/AHI/netCDF/test/B3_1000m.nc'
    gdal_sub(in_netcdf_file,chan_name,out_netcdf_file_B3_1000,tr,lats,lons)
    B3_1000_m=getData(out_netcdf_file_B3_1000,chan_name)

    in_netcdf_file=in_dir+'/'+splt_f_name[0]+'-P1S-ABOM_CREFL_B01-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc'
    if not os.path.isfile(in_netcdf_file): continue
    chan_name='channel_0001_corrected_reflectance'
    out_netcdf_file_B1_500='/short/er8/mab573/AHI/netCDF/test/B1_500m.nc'
    gdal_sub(in_netcdf_file,chan_name,out_netcdf_file_B1_500,tr,lats,lons)
    B1_1000_m=getData(out_netcdf_file_B1_500,chan_name)

    in_netcdf_file=in_dir+'/'+splt_f_name[0]+'-P1S-ABOM_CREFL_B02-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc'
    if not os.path.isfile(in_netcdf_file): continue
    chan_name='channel_0002_corrected_reflectance'
    out_netcdf_file_B2_500='/short/er8/mab573/AHI/netCDF/test/B2_500m.nc'
    gdal_sub(in_netcdf_file,chan_name,out_netcdf_file_B2_500,tr,lats,lons)
    B2_1000_m=getData(out_netcdf_file_B2_500,chan_name)

    in_netcdf_file=in_dir+'/'+splt_f_name[0]+'-P1S-ABOM_CREFL_B03-PRJ_GEOS141_500-HIMAWARI8-AHI.nc'
    if not os.path.isfile(in_netcdf_file): continue
    chan_name='channel_0003_corrected_reflectance'
    out_netcdf_file_B3_500='/short/er8/mab573/AHI/netCDF/test/B3_500m.nc'
    gdal_sub(in_netcdf_file,chan_name,out_netcdf_file_B3_500,tr,lats,lons)
    B3_500_m=getData(out_netcdf_file_B3_500,chan_name)

    in_netcdf_file=in_dir+'/'+splt_f_name[0]+'-P1S-ABOM_CREFL_B04-PRJ_GEOS141_1000-HIMAWARI8-AHI.nc'
    if not os.path.isfile(in_netcdf_file): continue
    chan_name='channel_0004_corrected_reflectance'
    out_netcdf_file_B4_500='/short/er8/mab573/AHI/netCDF/test/B4_500m.nc'
    gdal_sub(in_netcdf_file,chan_name,out_netcdf_file_B4_500,tr,lats,lons)
    B4_1000_m=getData(out_netcdf_file_B4_500,chan_name)


    # Generate the Ratio band (B3_1000/B3_500)

    R=B3_1000_m/B3_500_m

    # Sharpen the other bands

    B1_500_m=B1_1000_m/R
    B2_500_m=B2_1000_m/R
    B4_500_m=B4_1000_m/R

    # Colour correct green band
    green=(0.93*B2_500_m)+(0.07*B4_500_m)

    # Make some images

    #Blue band
    np.clip(B1_500_m, clip_level_min, clip_level_max, out=B1_500_m)
    arr = (B1_500_m - clip_level_min) / clip_level_max * 255.0
    band_scaled = gamma(arr, im_gamma_blue)
    out_image_blue = in_dir_t + '/BOA_blue.png'
    jb = Image.fromarray(band_scaled[::-1].astype(np.uint8), mode='L')
    jb.save(out_image_blue)

    #Red band
    np.clip(B3_500_m, clip_level_min, clip_level_max, out=B3_500_m)
    arr = (B3_500_m - clip_level_min) / clip_level_max * 255.0
    band_scaled = gamma(arr, im_gamma_red)
    out_image_red = in_dir_t + '/BOA_red.png'
    jr = Image.fromarray(band_scaled[::-1].astype(np.uint8), mode='L')
    jr.save(out_image_red)

    #Green band
    np.clip(green, clip_level_min, clip_level_max, out=green)
    arr = (green - clip_level_min) / clip_level_max * 255.0
    band_scaled = gamma(arr, im_gamma_green)
    out_image_green = in_dir_t + '/BOA_green.png'
    jg = Image.fromarray(band_scaled[::-1].astype(np.uint8), mode='L')
    jg.save(out_image_green)

    # Adjust the image and output
    imrgb = Image.merge('RGB', (jr,jg,jb))
    #imrgb.save(out_image_file)
    contrast=ContEnh.Contrast(imrgb,c_mid)
    imrgb_en=contrast.enhce(c_enh)
    imrgb_en.save(out_image_file)

    os.system("rm "+out_image_blue)
    os.system("rm "+out_image_green)
    os.system("rm "+out_image_red)
    os.system("rm "+out_netcdf_file_B3_1000)
    os.system("rm "+out_netcdf_file_B1_500)
    os.system("rm "+out_netcdf_file_B2_500)
    os.system("rm "+out_netcdf_file_B3_500)
    os.system("rm "+out_netcdf_file_B4_500)

