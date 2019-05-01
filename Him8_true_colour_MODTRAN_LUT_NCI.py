#!/usr/bin/env python

'''
This will perform Rayleigh correction of Himawari 8 AHI bands 1 - 4 so that these can be used to produce
a true-colour image. This uses a lookup table produced by M. Broomhall using MODTRAN5.
Once the Rayleigh correction is done a colour transform is done on band 2 (green band). This takes the band
centred on ~510 nm and infers a band centred on 532 nm. This creates images that are 'less brown' and resemble
more pleasing and closer to expected images. This transformation is taken from Ref (1)

(3) Miller, S., Schmit, T., Seaman, C., Lindsey, D., Gunshor, M., Sumida, Y., Hillger, D., 2016, 'A sight for sore eyes-
the return of true color to geostationary satellites', Bull. AMer. Meteor. Soc.

'''


import getopt, sys, math, os, glob, fnmatch
import numpy as np
import matplotlib
matplotlib.use('Agg')
import netCDF4
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import timeit
from PIL import Image, ImageEnhance, ImageFilter, ImageStat
import operator

##--------- Local programs/functions -------------
import cloud_height_correction_V2 as chc
import ContEnh


## Used to calculate run time
start = timeit.default_timer()

def usage():

 print "SYNOPSIS:\n This program produces true-colour images from Himawari8 (AHI) data using a MODTRAN5-derived lookup table"
 print "\nUSAGE: Him8_true_colour_Yi_LUT [OPTIONS] "
 print "\nOPTIONS:\n"
 print "\n-h, Help or usage"
 print "\n-a, Directory location of the ancillary AHI data"
 print "\n-i, Directory location of the AHI data to be processed"
 print "\n-o, Directory location to write images to"
 print "Base Usage:"
 print "python Him8_true_colour_Yi_LUT -a /Dir/sub_dir/... -i /Dir/sub_dir/...\n\n"
 sys.exit()

try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:a:r:o:", ["help", "in_dir","anc_dir","rtc_dir","out_dir"])

except getopt.GetoptError, err:

        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-a", "--anc_dir"):
        anc_dir = a
    elif o in ("-i", "--in_dir"):
        in_dir = a
    elif o in ("-r", "--rtc_dir"):
        rtc_dir = a
    elif o in ("-o", "--out_dir"):
        out_dir = a
    else:
        assert False, "unhandled option"


#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>********************

######################## Constant values #####################################
# Scale value used fro the Miller et al. green band transformation
F = 0.07

# Apply a gamma transform to each colour band individually
# This should make the darker blue bit preferentially brighter
# and brighten the whole image
im_gamma_blue=6.7
im_gamma_green=6.7
im_gamma_red=6.7
im_gamma=6.7
#im_gamma=10.
ir_gamma=0.32

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


### Band calibration coefficients
### These will convert the TOA Albedo? to radiance
### Units are W.m^-2.sr^-1.um (although this is a band equivalent value)
AHI_Cal=[0.0015588241,0.001611941,0.0019254997,0.0032324977,0.0129632414,0.0418124877]

## The spherical albedo is directly related to the solar zenith angle (for each band)
## but varies only by a very small amount and probably not by enough to require using an interpolated
## field in the calculation. These values are close to the maximum expected with the lookup table
## produced in June 2016 from MODTRAN. May need to look at these again if the LUT changes.
S = [0.14319, 0.11086, 0.04754, 0.01559]

#%%%%%%%%^^^^^^^^^^%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%
# Functions for doing clever stuff

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

# simple routine to read a variable from a netcdf file
def read_nc_array(file_name, var_name, selection_slice=None):
    """
    @brief Read a netcdf variable and return data and attributes
    """
    root_group = netCDF4.Dataset(file_name)
    if var_name not in root_group.variables.keys():
        root_group.close()
        raise Exception('Variable (%s) not found' % var_name)
    var = root_group.variables[var_name]
    attr_names = var.ncattrs()
    attrs = {}
    for attr_name in attr_names:
        attrs[attr_name] = var.getncattr(attr_name)
    if selection_slice is None:
        data = var[:]
    else:
        data = var[selection_slice]
    root_group.close()
    return data, attrs

def getData(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][0,:,:]
    nc.close()
    return v

def getNav(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][0,:,:]
    nc.close()
    return v

def getTime(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][:]
    nc.close()
    return v

def get_RTC(fileName, varName):

    nc = netCDF4.Dataset(fileName)
    if varName not in nc.variables:
        sys.exit("Var "+varName+" not in file "+fileName)

    v = nc.variables[varName][:,:]
    nc.close()
    return v

def Write_RTC_NetCDF(filename,Lp,Eg):

    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(Lp[:, 0]))
    nc.createDimension("y_dim", len(Lp[0, :]))
    v = nc.createVariable("Lp_0", float, ("x_dim", "y_dim"))
    v[:, :] = Lp
    v = nc.createVariable("Eg_0", float, ("x_dim", "y_dim"))
    v[:, :] = Eg
    nc.close()
    return

def EarthSunD(doy):
    ESD = 1 - 0.01672 * np.cos(np.radians(0.9856 * (doy - 4)))
    return ESD

def DayOfYear(year,month,day):
    d = datetime.date(year, 1, 1)
    t_d=datetime.date(year,month,day)
    doy=t_d-d
    return doy.days

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

#%%%%%%%%^^^^^^^^^^%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%^^^^^^^^^^^%%%%%%%%%%%
#### Main program that does the bizness ####
'''
#### This block of code makes TOA RGB images with the green colour correction
for ancillary_file in glob.glob(anc_dir+"/*-P1S-ABOM_GEOM_SENSOR-PRJ_GEOS141_2000-*.nc"):

    sen_zen = getNav(ancillary_file, 'sensor_zenith_angle')
    print 'Found an ancillary file'

for file_data in glob.glob(in_dir + "/*-P1S-ABOM_OBS_B0[1234]-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):
    if 'B01' in file_data:
        dat_file_split = file_data.split('-')
        solar_file = dat_file_split[0] + "-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"
        sol_zen = getNav(solar_file, 'solar_zenith_angle')
        TOA = getData(file_data, 'channel_0001_scaled_radiance')
        #TOA[sen_zen > 85.0] = 0.0
        sol_zen[sol_zen > 85.0] = 85.0
        band=TOA*AHI_Cal[0]
        arr = band / np.max(band) * 255
        band_scaled = gamma(arr, 3.0)
        plt.imsave(fname='TOA_blue.png', arr=band_scaled, dpi=300,cmap='gray')
    elif 'B03' in file_data:
        dat_file_split = file_data.split('-')
        solar_file = dat_file_split[0] + "-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"
        sol_zen = getNav(solar_file, 'solar_zenith_angle')
        TOA = getData(file_data, 'channel_0003_scaled_radiance')
        #TOA[sen_zen > 85.0] = 0.0
        sol_zen[sol_zen > 85.0] = 85.0
        band = TOA*AHI_Cal[2]
        arr = band / np.max(band) * 255
        band_scaled = gamma(arr, 3.0)
        plt.imsave(fname='TOA_red.png', arr=band_scaled, dpi=300,cmap='gray')
    elif 'B04' in file_data:
        dat_file_split = file_data.split('-')
        solar_file = dat_file_split[0] + "-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"
        sol_zen = getNav(solar_file, 'solar_zenith_angle')
        TOA = getData(file_data, 'channel_0004_scaled_radiance')
        #TOA[sen_zen > 85.0] = 0.0
        sol_zen[sol_zen > 85.0] = 85.0
        NIR_band = TOA*AHI_Cal[3]
    elif 'B02' in file_data:
        dat_file_split = file_data.split('-')
        solar_file = dat_file_split[0] + "-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"
        sol_zen = getNav(solar_file, 'solar_zenith_angle')
        TOA = getData(file_data, 'channel_0002_scaled_radiance')
        #TOA[sen_zen > 85.0] = 0.0
        sol_zen[sol_zen > 85.0] = 85.0
        Green_band = TOA*AHI_Cal[1]
    else:
        print "Nothing to see here"

green = (1 - F) * (Green_band) + F * (NIR_band)
arr = green / np.max(green) * 255
band_scaled = gamma(arr, 3.0)

#green = (1 - F) * (Green_band) + F * (NIR_band)
#arr = green / np.max(green) * 255
#band_scaled = log_bright(arr, 0.5)

plt.imsave(fname='TOA_green.png', arr=band_scaled, dpi=300,cmap='gray')

os.system("composite -compose CopyGreen TOA_green.png TOA_red.png red_green.png")
os.system("composite -compose CopyBlue TOA_blue.png red_green.png TOA_MODTRAN_RGB.png")
#os.system("convert Test_MODTRAN_RGB.png -transparent black AHI_BOA_RGB.png")

sys.exit()
'''


## Grab the ancillary data (these parameters are constant for all scenes) so it only needs doing once
for ancillary_file in glob.glob(anc_dir+"/*-P1S-ABOM_GEOM_SENSOR-PRJ_GEOS141_2000-*.nc"):

    sen_zen = getNav(ancillary_file, 'sensor_zenith_angle')
    print 'Found an ancillary file....\n'

for solar_file in glob.glob(in_dir + "/*-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):

    sol_zen = getNav(solar_file, 'solar_zenith_angle')
    print 'Found a solar file....\n'


# Calculate cloud height/temp correction to reduce the path radiance for high clouds
Lp_0_scale=chc.ray_scale(in_dir)

## Open fiddle create and blend out IR imagery
for file_data in glob.glob(in_dir + "/*-P1S-ABOM_OBS_B13-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):
    if 'B13' in file_data:

        out_path=check_outdir(file_data,out_dir)

        CT = getData(file_data, 'channel_0013_brightness_temperature')
        blend = (88.0 - sol_zen) / 10.0
        blend[blend < 0.0] = 0.0
        blend[blend > 1.0] = 1.0

        #ir=255.0-np.clip(2*(CT-180), 0.0, 255.0)
        CT=np.clip(CT,170.0,340.0)
 
        ir=200.0-np.clip(((CT-170.0)/170.0*200.0), 0.0, 200.0)

        ir = (ir * (1.0 - blend))
        ir[(sen_zen > 85.0) | (sol_zen < 78.0)] = 0.0

        # Create IR sandwich product

        #CT=np.clip(CT,180.0,240.0)
        #arr = CT / np.max(CT) * 255.0
        #arr[(sen_zen > 85.0) | (sol_zen > 88.0)] = 255.0
        
        #plt.imsave(fname=out_path+'/B13_BT_sandwich.png', arr=arr, dpi=300,cmap='spectral_r')
        #os.system("convert "+out_path+"/B13_BT_sandwich.png -transparent black "+out_path+"/B13_BT_sandwich_trans.png")
        #CT = getData(file_data, 'channel_0013_brightness_temperature')


for file_data in glob.glob(in_dir + "/*-P1S-ABOM_OBS_B01-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):
    if 'B01' in file_data:

        print 'Now processing band 1....'

        rtc_path = check_outdir(file_data, rtc_dir)
        print rtc_path

        fileName = rtc_path + '/Band1_2000_RTC_interpolated_values.nc'
        Lp_0 = get_RTC(fileName, 'Lp_0')
        Eg_0 = get_RTC(fileName, 'Eg_0')
        T_up = get_RTC(fileName, 'T_up')
        S = get_RTC(fileName, 'S')

        TOA = getData(file_data, 'channel_0001_scaled_radiance')

        ## Parse up file_data name to get date
        splt_file = file_data.split('/')
        f_name = splt_file[len(splt_file) - 1]
        splt_f_name = f_name.split('-')
        d_time = splt_f_name[0]
        s_year = d_time[0:4]
        s_month = d_time[4:6]
        s_day = d_time[6:8]
        esd = EarthSunD(DayOfYear(int(s_year), int(s_month), int(s_day)))

        A = (math.pi * (((TOA / AHI_Cal[0])/10.0) - (Lp_0*Lp_0_scale))) / (Eg_0 * T_up)
        band = A / (1 + (A * S))

        band[(sen_zen > 85.0) | (sol_zen > 88.0)] = 0.0
        blend = (sol_zen - 78.0) / 10.0
        blend[blend < 0.0] = 0.0
        blend[blend > 1.0] = 1.0
        blend2 = (sen_zen - 75.0) / 10.0
        blend2[blend2 < 0.0] = 0.0
        blend2[blend2 > 1.0] = 1.0

        band = (band * (1.0 - blend))
        band = (band * (1.0 - blend2))

        out_path=check_outdir(file_data,out_dir)

        int_band = np.around((band * 10000.0), decimals=0)
        filename = out_path + '/Blue.nc'
        nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
        nc.createDimension("x_dim", len(band[:, 0]))
        nc.createDimension("y_dim", len(band[0, :]))
        v = nc.createVariable("Blue", float, ("x_dim", "y_dim"))
        #v[:, :] = int_band
        v[:, :] = ir_scale
        nc.close()

        print 'Maximum blue reflectance value is ',np.max(band)

        np.clip(band, clip_level_min, clip_level_max, out=band)
        arr = (band-clip_level_min) / clip_level_max * 255.0
        #np.clip(band, clip_level_min, np.max(band), out=band)
        #arr = (band-clip_level_min) / np.max(band) * 255.0
        
        band_scaled = gamma(np.clip((arr+gamma(ir,ir_gamma)),0,255), im_gamma_blue)
        
        out_image_blue=out_path + '/BOA_blue.png'

        jb=Image.fromarray(band_scaled.astype(np.uint8),mode='L')
        #jb=Image.fromarray(scaled_gamma.astype(np.uint8),mode='L')
        #jb=Image.fromarray(arr.astype(np.uint8),mode='L')
        #jb=Image.fromarray(np.clip((arr+ir).astype(np.uint8),0,255),mode='L')
        jb.save(out_image_blue)
        
        #arr[(CT < 180.0) | (CT > 240.0)] = 0.0
        #cb=Image.fromarray(arr.astype(np.uint8),mode='L')
        #cb.save(out_path + '/blue_cloud.png')

        

for file_data in glob.glob(in_dir + "/*-P1S-ABOM_OBS_B03-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):
    if 'B03' in file_data:

        print '\nNow processing band 3....'

        rtc_path = check_outdir(file_data, rtc_dir)

        fileName = rtc_path + '/Band3_2000_RTC_interpolated_values.nc'

        Lp_0 = get_RTC(fileName, 'Lp_0')
        Eg_0 = get_RTC(fileName, 'Eg_0')
        T_up = get_RTC(fileName, 'T_up')
        S = get_RTC(fileName, 'S')

        TOA = getData(file_data, 'channel_0003_scaled_radiance')

        ## Parse up file_data name to get date
        splt_file = file_data.split('/')
        f_name = splt_file[len(splt_file) - 1]
        splt_f_name = f_name.split('-')
        d_time = splt_f_name[0]
        s_year = d_time[0:4]
        s_month = d_time[4:6]
        s_day = d_time[6:8]

        A = (math.pi * (((TOA / AHI_Cal[2])/10.0) - (Lp_0*Lp_0_scale))) / (Eg_0 * T_up)
        band = A / (1 + (A * S))
        
        ## This should blend the TOA into the BOA from 75 - 89 sol zenith

        band[(sen_zen > 85.0) | (sol_zen > 88.0)] = 0.0
        blend = (sol_zen - 78.0) / 10.0
        blend[blend < 0.0] = 0.0
        blend[blend > 1.0] = 1.0
        blend2 = (sen_zen - 75.0) / 10.0
        blend2[blend2 < 0.0] = 0.0
        blend2[blend2 > 1.0] = 1.0

        band = (band * (1.0 - blend))
        band = (band * (1.0 - blend2))

        out_path=check_outdir(file_data,out_dir)

        int_band = np.around((band * 10000.0), decimals=0)
        filename = out_path + '/Red.nc'
        nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
        nc.createDimension("x_dim", len(band[:, 0]))
        nc.createDimension("y_dim", len(band[0, :]))
        v = nc.createVariable("Blue", int, ("x_dim", "y_dim"))
        v[:, :] = int_band
        nc.close()

        print 'Maximum red reflectance value is ', np.max(band)

        np.clip(band, clip_level_min, clip_level_max, out=band)
        arr = (band-clip_level_min) / clip_level_max * 255.0
        #np.clip(band, clip_level_min, np.max(band), out=band)
        #arr = (band-clip_level_min) / np.max(band) * 255.0
        #band_scaled = gamma(arr, im_gamma_red)
        band_scaled = gamma(np.clip((arr+gamma(ir,ir_gamma)),0,255), im_gamma_red)

        out_image_red = out_path + '/BOA_red.png'

        jr=Image.fromarray(band_scaled.astype(np.uint8),mode='L')
        #jr=Image.fromarray(arr.astype(np.uint8),mode='L')
        #jr=Image.fromarray((arr+ir).astype(np.uint8),mode='L')
        jr.save(out_image_red)

        #arr[(CT < 180.0) | (CT > 240.0)] = 0.0
        #arr=gamma(arr,1.5)
        #cr=Image.fromarray(arr.astype(np.uint8),mode='L')
        #cr_imrgb = Image.merge('RGB', (cr,cr,cr))
        #contrast=ImageEnhance.Contrast(cr_imrgb)
        #cr_en=contrast.enhance(2.5)
        #cr_en.save(out_path + '/red_cloud.png')


for file_data in glob.glob(in_dir + "/*-P1S-ABOM_OBS_B02-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):
    if 'B02' in file_data:

        print 'Now processing band 2....\n'

        fileName = rtc_path + '/Band2_2000_RTC_interpolated_values.nc'
        Lp_0 = get_RTC(fileName, 'Lp_0')
        Eg_0 = get_RTC(fileName, 'Eg_0')
        T_up = get_RTC(fileName, 'T_up')
        S = get_RTC(fileName, 'S')

        TOA = getData(file_data, 'channel_0002_scaled_radiance')

        ## Parse up file_data name to get date
        splt_file = file_data.split('/')
        f_name = splt_file[len(splt_file) - 1]
        splt_f_name = f_name.split('-')
        d_time = splt_f_name[0]
        s_year = d_time[0:4]
        s_month = d_time[4:6]
        s_day = d_time[6:8]

        A = (math.pi * (((TOA / AHI_Cal[1])/10.0) - (Lp_0*Lp_0_scale))) / (Eg_0 * T_up)
        band = A / (1 + (A * S))

        ## This should blend the TOA into the BOA from 75 - 89 sol zenith

        out_path = check_outdir(file_data, out_dir)

        filename = out_path+'/Temp_green.nc'
        nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
        nc.createDimension("x_dim", len(band[:, 0]))
        nc.createDimension("y_dim", len(band[0, :]))
        v = nc.createVariable("Temp_green", float, ("x_dim", "y_dim"))
        v[:, :] = band
        nc.close()

for file_data in glob.glob(in_dir + "/*-P1S-ABOM_OBS_B04-PRJ_GEOS141_2000-HIMAWARI8-AHI.nc"):
    if 'B04' in file_data:

        print 'Now processing band 4....\n'

        fileName = rtc_path + '/Band4_2000_RTC_interpolated_values.nc'
        Lp_0 = get_RTC(fileName, 'Lp_0')
        Eg_0 = get_RTC(fileName, 'Eg_0')
        T_up = get_RTC(fileName, 'T_up')
        S = get_RTC(fileName, 'S')


        TOA = getData(file_data, 'channel_0004_scaled_radiance')

        ## Parse up file_data name to get date
        splt_file = file_data.split('/')
        f_name = splt_file[len(splt_file) - 1]
        splt_f_name = f_name.split('-')
        d_time = splt_f_name[0]
        s_year = d_time[0:4]
        s_month = d_time[4:6]
        s_day = d_time[6:8]

        A = (math.pi * (((TOA / AHI_Cal[3])/10.0) - (Lp_0*Lp_0_scale))) / (Eg_0 * T_up)
        band = A / (1 + (A * S))


        ## This should blend the TOA into the BOA from 75 - 89 sol zenith

        out_path = check_outdir(file_data, out_dir)

        filename = out_path+'/Temp_NIR.nc'
        nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
        nc.createDimension("x_dim", len(band[:, 0]))
        nc.createDimension("y_dim", len(band[0, :]))
        v = nc.createVariable("Temp_NIR", float, ("x_dim", "y_dim"))
        v[:, :] = band
        nc.close()

print 'Now adjusting the green band....\n'
fileName = out_path+'/Temp_green.nc'
band2 = get_RTC(fileName, 'Temp_green')

green = ((1 - F) * (band2)) + (F * (band))

green[(sen_zen > 85.0) | (sol_zen > 88.0)] = 0.0
blend = (sol_zen - 78.0) / 10.0
blend[blend < 0.0] = 0.0
blend[blend > 1.0] = 1.0
blend2 = (sen_zen - 75.0) / 10.0
blend2[blend2 < 0.0] = 0.0
blend2[blend2 > 1.0] = 1.0

green = (green * (1.0 - blend))
green = (green * (1.0 - blend2))

out_path=check_outdir(file_data,out_dir)

int_green = np.around((green * 10000.0), decimals=0)
filename = out_path + '/Green.nc'
nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
nc.createDimension("x_dim", len(green[:, 0]))
nc.createDimension("y_dim", len(green[0, :]))
v = nc.createVariable("Green", int, ("x_dim", "y_dim"))
v[:, :] = int_green
nc.close()

print 'Maximum green reflectance value is ',np.max(green)

np.clip(green,clip_level_min, clip_level_max, out=green)
arr = (green-clip_level_min) / clip_level_max * 255.0
#np.clip(green, clip_level_min, np.max(green), out=green)
#arr = (green-clip_level_min) / np.max(green) * 255.0

#band_scaled = gamma(arr, im_gamma_green)
band_scaled = gamma(np.clip((arr+gamma(ir,ir_gamma)),0,255), im_gamma_green)

out_image_green=out_path + '/BOA_green.png'

jg=Image.fromarray(band_scaled.astype(np.uint8),mode='L')
#jg=Image.fromarray(arr.astype(np.uint8),mode='L')
#jg=Image.fromarray((arr+ir).astype(np.uint8),mode='L')
jg.save(out_image_green)


#arr[(CT < 180.0) | (CT > 240.0)] = 0.0
#cg=Image.fromarray(arr.astype(np.uint8),mode='L')
#cg.save(out_path + '/green_cloud.png')

#c_imrgb = Image.merge('RGB', (cr,cg,cb))

#contrast=ImageEnhance.Contrast(c_imrgb)
#c_imrgb_en=contrast.enhance(2.5)

#c_imrgb_en.save(out_path + '/RGB_cloud.png')

in_dir, f_name = os.path.split(file_data)
splt_f_name=f_name.split('-')
#out_image_file=out_path+"/"+splt_f_name[0]+"_2000-HIMAWARI8-AHI_Ray_BOA_RGB.jpg"
out_image_file=out_path+"/"+splt_f_name[0]+"_2000-HIMAWARI8-AHI_Ray_BOA_RGB.png"

#os.system("composite -compose CopyGreen "+out_image_green+" "+out_image_red+" "+out_path+"/red_green.png")
#os.system("composite -compose CopyBlue "+out_image_blue+" "+out_path+"/red_green.png "+out_image_file)
#os.system("convert "+out_image_file+" -sigmoidal-contrast 2,50% "+out_image_file)
#os.system("convert "+out_image_file+" -contrast-stretch 1%x1% "+out_image_file)

imrgb = Image.merge('RGB', (jr,jg,jb))

#imrgb.save(out_image_file)
#os.system("mogrify -transparent-color black -transparent black "+out_image_file)



#imrgb_ga=do_gamma(imrgb, im_gamma)

contrast=ContEnh.Contrast(imrgb,c_mid)

imrgb_en=contrast.enhce(c_enh)

#sharp=ImageEnhance.Sharpness(imrgb_en)
#imrgb_en_sh=sharp.enhance(s_enh)

#imrgb_en.save(out_image_file,quality=95)
imrgb_en.save(out_image_file)

#mean = int(ImageStat.Stat(imrgb_ga.convert("L")).mean[0] + 0.5)
#print 'Mean ', mean


os.system("rm "+out_path+"/Temp_NIR.nc")
os.system("rm "+out_path+"/Temp_green.nc")
#os.system("rm "+out_path+"/BOA_red.png")
#os.system("rm "+out_path+"/BOA_green.png")
#os.system("rm "+out_path+"/BOA_blue.png")
os.system("rm "+out_path+"/red_green.png")

stop = timeit.default_timer()
print 'Program run time is ', stop - start


sys.exit()
