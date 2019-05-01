
import getopt, sys, math, os, glob, fnmatch
import numpy as np
import matplotlib
matplotlib.use('Agg')
import netCDF4
from pylab import *
from scipy import *
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
import timeit
import pickle

## Import local modules
#sys.path.insert(0, '/flurry/home/mbroom/src/AHI_code/True_colour/')
#sys.path.insert(0, '/flurry/home/mbroom/src/AHI_code/True_colour/AHI_look_up_generation/')
#sys.path.insert(0, '/flurry/home/mbroom/src/AHI_code/True_colour/new_code/')

import read_MODTRAN_lut_V2

def usage():

 print "SYNOPSIS:\n This program produces interpolated RTC data from MODTRAN LUTs"
 print "\nUSAGE: Create_2000m_interp_RTC_NCI [OPTIONS] "
 print "\nOPTIONS:\n"
 print "\n-h, Help or usage"
 print "\n-a, Directory location of the ancillary AHI data"
 print "\n-i, Directory location of the AHI data to be processed"
 print "\n-o, Output location for data/products to be written to"
 print "Base Usage:"
 print "python Create_2000m_interp_RTC_NCI -a /Dir/sub_dir/... -i /Dir/sub_dir/... -o /Dir/sub_dir/ \n\n"
 sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:a:o:", ["help", "in_dir", "anc_dir", "out_dir"])

except getopt.GetoptError, err:

    # print help information and exit:
    print str(err)  # will print something like "option -a not recognized"
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
    elif o in ("-o", "--out_dir"):
        out_dir = a
    else:
        assert False, "unhandled option"


def Write_RTC_NetCDF(filename,Lp,Eg,Tup,S):

    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(Lp[:, 0]))
    nc.createDimension("y_dim", len(Lp[0, :]))
    v = nc.createVariable("Lp_0", float, ("x_dim", "y_dim"))
    v[:, :] = Lp
    v = nc.createVariable("Eg_0", float, ("x_dim", "y_dim"))
    v[:, :] = Eg
    v = nc.createVariable("T_up", float, ("x_dim", "y_dim"))
    v[:, :] = Tup
    v = nc.createVariable("S", float, ("x_dim", "y_dim"))
    v[:, :] = S
    nc.close()
    return

def Write_T_up_NetCDF(filename,T_up):

    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(T_up[:, 0]))
    nc.createDimension("y_dim", len(T_up[0, :]))
    v = nc.createVariable("T_up", float, ("x_dim", "y_dim"))
    v[:, :] = T_up
    return

def Write_S_NetCDF(filename,S):

    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(S[:, 0]))
    nc.createDimension("y_dim", len(S[0, :]))
    v = nc.createVariable("S", float, ("x_dim", "y_dim"))
    v[:, :] = S
    return

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

out_file_base='/short/er8/mab573/AHI/RTC/'

## Used to calculate run time
start = timeit.default_timer()


for solar_file in glob.glob(in_dir + "/*-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-*.nc"):

    for ancillary_file in glob.glob(anc_dir+"/*-P1S-ABOM_GEOM_SENSOR-PRJ_GEOS141_2000-*.nc"):

        out_path = check_outdir(solar_file, out_file_base)

        print 'Processing band 1 ...\n'
        Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 0)
        file_name=out_path+'/Band1_2000_RTC_interpolated_values.nc'
        Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

        print 'Processing band 2 ...\n'
        Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 1)
        file_name = out_path + '/Band2_2000_RTC_interpolated_values.nc'
        Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

        print 'Processing band 3 ...\n'
        Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 2)
        file_name = out_path + '/Band3_2000_RTC_interpolated_values.nc'
        Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

        print 'Processing band 4 ...\n'
        Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 3)
        file_name = out_path + '/Band4_2000_RTC_interpolated_values.nc'
        Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

stop = timeit.default_timer()
print 'Completed RTC interpolation to 2000m resolution in', stop - start


'''
### When running this initially on a new LUT, run this block to interpolate T_up and S.
### All other situations use the above block as T_up and S only needs doing once.

for ancillary_file in glob.glob(anc_dir+"/*-P1S-ABOM_GEOM_SENSOR-PRJ_GEOS141_2000-*.nc"):

    for solar_file in glob.glob(in_dir + "/*-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-*.nc"):

        print 'Processing band 1 ...\n'
        S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 0)
        #file_name=anc_dir+'/Band1_2000m_T_up_interpolated_values.nc'
        #Write_T_up_NetCDF(file_name, T_up)
        file_name = anc_dir + '/Band1_2000m_S_interpolated_values.nc'
        Write_S_NetCDF(file_name, S)

        print 'Processing band 2 ...\n'
        S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 1)
        #file_name = anc_dir + '/Band2_2000m_T_up_interpolated_values.nc'
        #Write_T_up_NetCDF(file_name, T_up)
        file_name = anc_dir + '/Band2_2000m_S_interpolated_values.nc'
        Write_S_NetCDF(file_name, S)

        print 'Processing band 3 ...\n'
        S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 2)
        #file_name = anc_dir + '/Band3_2000m_T_up_interpolated_values.nc'
        #Write_T_up_NetCDF(file_name, T_up)
        file_name = anc_dir + '/Band3_2000m_S_interpolated_values.nc'
        Write_S_NetCDF(file_name, S)

        print 'Processing band 4 ...\n'
        S= read_MODTRAN_lut_V2.Interp_LUT_stuff(solar_file, ancillary_file, 3)
        #file_name = anc_dir + '/Band4_2000m_T_up_interpolated_values.nc'
        #Write_T_up_NetCDF(file_name, T_up)
        file_name = anc_dir + '/Band4_2000m_S_interpolated_values.nc'
        Write_S_NetCDF(file_name, S)

stop = timeit.default_timer()
print 'Completed RTC interpolation to 2000m resolution in', stop - start
'''

