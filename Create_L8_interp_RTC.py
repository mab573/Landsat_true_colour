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
# from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
import pickle

from osgeo import gdal
from osgeo import osr

import read_MODTRAN_lut_L8

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

def Write_NetCDF(filename,data):

    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(data[:, 0]))
    nc.createDimension("y_dim", len(data[0, :]))
    v = nc.createVariable("data", float, ("x_dim", "y_dim"))
    v[:, :] = data
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


#PROCESS THE RED BAND (BAND4). OPEN, ATMOSPERIC COMPENSATION AND WRITE OUT RESULTANT BAND
def open_geotiff(filename):
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    (X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
    srs_wkt = ds.GetProjection()
    Nx = ds.RasterXSize
    Ny = ds.RasterYSize
    ary = []
    ary = ds.GetRasterBand(1).ReadAsArray()
    return ary

def band8_create_raster(filename):

    in_dir, f_name = os.path.split(filename)
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    dsProj = ds.GetProjection()
    dsTrans = ds.GetGeoTransform()
    x = ds.RasterXSize 
    y = ds.RasterYSize
    ary = ds.GetRasterBand(1).ReadAsArray()

    outputfile = in_dir+'/temp.tif'
    driver= gdal.GetDriverByName('GTiff')
    output = driver.Create(outputfile,x*2,y*2,1,gdal.GDT_UInt16)
    output.SetGeoTransform(dsTrans)
    output.SetProjection(dsProj)

    gdal.ReprojectImage(ds,output,dsProj,dsProj,gdal.GRA_Bilinear)

    return outputfile
    

def generate_RTC_rasters(VZA_file,SZA_file,VA_file,SA_file,VZA_file_B8,SZA_file_B8,VA_file_B8,SA_file_B8):

    # Open solar zenith, Satellite view and relative azimuth
    #
    ## Temporary hard coded path
    path='/short/er8/mab573/Landsat/'
    out_file_base=path+'RTC/'
    VZA=open_geotiff(VZA_file_B8)
    SZA=open_geotiff(SZA_file_B8)
    VA=open_geotiff(VA_file_B8)
    SA=open_geotiff(SA_file_B8)
    #RA=open_geotiff(path+'RELATIVE-AZIMUTH.tif')

    # This should make the RA compatible with my table (ra 0 -180)
    RA = abs(abs(abs(VA - SA)-180.0)-180.0)

    #Write_NetCDF(path+'test_RA.nc',RA)
    #Write_NetCDF(path+'test_SA.nc',SA)
    #Write_NetCDF(path+'test_VA.nc',VA)
    #Write_NetCDF(path+'test_SZA.nc',SZA)
    #Write_NetCDF(path+'test_VZA.nc',VZA)
    #out_path = check_outdir(solar_file, out_file_base)

    #print 'Processing band 2 ...\n'
    Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_L8.Interp_LUT_stuff(SZA,VZA,RA, 1)
    file_name=out_file_base+'/Band2_L8_RTC_interpolated_values.nc'
    Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

    #print 'Processing band 3 ...\n'
    Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_L8.Interp_LUT_stuff(SZA,VZA,RA, 2)
    file_name = out_file_base + '/Band3_L8_RTC_interpolated_values.nc'
    Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

    #print 'Processing band 4 ...\n'
    Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_L8.Interp_LUT_stuff(SZA,VZA,RA, 3)
    file_name = out_file_base + '/Band4_L8_RTC_interpolated_values.nc'
    Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)
  
    ## Need to generate higher res resampled rasters rasters for band8

    #VZA=open_geotiff(VZA_file_B8)
    #SZA=open_geotiff(SZA_file_B8)
    #VA=open_geotiff(VA_file_B8)
    #SA=open_geotiff(SA_file_B8) 

    #RA = abs(abs(abs(VA - SA)-180.0)-180.0)

    #print 'Processing band 8 ...\n'
    Lp_0, Eg_0, T_up,S= read_MODTRAN_lut_L8.Interp_LUT_stuff(SZA,VZA,RA, 7)
    file_name = out_file_base + '/Band8_L8_RTC_interpolated_values.nc'
    Write_RTC_NetCDF(file_name, Lp_0, Eg_0, T_up,S)

    #print 'RTC band8 shape ',VZA.shape


    return



