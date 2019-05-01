#!/usr/bin/env python

## Information 24/06/16 ##
## This program has been modified to no longer do interpolations for transmittance and spherical albedo
## The transmittance only has to be done once and the spherical albedo now uses a single band equivalent value
## as the actual value over all SZA changes only marginally and there is no discernable effect from using a single
## value.

import sys
import glob
import pickle
import os.path
import argparse
import datetime
import netCDF4
import numpy as np
import numpy.ma as ma
import netCDF4
from scipy import interpolate
import matplotlib.pyplot as plt
from osgeo import gdal
from osgeo import osr


import three_D_interpolate_V2

#def open_tif(tif_name)
#    #Open and save image of green band
#    ds = gdal.Open(tif_name, gdal.GA_ReadOnly)
#    (X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
#    srs_wkt = ds.GetProjection()
#    Nx = ds.RasterXSize
#    Ny = ds.RasterYSize
#    ary = []
#    ary = ds.GetRasterBand(1).ReadAsArray()
#    return ary

def Interp_LUT_stuff(sol_zen, sen_zen, rel_az, band):

    def getMODLut_Lp(pickleFile):
        """
        Read the MODTRAN Lp LUT
        """
        with open(pickleFile, "r") as f:
            d = pickle.load(f)

        ## Pickled these as string lists which was a bit of a mistake.
        solin = d['sza']
        satin = d['vza']
        rain = d['raz']

        solzen = [float(i) for i in solin]
        satzen = [float(i) for i in satin]
        raz = [float(i) for i in rain]
        Lp_0 = d['Lp_0']
        return solzen, satzen, raz, Lp_0

    def getMODLut_Eg(pickleFile):
        """
        Read the MODTRAN Lp LUT
        """
        with open(pickleFile, "r") as f:
            d = pickle.load(f)

        ## Pickled these as string lists which was a bit of a mistake.
        solin = d['sza']
        satin = d['vza']
        rain = d['raz']

        solzen = [float(i) for i in solin]
        satzen = [float(i) for i in satin]
        raz = [float(i) for i in rain]
        Eg_0 = d['Eg_0']
        return solzen, satzen, raz, Eg_0

    def getMODLut_Tup(pickleFile):
        """
        Read the MODTRAN Lp LUT
        """
        with open(pickleFile, "r") as f:
            d = pickle.load(f)

        ## Pickled these as string lists which was a bit of a mistake.
        solin = d['sza']
        satin = d['vza']
        rain = d['raz']

        solzen = [float(i) for i in solin]
        satzen = [float(i) for i in satin]
        raz = [float(i) for i in rain]
        T_up = d['T_up']
        return solzen, satzen, raz, T_up

    def getMODLut_S(pickleFile):
        """
        Read the MODTRAN Lp LUT
        """
        with open(pickleFile, "r") as f:
            d = pickle.load(f)

        ## Pickled these as string lists which was a bit of a mistake.
        solin = d['sza']
        satin = d['vza']
        rain = d['raz']

        solzen = [float(i) for i in solin]
        satzen = [float(i) for i in satin]
        raz = [float(i) for i in rain]
        S = d['S']    #made a mistake when I pickled this so the name isn't right but it doesn't matter
        return solzen, satzen, raz, S

    def getNav(fileName, varName):

        nc = netCDF4.Dataset(fileName)
        if varName not in nc.variables:
            sys.exit("Var "+varName+" not in file "+fileName)

        v = nc.variables[varName][0,:,:]
        nc.close()
        return v

    ## LUT with Landsat8 OLI specific bands (all of them)
    LUT_path='/home/573/mab573/Landsat/lookup_tables_code/'
    pickleFile_S=LUT_path+'Landsat-08_S.dat'
    pickleFile_T_up =LUT_path+'Landsat-08_T_up.dat'
    pickleFile_Eg_0 =LUT_path+'Landsat-08_Eg_0.dat'
    pickleFile_Lp_0 =LUT_path+'Landsat-08_Lp_0.dat'

    solzen, satzen, raz, Lp_0 = getMODLut_Lp(pickleFile_Lp_0)
    solzen, satzen, raz, Eg_0 = getMODLut_Eg(pickleFile_Eg_0)

    solzen, satzen, raz, T_up = getMODLut_Tup(pickleFile_T_up)
    solzen, satzen, raz, S = getMODLut_S(pickleFile_S)

    solzen=np.asarray(solzen,dtype=int)
    satzen=np.asarray(satzen,dtype=int)
    raz=np.asarray(raz,dtype=int)


    ### These are all tif files
    #sol_zen =open_tif(solar_zen) 
    #rel_az = open_tif(rel_az)
    #sen_ele = open_tif(sat_el)
    #sen_zen=90.-sen_ele

    ## Check to see how rel_az is defined in my LUT and with GA L8
    ## I may need to apply a correction
    #rel_az = abs(abs(abs(sensor_az - solar_az)-180.0)-180.0)

    sol_zen[sol_zen > 89.0] = 89.0
    sen_zen[sen_zen > 89.0] = 89.0

    max_raz = 180.0
    max_solzen = 89.0
    max_satzen = 89.0

    rel_az[rel_az >= max_raz] = max_raz - 0.001
    sol_zen[sol_zen >= max_solzen] = max_solzen - 0.001
    sen_zen[sen_zen >= max_satzen] = max_satzen - 0.001

    ## For each pixel, identify the nearest entries from the table

    low_sol_zen=np.zeros((len(sol_zen[:,0]),len(sol_zen[0,:])),dtype=int)
    low_relaz=np.zeros((len(sol_zen[:,0]),len(sol_zen[0,:])),dtype=int)
    low_sen_zen=np.zeros((len(sol_zen[:,0]),len(sol_zen[0,:])),dtype=int)

    Lp_interp=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:]),4),dtype=float)
    Eg_interp = np.zeros((len(sen_zen[:, 0]), len(sen_zen[0, :]), 4), dtype=float)
    T_up_interp = np.zeros((len(sen_zen[:, 0]), len(sen_zen[0, :]), 4), dtype=float)
    S_interp = np.zeros((len(sen_zen[:, 0]), len(sen_zen[0, :]), 4), dtype=float)


    dp_x=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:])),dtype=float)
    dp_y=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:])),dtype=float)
    dp_z=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:])),dtype=float)

    ## We know the step so just need to divide by step and floor it
    ## If the step changes then these will need to be changed

    low_sol_zen=np.floor(sol_zen/5)
    low_relaz=np.floor(rel_az/10)
    low_sen_zen=np.floor(sen_zen/5)

    low_sol_zen=low_sol_zen.astype(int)
    low_relaz=low_relaz.astype(int)
    low_sen_zen=low_sen_zen.astype(int)

    ## Calculate distances to the table entries for each point in the real-world data

    dp_x = (sol_zen - solzen[low_sol_zen])/5.0
    dp_y = (rel_az - raz[low_relaz])/10.0
    dp_z = (sen_zen - satzen[low_sen_zen])/5.0

    del rel_az
    del sol_zen
    del sen_zen

    Lp_interp = three_D_interpolate_V2.Interp_3D(dp_x,dp_y,dp_z,low_sol_zen,low_relaz,low_sen_zen,Lp_0,band)
    Eg_interp = three_D_interpolate_V2.Interp_3D(dp_x, dp_y, dp_z, low_sol_zen, low_relaz, low_sen_zen, Eg_0, band)

    S_interp = three_D_interpolate_V2.Interp_3D(dp_x, dp_y, dp_z, low_sol_zen, low_relaz, low_sen_zen, S, band)
    T_up_interp = three_D_interpolate_V2.Interp_3D(dp_x, dp_y, dp_z, low_sol_zen, low_relaz, low_sen_zen, T_up, band)

    return Lp_interp, Eg_interp, T_up_interp,S_interp
