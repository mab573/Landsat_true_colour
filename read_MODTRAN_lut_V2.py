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

## Import local modules
#sys.path.insert(0, '/flurry/home/mbroom/src/AHI_code/True_colour/')
#sys.path.insert(0, '/flurry/home/mbroom/src/AHI_code/True_colour/AHI_look_up_generation/')
#sys.path.insert(0, '/flurry/home/mbroom/src/AHI_code/True_colour/new_code/')

import three_D_interpolate_V2

def Interp_LUT_stuff(data_file,ancillary_file,band):

    print data_file
    print ancillary_file

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

    ## The first time this is run the interpolation should be run for T_up and S
    ## T_up will be constant over the globe and not change
    ## The Spherical albedo will change but only very marginally. Using a single number
    ## for this (estimated from the S interpolated field) is sufficient.

    #pickleFile_S='/short/er8/mab573/seven_sc/AHI-08_S.dat'
    #pickleFile_T_up = '/short/er8/mab573/seven_sc/AHI-08_T_up.dat'
    #pickleFile_Eg_0 = '/short/er8/mab573/seven_sc/AHI-08_Eg_0.dat'
    #pickleFile_Lp_0 = '/short/er8/mab573/seven_sc/AHI-08_Lp_0.dat'

    print "Definately reading the new LUT"

    pickleFile_S='/short/er8/mab573/seven_sc/AHI-08_TROP_S.dat'
    pickleFile_T_up = '/short/er8/mab573/seven_sc/AHI-08_TROP_T_up.dat'
    pickleFile_Eg_0 = '/short/er8/mab573/seven_sc/AHI-08_TROP_Eg_0.dat'
    pickleFile_Lp_0 = '/short/er8/mab573/seven_sc/AHI-08_TROP_Lp_0.dat'

    #pickleFile_T_up = '/short/er8/mab573/seven_sc/AHI-08_T_up_new.dat'
    #pickleFile_Eg_0 = '/short/er8/mab573/seven_sc/AHI-08_Eg_0_new.dat'
    #pickleFile_Lp_0 = '/short/er8/mab573/seven_sc/AHI-08_Lp_0_new.dat'

    solzen, satzen, raz, Lp_0 = getMODLut_Lp(pickleFile_Lp_0)
    solzen, satzen, raz, Eg_0 = getMODLut_Eg(pickleFile_Eg_0)

    solzen, satzen, raz, T_up = getMODLut_Tup(pickleFile_T_up)
    solzen, satzen, raz, S = getMODLut_S(pickleFile_S)

    solzen=np.asarray(solzen,dtype=int)
    satzen=np.asarray(satzen,dtype=int)
    raz=np.asarray(raz,dtype=int)

    ## Retrieve all of the important navigational data
    solar_az = getNav(data_file, 'solar_azimuth_angle')
    sol_zen = getNav(data_file, 'solar_zenith_angle')
    #lat = getNav(ancillary_file, 'lat')
    #lon = getNav(ancillary_file, 'lon')
    sensor_az = getNav(ancillary_file, 'sensor_azimuth_angle')
    sen_zen = getNav(ancillary_file, 'sensor_zenith_angle')

    '''
    dir ='/g/ns/oeb/scratch1/mbroom/data/AHI/netCDF/'
    for f in glob.glob(dir + "/*-P1S-ABOM_GEOM_SOLAR-PRJ_GEOS141_2000-*.nc"):
        sol_az = getNav(f, 'solar_azimuth_angle')
        sol_zen = getNav(f, 'solar_zenith_angle')

    dir ='/g/ns/oeb/scratch1/mbroom/data/AHI/ancillary/'
    for f in glob.glob(dir+"/*-P1S-ABOM_GEOM_SENSOR-PRJ_GEOS141_2000-*.nc"):

        lat = getNav(f, 'lat')
        lon = getNav(f, 'lon')
        sen_az = getNav(f, 'sensor_azimuth_angle')
        sen_zen = getNav(f, 'sensor_zenith_angle')
    '''

    sen_zen_um=sen_zen.filled(fill_value=90.0)
    #sen_az_um=sen_az.filled(fill_value=360.0)

    ## This should convert actual angle differences (up to 360) to relative difference (up to 180)
    #rel_az =abs(sen_az_um-sol_az)%180
    #rel_az = abs(sen_az_um%180 - sol_az%180)
    data=np.zeros((8),dtype=float)

    # To use the satellite az angles do
    # raz = (360+360-(180-senaz+solaz)) % 360
    # raz = raz > 180 ? 360-ra : raz

    # set angles higher or equal to the max value to just under the max value
    # This makes it work with the rest of the table entry selection code

    '''
    sol_az_1=sol_az%-180
    sol_az_1[sol_az < 180]=0.0
    sol_az_1=sol_az_1*-1.0
    sol_az_2=sol_az%180
    sol_az_2[sol_az >= 180] = 0.0
    sol_az=sol_az_2+sol_az_1

    sen_az_1 = sen_az % -180
    sen_az_1[sen_az < 180] = 0.0
    sen_az_1 = sen_az_1 * -1.0
    sen_az_2 = sen_az % 180
    sen_az_2[sen_az >= 180] = 0.0
    sen_az = sen_az_2 + sen_az_1

    rel_az = abs(sen_az - sol_az)

    del sen_az, sen_az_2, sen_az_1, sol_az, sol_az_2, sol_az_1
    '''

    #sol_az = abs(solar_az - 180.0)
    #sen_az = abs(sensor_az - 180.0)
    #rel_az = abs(sen_az - sol_az)

    #sol_az = abs(abs(solar_az - 180.0) - 180.0)
    #sen_az = abs(abs(sensor_az - 180.0) - 180.0)
    #rel_az = abs(sen_az - sol_az)

    rel_az = abs(abs(abs(sensor_az - solar_az)-180.0)-180.0)

    #rel_az = abs(abs(sensor_az - solar_az) - 180.0)
    
    filename = 'obs_geometry_20160126_0400.nc'
    nc = netCDF4.Dataset(filename, "w", format='NETCDF4')
    nc.createDimension("x_dim", len(rel_az[:, 0]))
    nc.createDimension("y_dim", len(rel_az[0, :]))
    v = nc.createVariable("solar_zen", float, ("x_dim", "y_dim"))
    v[:, :] = sol_zen
    v = nc.createVariable("sensor_zen", float, ("x_dim", "y_dim"))
    v[:, :] = sen_zen_um
    v = nc.createVariable("rel_az", float, ("x_dim", "y_dim"))
    v[:, :] = rel_az
    nc.close()

    sys.exit()
    

    del solar_az
    del sensor_az
    #del sol_az
    #del sen_az

    ## LUT entries only go to 89. This should prevent the calculated distance from becoming weird
    #sol_zen[sol_zen > 89.0]=89.0
    #sen_zen_um[sen_zen_um > 89.0] = 89.0

    sol_zen[sol_zen > 89.0] = 89.0
    sen_zen_um[sen_zen_um > 89.0] = 89.0

    max_raz = 180.0
    max_solzen = 89.0
    max_satzen = 89.0

    rel_az[rel_az >= max_raz] = max_raz - 0.001
    #sol_az[sol_az >= max_raz] = np.max(raz) - 0.001
    sol_zen[sol_zen >= max_solzen] = max_solzen - 0.001
    sen_zen_um[sen_zen_um >= max_satzen] = max_satzen - 0.001

    ## For each pixel, identify the nearest entries from the table

    low_sol_zen=np.zeros((len(sol_zen[:,0]),len(sol_zen[0,:])),dtype=int)
    low_relaz=np.zeros((len(sol_zen[:,0]),len(sol_zen[0,:])),dtype=int)
    low_sen_zen=np.zeros((len(sol_zen[:,0]),len(sol_zen[0,:])),dtype=int)

    Lp_interp=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:]),4),dtype=float)
    Eg_interp = np.zeros((len(sen_zen[:, 0]), len(sen_zen[0, :]), 4), dtype=float)
    #T_up_interp = np.zeros((len(sen_zen[:, 0]), len(sen_zen[0, :]), 4), dtype=float)
    #S_interp = np.zeros((len(sen_zen[:, 0]), len(sen_zen[0, :]), 4), dtype=float)


    dp_x=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:])),dtype=float)
    dp_y=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:])),dtype=float)
    dp_z=np.zeros((len(sen_zen[:,0]),len(sen_zen[0,:])),dtype=float)

    ## We know the step so just need to divide by step and floor it
    ## If the step changes then these will need to be changed

    low_sol_zen=np.floor(sol_zen/5)
    low_relaz=np.floor(rel_az/10)
    #low_relaz = np.floor(sol_az / 10)
    low_sen_zen=np.floor(sen_zen_um/5)

    low_sol_zen=low_sol_zen.astype(int)
    low_relaz=low_relaz.astype(int)
    low_sen_zen=low_sen_zen.astype(int)

    ## Calculate distances to the table entries for each point in the real-world data

    dp_x = (sol_zen - solzen[low_sol_zen])/5.0
    #dp_y = (sol_az - raz[low_relaz])/10.0
    dp_y = (rel_az - raz[low_relaz])/10.0
    dp_z = (sen_zen_um - satzen[low_sen_zen])/5.0

    #print 'Angles'
    #print sol_zen[2000, 3600], rel_az[2000, 3600], sen_zen_um[2000, 3600]
    #print dp_x[2000, 3600],dp_y[2000, 3600],dp_z[2000, 3600]
    #print sol_zen[1570, 3950], rel_az[1570, 3950], sen_zen_um[1570, 3950]
    #print dp_x[1570, 3950], dp_y[1570, 3950], dp_z[1570, 3950]

    del rel_az
    #del sol_az
    del sol_zen
    del sen_zen_um

    Lp_interp = three_D_interpolate_V2.Interp_3D(dp_x,dp_y,dp_z,low_sol_zen,low_relaz,low_sen_zen,Lp_0,band)
    Eg_interp = three_D_interpolate_V2.Interp_3D(dp_x, dp_y, dp_z, low_sol_zen, low_relaz, low_sen_zen, Eg_0, band)

    S_interp = three_D_interpolate_V2.Interp_3D(dp_x, dp_y, dp_z, low_sol_zen, low_relaz, low_sen_zen, S, band)
    T_up_interp = three_D_interpolate_V2.Interp_3D(dp_x, dp_y, dp_z, low_sol_zen, low_relaz, low_sen_zen, T_up, band)


    return Lp_interp, Eg_interp, T_up_interp,S_interp
    #return T_up_interp, S_interp
