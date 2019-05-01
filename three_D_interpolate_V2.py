#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

## This should interpolate for a space inside a cube. It is designed to produce values
## from Yi's lookup table, which is indexed by solar zenith, relative azimuth and view zenith
## This is based on trilinear interpolation methods from https://en.wikipedia.org/wiki/Trilinear_interpolation

def Interp_3D(dp_x,dp_y,dp_z,xx_0,yy_0,zz_0,ra,band):

    ## Hopefully this creates a bunch of arrays with the right data in it
    data_0 = ra[xx_0, yy_0, zz_0, band]
    data_1 = ra[xx_0+1, yy_0, zz_0, band]
    data_2 = ra[xx_0, yy_0, zz_0+1, band]
    data_3 = ra[xx_0+1, yy_0, zz_0+1, band]
    data_4 = ra[xx_0, yy_0+1, zz_0, band]
    data_5 = ra[xx_0+1, yy_0+1, zz_0, band]
    data_6 = ra[xx_0, yy_0+1, zz_0+1, band]
    data_7 = ra[xx_0+1, yy_0+1, zz_0+1, band]

    #print xx_0.shape, data_0.shape

    #print 'Lp data at point'
    #print xx_0[1570, 3950], yy_0[1570, 3950], zz_0[1570, 3950]
    #print data_0[1570, 3950],data_1[1570, 3950],data_2[1570, 3950],data_3[1570, 3950],data_4[1570, 3950],data_5[1570, 3950],data_6[1570, 3950],data_7[1570, 3950]

    ## Calculate values in the x-direction
    c_00 = data_0*(1 - dp_x) + data_1*dp_x
    c_01 = data_2*(1 - dp_x) + data_3*dp_x
    c_10 = data_4*(1 - dp_x) + data_5*dp_x
    c_11 = data_6*(1 - dp_x) + data_7*dp_x

    #print c_00[1570, 3950],c_01[1570, 3950],c_10[1570, 3950],c_11[1570, 3950]

    ## Calculate the values in the y-direction
    c_0 = c_00*(1 - dp_y) + c_10*dp_y
    c_1 = c_01*(1 - dp_y) + c_11*dp_y

    #print c_0[1570, 3950],c_1[1570, 3950]

    ## Calculate value in the z-direction
    c = c_0*(1 - dp_z) + c_1*dp_z

    #print c[3950, 1570]


    return c
'''
def Interp_3D(x_point,y_point,z_point,dp_x,dp_y,dp_z,data):

    ## Calculate the distance between the point and the smallest nearby table entry
    x_d = (x_point - dp_x[0]) / (dp_x[1] - dp_x[0])
    y_d = (y_point - dp_y[0]) / (dp_y[1] - dp_y[0])
    z_d = (z_point - dp_z[0]) / (dp_z[1] - dp_z[0])

    ## Calculate values in the x-direction
    c_00 = data[0]*(1 - x_d) + data[1]*x_d
    c_01 = data[2]*(1 - x_d) + data[3]*x_d
    c_10 = data[4]*(1 - x_d) + data[5]*x_d
    c_11 = data[6]*(1 - x_d) + data[7]*x_d


    ## Calculate the values in the y-direction
    c_0 = c_00*(1 - y_d) + c_10*y_d
    c_1 = c_01*(1 - y_d) + c_11*y_d

    ## Calculate value in the z-direction
    c = c_0*(1 - z_d) + c_1*z_d

    return c
'''

