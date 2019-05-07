#!/usr/bin/env python 

import numpy as np

import rasterio 
from rasterio.warp import transform_bounds, reproject, Resampling
from rasterio.crs import CRS
from affine import Affine
from rasterio import shutil as rio_shutil 
from rasterio.vrt import WarpedVRT
import os
from os.path import join as pjoin, basename, dirname 
import h5py
import tempfile 
from wagl.hdf5 import H5CompressionFilter
from wagl.longitude_latitude_arrays import _create_lon_lat_grids
from wagl.satellite_solar_angles import calculate_angles 
from wagl.acquisition import acquisitions
from wagl.geobox import GriddedGeoBox 
from wagl.data import write_img 
import pyproj
from osgeo import gdal 
from osgeo import gdalconst
import fnmatch
import tarfile 
import logging 

from memory_profiler import profile 
import time 

import read_MODTRAN_lut_L8 

_log = logging.getLogger(__name__)


SAT_SOLAR_BANDS = ['SATELLITE-VIEW', 'SATELLITE-AZIMUTH', 'SOLAR-ZENITH', 'SOLAR-AZIMUTH']
LANDSAT_BANDS = ['B{}'.format(i) for i in range(2,5)]


def get_band_scale_offset(mtl_file): 
    """
    Returns the scale and offset factors from mtl file for reflectance 
    and radiance for each bands. 
    """

    rad_scale_tags=['RADIANCE_MULT_BAND_{}'.format(i) for i in range(1, 12)]
    rad_offset_tags = ['RADIANCE_ADD_BAND_{}'.format(i) for i in range(1, 12)] 
    ref_scale_tags = ['REFLECTANCE_MULT_BAND_{}'.format(i) for i in range(1, 10)]
    ref_offset_tags = ['REFLECTANCE_ADD_BAND_{}'.format(i) for i in range(1, 10)]
    
    tags = rad_scale_tags + rad_offset_tags + ref_scale_tags + ref_offset_tags

    with open(mtl_file, 'r') as fid:
        meta_file_lines = fid.readlines()
    
    meta_file_dict = {line.split('=')[0].strip(): float(line.split('=')[1].strip().rstrip())
                      for line in meta_file_lines for tag in tags  if fnmatch.fnmatch(line, '*{}*'.format(tag))}
   
    return meta_file_dict   


def unpack(tar_file, out_dir): 
    """
    Unpacks the tar or tar gz  file 
    """
    if tar_file.endswith("tar.gz"): 
        
        with tarfile.open(tar_file, "r:gz") as tar:
            if out_dir:
                return tar.extractall(out_dir)
            return tar.extractall()

    with  tarfile.open(tar_file, "r:") as tar:
        if out_dir:
            return tar.extractall(out_dir)
        return tar.extractall()


def write_sat_solar(granule, acq, out_dir,
                    tle_path='/g/data/v10/eoancillarydata/sensor-specific', acq_parser_hint=''):
    """
    This function will compute satellite solar (view and azimuth) angles using level 1 dataset.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        lat_lon_h5 = pjoin(tmp_dir, 'latitude_longitude.h5')
        sat_sol_h5 = pjoin(tmp_dir, 'satellite_solar.h5')
        _create_lon_lat_grids(acq, lat_lon_h5, H5CompressionFilter.LZF, {}) 
        
        with h5py.File(lat_lon_h5, 'r') as lon_lat_fid: 
            with h5py.File(sat_sol_h5, 'w') as sat_sol_fid: 
                lon_lat_grp = lon_lat_fid['LONGITUDE-LATITUDE']
                calculate_angles(acq, lon_lat_grp, sat_sol_fid, H5CompressionFilter.LZF, {}, 
                                 tle_path)
        
        with h5py.File(sat_sol_h5, 'r') as fid: 
            dset_group = fid['SATELLITE-SOLAR']
            for band in SAT_SOLAR_BANDS: 
                dset = dset_group[band]
                geo_box = GriddedGeoBox.from_h5_dataset(dset)
                out_name = pjoin(out_dir, '{}_{}.TIF'.format(granule, band))
                write_img(dset, out_name, geobox=geo_box, options={}, nodata=dset.attrs.get('no_data_value'))


def get_grid_options(filename, subset_coords=None):
    with rasterio.open(filename) as raster:
        if subset_coords is None:
            affine = raster.transform
            return {'resampling': Resampling.nearest,
                    'transform': raster.transform,
                    'crs': raster.crs,
                    'height': raster.shape[1],
                    'width': raster.shape[0], 
                    'nodata': raster.nodata
                   }
        else:
            ul_lat, ul_lon, lr_lat, lr_lon = subset_coords
            wgs = CRS.from_epsg(4326)
            left, bottom, right, top = transform_bounds(wgs, raster.crs,
                                                        ul_lon, lr_lat, lr_lon, ul_lat)
            
            start = ~raster.transform * (left, top)
            finish = ~raster.transform * (right, bottom)

            start = [int(x) for x in reversed(start)]
            finish = [int(x) for x in reversed(finish)]

            if start[0] < 0:
                start[0] = 0
            if start[1] < 0:
                start[1] = 0
            if finish[0] > raster.shape[0]:
                finish[0] = raster.shape[0]
            if finish[1] > raster.shape[1]:
                finish[1] = raster.shape[1]

            window = ((start[0], finish[0]), (start[1], finish[1]))

            data = raster.read(1, window=window)
            affine = rasterio.windows.transform(window, raster.transform)
          
            return {'resampling': Resampling.nearest,
                    'transform': affine,
                    'crs': raster.crs, 
                    'height': data.shape[1], 
                    'width': data.shape[0], 
   'nodata': raster.nodata
                    }


def subset(file_path, grid_specs): 
    with rasterio.open(file_path) as src: 
        with WarpedVRT(src, **grid_specs) as vrt:
            data = vrt.read(1)
            return {'data': data, 'profile': vrt.profile, 'shape': data.shape}


def generate_rtc_raster(vza, sza, va, sa):

    ra = abs(abs(abs(va - sa) - 180.0) - 180.0)
    
    lut_dict = {}

    for band in [2, 3, 4, 8]: 
        lp_0, eg_0, t_up, s = read_MODTRAN_lut_L8.Interp_LUT_stuff(sza, vza, ra, band - 1)
        lut_dict['B{}'.format(band)] = {}
        lut_dict['B{}'.format(band)]['Lp_0'] = lp_0 
        lut_dict['B{}'.format(band)]['Eg_0'] = eg_0 
        lut_dict['B{}'.format(band)]['T_up'] = t_up
        lut_dict['B{}'.format(band)]['S'] = s
    
    return lut_dict

    
def get_data(out_dir, subset_coords=None, upsample=True): 
    
    tiff_files = [f for f in os.listdir(out_dir) if f.endswith('.TIF')]
    bands = LANDSAT_BANDS + SAT_SOLAR_BANDS
        
    if upsample:
        bands.append('B8')
        grid_spec_file = pjoin(out_dir, fnmatch.filter(tiff_files, '*B8.TIF')[0])
    else: 
        grid_spec_file = pjoin(out_dir, fnmatch.filter(tiff_files, '*{}.TIF'.format(bands[0]))[0])

    grid_specs = get_grid_options(grid_spec_file, subset_coords)
    file_paths = [pjoin(out_dir, fnmatch.filter(tiff_files, '*{}.TIF'.format(band))[0]) for band in bands]
    
    return {bands[idx]: subset(file_path, grid_specs) for idx, file_path in  enumerate(file_paths)}


def normalize_data(data, scale_offset_dict):
    return {key: data[key]['data'] * scale_offset_dict['RADIANCE_MULT_BAND_{}'.format(key[1:])] + 
                              scale_offset_dict['RADIANCE_ADD_BAND_{}'.format(key[1:])]
            for key in data
            if key.startswith('B')}


def main(level1, out_dir, subset_coords=None, upsample=True): 

    #write_sat_solar(level1, out_dir)
    #unpack(level1, out_dir)
    
    files = [f for f in os.listdir(out_dir)]
    mtl_file = pjoin(out_dir, fnmatch.filter(files, '*MTL.txt')[0])
    scale_offset_dict = get_band_scale_offset(mtl_file)
    print(scale_offset_dict)
        
    all_data = get_data(out_dir, subset_coords, upsample)
    #generate_rtc_raster(all_data['SATELLITE-VIEW']['data'][0], all_data['SOLAR-ZENITH']['data'][0],
    #                    all_data['SATELLITE-AZIMUTH']['data'][0], all_data['SOLAR-AZIMUTH']['data'][0])
    landsat_band_data = {key: all_data[key]['data'][0] * scale_offset_dict['RADIANCE_MULT_BAND_{}'.format(key[1:])] + 
                              scale_offset_dict['RADIANCE_ADD_BAND_{}'.format(key[1:])] for key in all_data.keys() if fnmatch.fnmatch(key, '*B*')}

    

if __name__ == '__main__':
    upsample = True
    subset_coords = (-36.5, 144.5, -37.0, 145.0)
    level1 = '//data/da82/AODH/USGS/L1/Landsat/C1/093_085/LC80930852019051/LC08_L1TP_093085_20190220_20190222_01_T1.tar'
    out_dir = '/g/data/u46/users/pd1813/Landsat_true_colour/test'
    start_time = time.clock()
    print(time.clock() - start_time)
    main(level1, out_dir, subset_coords, upsample)
    #for key in data.keys():
    #    print(data[key]['data'].shape)




