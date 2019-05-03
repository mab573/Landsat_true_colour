#!/usr/bin/env python 

import numpy as np

import rasterio 
from rasterio.warp import transform_bounds, reproject, Resampling
from rasterio.crs import CRS
from affine import Affine

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

_log = logging.getLogger(__name__)


SAT_SOLAR_BANDS = ['SATELLITE-VIEW', 'SATELLITE-AZIMUTH', 'SOLAR-ZENITH', 'SOLAR-AZIMUTH']
LANDSAT_BANDS = ['B{}'.format(i) for i in range(1,12)]


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


def write_sat_solar(level1, out_dir, tle_path='/g/data/v10/eoancillarydata/sensor-specific', acq_parser_hint=''):
    """
    This function will compute satellite solar (view and azimuth) angles using level 1 dataset.
    """
    container = acquisitions(level1, acq_parser_hint)
 
    acq = (container.get_acquisitions(container.supported_groups[0], container.granules[0]))[0]
    
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
                out_name = pjoin(out_dir, '{}_{}.TIF'.format(container.granules[0], band))
                write_img(dset, out_name, geobox=geo_box, nodata=dset.attrs.get('no_data_value'))


def subset_img(filename, subset_coords=None, down_sample_factor=None):
    with rasterio.open(filename) as raster:
        if subset_coords is None:
            data = raster.read(1)
            affine = raster.transform

        else:
            ul_lat, ul_lon, lr_lat, lr_lon = subset_coords

            wgs = CRS.from_epsg(4326)

            left, bottom, right, top = transform_bounds(wgs, raster.crs,
                                                        ul_lon, lr_lat, lr_lon, ul_lat)

            start = ~raster.transform * (left, top)
            finish = ~raster.transform * (right, bottom)

            start = [int(x) for x in reversed(start)]
            finish = [int(x) for x in reversed(finish)]

            print('raster shape: ', raster.shape)
            print('start before clipping', start)
            print('finish before clipping', finish)
            if start[0] < 0:
                start[0] = 0
            if start[1] < 0:
                start[1] = 0
            if finish[0] > raster.shape[0]:
                finish[0] = raster.shape[0]
            if finish[1] > raster.shape[1]:
                finish[1] = raster.shape[1]

            print('start after clipping', start)
            print('finish after clipping', finish)
            window = ((start[0], finish[0]), (start[1], finish[1]))

            data = raster.read(1, window=window)
            affine = rasterio.windows.transform(window, raster.transform)

        if down_sample_factor is None:
            return {
                'crs': raster.crs,
                'transform': affine,
                'data': data
            }

        factor = down_sample_factor
        new_data = np.empty(shape=(round(data.shape[0] * factor),
                                   round(data.shape[1] * factor)),
                            dtype=data.dtype)

        new_affine = Affine(affine.a / factor, affine.b, affine.c,
                            affine.d, affine.e / factor, affine.f)

        reproject(data, new_data, src_transform=affine, dst_transform=new_affine,
                  src_crs=raster.crs, dst_crs=raster.crs, resampling=Resampling.nearest)

        return {
            'crs': raster.crs,
            'transform': new_affine,
            'data': new_data
        }





def main(leve1, out_dir, subset_coords=None, Bands=None): 
    # write_sat_solar(leve1, out_dir)
    
    # unpack(level1, out_dir)
    
    if Bands: 
        bands = Bands 
    else: 
        bands = LANDSAT_BANDS
    
    tiff_files = [f for f in os.listdir(out_dir)]

    for band in bands + SAT_SOLAR_BANDS: 
        try:
            print(band)
            tiff_file = fnmatch.filter(tiff_files, '*{}.TIF'.format(band))[0]
            if band is not 'B8':
                down_sampling_factor = 2 
            else: 
                down_sampling_factor = None

            dicts = subset_img(pjoin(out_dir, tiff_file), subset_coords, down_sampling_factor)
            print(dicts['transform'])
            print(dicts['data'].shape)
        except IndexError: 
            pass 
        



if __name__ == '__main__':
    Bands = ['B2', 'B3', 'B4', 'B8']
    subset_coords = (-35.5, 144.0, -37.0, 145.0)
    level1 = '/g/data/da82/AODH/USGS/L1/Landsat/C1/093_085/LC80930852019051/LC08_L1TP_093085_20190220_20190222_01_T1.tar'
    out_dir = '/g/data/u46/users/pd1813/Landsat_true_colour/test'
    
    main(level1, out_dir, subset_coords, Bands)






