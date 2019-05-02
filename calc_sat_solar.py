#!/usr/bin/env python 

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


SAT_SOLAR_BANDS = ['SATELLITE-VIEW', 'SATELLITE-AZIMUTH', 'SOLAR-ZENITH', 'SOLAR-AZIMUTH']


def calc_lat_lon_grids(level1, out_dir, tle_path='/g/data/v10/eoancillarydata/sensor-specific', acq_parser_hint=''):
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


if __name__ == '__main__': 

    level1 = '/g/data/da82/AODH/USGS/L1/Landsat/C1/093_085/LC80930852019051/LC08_L1TP_093085_20190220_20190222_01_T1.tar'
    out_dir = '/g/data/u46/users/pd1813/Landsat_true_colour/test'
    calc_lat_lon_grids(level1, out_dir)

