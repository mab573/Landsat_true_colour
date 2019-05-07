import os
import os.path
import click

from glob import glob

from wagl.acquisition import acquisitions

from calc_sat_solar import unpack, write_sat_solar
from calc_sat_solar import get_band_scale_offset, get_data
from calc_sat_solar import normalize_data
from calc_sat_solar import generate_rtc_raster

from Simple_Pan_Sharpen import Simple_Pan_Sharpen as pan
from Landsat8_atmospheric_correction import Landsat_ATCOR

import numpy as np

from PIL import Image


MAX_REFL = 10000


def gamma(ary,bright):
    ary_scaled=((ary / 255.0) ** (1.0 /bright))*255.0
    return ary_scaled

def png_band(band):
    return gamma(np.clip(band / np.max(band) * 255.0, 0, 255), 1.0)
    

def atcor(rtc_data, radiance):
    Lp_0 = rtc_data['Lp_0']
    Eg_0 = rtc_data['Eg_0']
    T_up = rtc_data['T_up']
    S = rtc_data['S']
    A = (np.pi*((radiance/10.0)-(Lp_0)))/((Eg_0)*T_up)
    return np.around(((A)/(1+(A*S)))*10000.0, decimals=0)


@click.command()
@click.option('--level1', type=click.Path(exists=True), required=True,
              help='location of level1 tar file')
@click.option('--outdir', type=click.Path(), required=True,
              help='output directory')
@click.option('--extent', type=(float, float, float, float), default=(None, None, None, None),
              help='extent to subset in UL-lat UL-lon LR-lat LR-lon format')
@click.option('--sharpen', is_flag=True, default=False,
              help='whether to pan sharpen')
@click.option('--ac', default=False)
@click.option('--cleanup', is_flag=True, default=False,
              help='whether to clean up working directory')
def main(level1, outdir, extent, sharpen, ac, cleanup):
    acqs = acquisitions(level1)
    assert len(acqs.granules) == 1, 'cannot handle multi-granule datasets'
    granule = acqs.granules[0]
    acqs = {group: acqs.get_acquisitions(group=group, granule=granule, only_supported_bands=False)
            for group in acqs.groups}

    extracted = os.path.join(outdir, 'extracted')

    if not os.path.exists(os.path.join(extracted, '{}_{}'.format(granule, 'SOLAR-ZENITH.TIF'))):
        print('extracting... ', end='')
        os.makedirs(extracted, exist_ok=True)

        unpack(level1, extracted)
        write_sat_solar(granule, acqs['RES-GROUP-1'][0], extracted)
        print('done!')
    else:
        print('found extracted files')

    scale_offset_dict = get_band_scale_offset(glob(os.path.join(extracted, '*MTL.txt'))[0])
    data = get_data(extracted, extent, sharpen)
    
    rtc_data = generate_rtc_raster(data['SATELLITE-VIEW']['data'], data['SOLAR-ZENITH']['data'],
                                   data['SATELLITE-AZIMUTH']['data'], data['SOLAR-AZIMUTH']['data'])

    radiance = normalize_data(data, scale_offset_dict)
    
    if ac: 
        rho_r = atcor(rtc_data['B4'], radiance['B4'])
        rho_g = atcor(rtc_data['B3'], radiance['B3'])
        rho_b = atcor(rtc_data['B2'], radiance['B2'])
        rho_p = atcor(rtc_data['B8'], radiance['B8'])

    if sharpen:
        visible_bands = [rho_b, rho_g, rho_r, rho_p]
        #visible_bands = pan(radiance['B2'], radiance['B3'], radiance['B4'], radiance['B8'])
    else: 
        visible_bands = [radiance['B2'], radiance['B3'], radiance['B4']]
        #visible_bands = [rho_b, rho_g, rho_r]
       
    png_bands = [png_band(band) for band in visible_bands]
    
    jr = Image.fromarray(png_bands[2].astype(np.uint8))
    jg = Image.fromarray(png_bands[1].astype(np.uint8))
    jb = Image.fromarray(png_bands[0].astype(np.uint8))


    imrgb = Image.merge('RGB', (jr, jg, jb ))
    imrgb.save(os.path.join(extracted, 'corrected_sharpen_rgb.png'))


if __name__ == '__main__':
    main()
