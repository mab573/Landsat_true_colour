import os
import os.path
import click

from glob import glob

from wagl.acquisition import acquisitions

from calc_sat_solar import unpack, write_sat_solar
from calc_sat_solar import get_band_scale_offset, get_data
from calc_sat_solar import normalize_data
from calc_sat_solar import generate_rtc_raster


MAX_REFL = 10000

@click.command()
@click.option('--level1', type=click.Path(exists=True), required=True,
              help='location of level1 tar file')
@click.option('--outdir', type=click.Path(), required=True,
              help='output directory')
@click.option('--extent', type=(float, float, float, float), default=(None, None, None, None),
              help='extent to subset in UL-lat UL-lon LR-lat LR-lon format')
@click.option('--sharpen', is_flag=True, default=False,
              help='whether to pan sharpen')
@click.option('--cleanup', is_flag=True, default=False,
              help='whether to clean up working directory')
def main(level1, outdir, extent, sharpen, cleanup):
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

    rtc = generate_rtc_raster(data['SATELLITE-VIEW']['data'], data['SOLAR-ZENITH']['data'],
                              data['SATELLITE-AZIMUTH']['data'], data['SOLAR-AZIMUTH']['data'])

    radiance = normalize_data(data, scale_offset_dict)

    for key in radiance:
        print(key)
        print(radiance[key])
        print()




if __name__ == '__main__':
    main()
