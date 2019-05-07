import click

from wagl.acquisition import acquisitions


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
    for key in acqs:
        print(key)
        print(acqs[key])
        example = acqs[key][0]
        print(example)
        print(example.data())

if __name__ == '__main__':
    main()
