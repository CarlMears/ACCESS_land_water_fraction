

import pytiff

file = 'C:/Users/carl/Dropbox/docs/ACCESS/PROBAV_LC100_global_v3.0.1_2015-base_SeasonalWater-CoverFraction-layer_EPSG-4326.tif'
with pytiff.Tiff(file) as handle:
    part = handle[0:10000.0,0:10000.0]

print()
