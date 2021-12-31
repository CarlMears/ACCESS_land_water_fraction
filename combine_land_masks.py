
import numpy as np 
import xarray as xr
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import os
from rss_plotting.global_map import plot_global_map

land_path = Path('L:/access/land_water')
#read in map from Hansen Data
hansen_file = land_path / 'land_fraction_1440_721_30km.v3.nc'
land_hansen_xr = xr.open_dataset(hansen_file)
land_hansen = land_hansen_xr['land_fraction'].values
bad = ~np.isfinite(land_hansen)
land_hansen[bad] = 0.0

#read in NP NSIDC data
nsidc_NP_file = land_path / 'land_fraction_1440_721_30km.from_nsidc_3km_mask.north.v4.nc'
land_nsidc_NP_xr = xr.open_dataset(nsidc_NP_file)
land_nsidc_NP = land_nsidc_NP_xr['land_fraction'].values

#read in SP NSIDC data
nsidc_SP_file = land_path / 'land_fraction_1440_721_30km.from_nsidc_3km_mask.south.v3.nc'
land_nsidc_SP_xr = xr.open_dataset(nsidc_SP_file)
land_nsidc_SP = land_nsidc_SP_xr['land_fraction'].values

mask = np.zeros((721,1440))

mask[621:721,:] = 1.0
for ilat in range(610,621):
    mask[ilat,:] = (ilat-610)/10.0
mask[590:721,1220:1300] = 1.0
mask[590:634,1335:1400] = 0.0

mask[0:120,:] = 1.0
mask[85:120,1280:1320] = 0.0

plot_global_map(mask,vmin=0.0,vmax=1.0)

lats = land_hansen_xr['Latitude'].values
lons = land_hansen_xr['Longitude'].values

land_nsidc_combined = np.copy(land_nsidc_NP)
land_nsidc_combined[0:360,:] = land_nsidc_SP[0:360,:]

plot_global_map(land_nsidc_combined,vmin=0.0,vmax=1.0)
bad = ~np.isfinite(land_nsidc_combined)
land_nsidc_combined[bad] = 1.0
plot_global_map(land_nsidc_combined,vmin=0.0,vmax=1.0)


land_combined = mask*land_nsidc_combined + (1.0-mask)*land_hansen
plot_global_map(land_combined,vmin=0.0,vmax=1.0)

lon_array = np.arange(0.0,360.0,0.25)
lat_array = np.arange(-90.0,90.1,0.25)
land_combined_xr = xr.Dataset(
    data_vars = dict(
    land_fraction = (["Latitude","Longitude"],land_combined),
                    ),
                    coords = {
                        "Latitude":lat_array,
                        "Longitude":lon_array
                        }
                    )
plot_global_map(land_combined_xr['land_fraction'],vmin=0.0,vmax=1.0)
plt.show()
output_file = land_path / 'land_fraction_1440_721_30km.combined_hansen_nsidc.nc'
land_combined_xr.to_netcdf(output_file)

print()




