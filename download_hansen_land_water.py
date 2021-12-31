import sys
sys.path.append("m:/job_access/python/rebinning")
sys.path.append("m:/job_access/python/projections")

import os
import numpy as np 
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from libtiff import TIFF
  
from bin_ndarray import bin_ndarray
from gauss_kruger import local_geogauss

import requests


for corner_lat in range(0,5,10):
    for corner_lon in range(0,351,10):

        if corner_lat >= 0:
            lat_string = f"{corner_lat:02d}N"
        else:
            lat_abs = abs(corner_lat)
            lat_string = f"{lat_abs:02d}S"

        if corner_lon >= 180:
            corner_lon = corner_lon - 360

        if corner_lon < 0:
            lon_string = f"{abs(corner_lon):03d}W"
        else:
            lon_string = f"{corner_lon:03d}E"

        url = f"https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_datamask_{lat_string}_{lon_string}.tif"
        outfile = f'L:/access/hansen_land_mask/Hansen_GFC2015_datamask_{lat_string}_{lon_string}.tif'
        ncfile = f'L:/access/hansen_land_mask/Hansen_GFC2015_datamask_{lat_string}_{lon_string}.regrid2.nc'
        print(f"Retrieving Land Mask from {url}")

        r = requests.get(url)
        with open(outfile,"wb") as f:
            f.write(r.content)

        print(f"Writing temporary file {outfile}")
        land_mask_fyl = TIFF.open(outfile)
        land_mask = land_mask_fyl.read_image()
        land_mask[land_mask == 2] = 0.0

        # reduce the size of the array to 1000x1000 resolution = 10.0/1000.0 = 0.01 degrees
        land_mask2 = bin_ndarray(land_mask, (1000,1000), operation='mean')
        #first index is latitude 
        #second index is longitude
        dlat_array = (0.5 + np.arange(0,1000))*(10.0/1000.0)
        dlon_array = (0.5 + np.arange(0,1000))*(10.0/1000.0)
        lat_array = corner_lat - dlat_array
        lon_array = corner_lon + dlon_array

        land_mask_DS = xr.Dataset(
                            data_vars = dict(
                            land_fraction = (["Latitude","Longitude"],land_mask2),
                            ),
                            coords = {
                                "Latitude":lat_array,
                                "Longitude":lon_array
                            }
                        )
        print(f"Saving regridded data to {ncfile}")
        encoding = {"land_fraction":{'zlib' : True, 'complevel': 9 }}
        
        land_mask_DS.to_netcdf(ncfile,encoding=encoding)
        #print(f"removing {outfile}")
        #os.remove(outfile)
        print()