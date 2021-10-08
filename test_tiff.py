import sys
sys.path.append("c:/job_access/python/rebinning")
sys.path.append("c:/job_access/python/projections")

import numpy as np 
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from libtiff import TIFF
import pyproj as proj
  
from bin_ndarray import bin_ndarray
from gauss_kruger import local_geogauss

# def bin_ndarray(ndarray, new_shape, operation='sum'):
#     """
#     Bins an ndarray in all axes based on the target shape, by summing or
#         averaging.

#     Number of output dimensions must match number of input dimensions and 
#         new axes must divide old ones.

#     Example
#     -------
#     >>> m = np.arange(0,100,1).reshape((10,10))
#     >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
#     >>> print(n)

#     [[ 22  30  38  46  54]
#      [102 110 118 126 134]
#      [182 190 198 206 214]
#      [262 270 278 286 294]
#      [342 350 358 366 374]]

#     """
#     operation = operation.lower()
#     if not operation in ['sum', 'mean']:
#         raise ValueError("Operation not supported.")
#     if ndarray.ndim != len(new_shape):
#         raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
#                                                            new_shape))
#     compression_pairs = [(d, c//d) for d,c in zip(new_shape,
#                                                   ndarray.shape)]
#     flattened = [l for p in compression_pairs for l in p]
#     ndarray = ndarray.reshape(flattened)
#     for i in range(len(new_shape)):
#         op = getattr(ndarray, operation)
#         ndarray = op(-1*(i+1))
#     return ndarray

corner_lat = 40
corner_lon = 360-130

if corner_lat > 0:
    lat_string = f"{corner_lat:02d}N"
else:
    if corner_lat == 0:
        lat_string = "0"
    else:
        lat_abs = abs(corner_lat)
        lat_string = f"{lat_abs:02d}S"

if corner_lon > 180:
    corner_lon = corner_lon - 360

if corner_lon < 0:
    lon_string = f"{abs(corner_lon):03d}W"
else:
    lon_string = f"{corner_lon:03d}E"

nc_land_fraction_file = f'L:/access/hansen_land_mask/Hansen_GFC2015_datamask_{lat_string}_{lon_string}.regrid.nc'
print(nc_land_fraction_file)
  
ds = xr.open_dataset(nc_land_fraction_file)

land_frac = ds.land_fraction.values
land_frac_lons = ds.Longitude.values
land_frac_lats = ds.Latitude.values

print(land_frac.shape)

#second index is longitude

latitude0 = 38.0



for longitude0 in [-122.0,-122.25,-122.50]:

    ilat0 = np.floor((corner_lat - latitude0)*400).astype(np.int32)
    ilat_test = land_frac_lats[ilat0-200:ilat0+201]
    ilon0 = np.floor((longitude0 - corner_lon)*400).astype(np.int32)
    ilon_test = land_frac_lons[ilon0-200:ilon0+201]

    submask = land_frac[ilat0-200:ilat0+201,ilon0-200:ilon0+201]
    plt.imshow(submask)


    #construct a local projection, centered at the center of the target gaussian

    crs_wgs = proj.Proj(init='epsg:4326')  # assuming you're using WGS84 geographic
    cust = proj.Proj(f"+proj=aeqd +lat_0={latitude0} +lon_0={longitude0} +datum=WGS84 +units=m")


    lats = latitude0 - np.arange(-200,201)*10.0/4000.0
    lons = longitude0 + np.arange(-200,201)*10.0/4000.0
    print(ilat0,ilon0)

    latv,lonv = np.meshgrid(ilat_test,ilon_test,indexing='ij')
    xv, yv = proj.transform(crs_wgs, cust, lonv, latv)   

    #xv,yv = local_geogauss(latitude0,longitude0,latv,lonv)
    yv = yv/1000.0
    xv = xv/1000.0

    sigmax = 20.0
    sigmay = 20.0
    gauss = np.exp(-((np.square(xv)/(2.0*np.square(sigmax))) + 
                   + (np.square(yv)/(2.0*np.square(sigmay)))))

    cell_width_temp  =  xv[1:400,2:401] - xv[1:400,0:399]
    cell_height_temp =  yv[0:399,1:400] - yv[2:401,1:400]

    cell_width = np.zeros((401,401))
    cell_height = np.zeros((401,401))

    cell_width[1:400,1:400] = cell_width_temp
    cell_height[1:400,1:400] = cell_height_temp

    cell_width[0,:] = cell_width[1,:]
    cell_width[400,:] = cell_width[399,:]
    cell_width[:,0] = cell_width[:,1]
    cell_width[:,400] = cell_width[:,399]

    cell_height[0,:] = cell_height[1,:]
    cell_height[400,:] = cell_height[399,:]
    cell_height[:,0] = cell_height[:,1]
    cell_height[:,400] = cell_height[:,399]

    cell_area = cell_height*cell_width

    wted_tot = np.sum(submask*gauss*cell_area)/np.sum(gauss*cell_area)
    print(f"Gaussing Weighted Land Fraction = {wted_tot}")
    fig = plt.imshow(submask*gauss*cell_area)
    png_file = f'C:/job_access/python/land_water/plots/wted_land_frac_{latitude0:06.2f}_{longitude0:06.2f}.png'
    plt.savefig(png_file)
plt.show()
print()




plt.show()




print()