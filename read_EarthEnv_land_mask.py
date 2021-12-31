import numpy as np
import zarr
from geotiff import GeoTiff
import pyproj as proj
import warnings
  
import sys
sys.path.append("m:/job_access/python/rebinning")
sys.path.append("m:/job_access/python/projections")

from bin_ndarray import bin_ndarray
from gauss_kruger import local_geogauss

warnings.filterwarnings('ignore')


def read_EarthEnv_land_mask(*,filename = 'L:/access/land_water/EarthEnv/Consensus_reduced_class_12.open_water.tif',
                            center_lon,center_lat,delta_lon,delta_lat):

    geo_tiff = GeoTiff(filename)

    lower_lon = center_lon - delta_lon
    upper_lon = center_lon + delta_lon

    upper_lat = center_lat + delta_lat
    lower_lat = center_lat - delta_lat

    area_box = ((lower_lon,lower_lat),(upper_lon,upper_lat))
    
    lf = geo_tiff.read_box(area_box)
    num_lats,num_lons = lf.shape

    lat_arr_temp = lower_lat + (1.0/120.0)*np.arange(num_lats,0,-1)
    lat_arr = np.zeros((num_lats,num_lons))
    for i in range(0,num_lons):
        lat_arr[:,i] = lat_arr_temp

    lon_arr_temp = lower_lon + (1.0/120.0)*np.arange(0,num_lons)
    lon_arr = np.zeros((num_lats,num_lons))
    for i in range(0,num_lats):
        lon_arr[i,:] = lon_arr_temp

    return lf,lat_arr,lon_arr

def calc_gaussian_LF_EarthEnv_land_mask(*,lat,lon,gausian_radius):

    # Round the lat/lon to the nearest degree for loading.

    center_lat = int(np.round(lat))
    center_lon = int(np.round(lon))

    delta_lon = 4.0
    delta_lat = 4.0

    lf,lat_arr,lon_arr = read_EarthEnv_land_mask(filename = 'L:/access/land_water/EarthEnv/Consensus_reduced_class_12.open_water.tif',
                            center_lon=center_lon,center_lat=center_lat,delta_lon=delta_lon,delta_lat=delta_lat)
    
    crs_wgs = proj.Proj(init='epsg:4326')  # assuming you're using WGS84 geographic
    cust = proj.Proj(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m")

    #apply the local projection to obtain the mesh points in meters.
    xv, yv = proj.transform(crs_wgs, cust, lon_arr, lat_arr)   

    #convert to km
    yv = yv/1000.0
    xv = xv/1000.0

    #target size in km
    sigmax = 30.0
    sigmay = 30.0

    # calculate the footprint - don't car about normailization -- normalization is
    # done "by hand" later.

    gauss = np.exp(-((np.square(xv)/(2.0*np.square(sigmax))) + 
                   + (np.square(yv)/(2.0*np.square(sigmay)))))

    # calculate the adjustments for changes in grid size
    # The size of the lat/lon cells vary with latitude.
    # We'll use the information from the projected mesh to estimate this
    # with high accuracy.
2
    cell_width_temp  =  xv[:,2:-1] - xv[:,0:-2]
    cell_height_temp =  yv[0:-2,:] - yv[2:-1,:]

    cell_width = np.zeros((2*half_numy+1,2*half_numx+1))
    cell_height = np.zeros((2*half_numy+1,2*half_numx+1))

    cell_width[1:2*half_numy,1:2*half_numx] = cell_width_temp
    cell_height[1:2*half_numy,1:2*half_numx] = cell_height_temp

    cell_width[0,:] = cell_width[1,:]
    cell_width[2*half_numy,:] = cell_width[2*half_numy-1,:]
    cell_width[:,0] = cell_width[:,1]
    cell_width[:,2*half_numx] = cell_width[:,2*half_numx-1]

    cell_height[0,:] = cell_height[1,:]
    cell_height[2*half_numy,:] = cell_height[2*half_numy-1,:]
    cell_height[:,0] = cell_height[:,1]
    cell_height[:,2*half_numx] = cell_height[:,2*half_numx-1]

    cell_area = cell_height*cell_width
                              
    print()

if __name__ == '__main__':

    read_EarthEnv_land_mask()
    print()
