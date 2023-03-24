
import numpy as np 
import xarray as xr
import matplotlib.pyplot as plt


from libtiff import TIFF
import pyproj as proj
import warnings

import sys
sys.path.append("m:/job_access/python/resample_wts/AMSR2")  #contains the target footprint pattern code
from AMSR2_Antenna_Gain import target_gain


warnings.filterwarnings('ignore')
#turn off the warnings for proj


def calc_footprint_LF_from_nsidc_3km(*,lat,lon,land,grid_lats,grid_lons,diameter_in_km):

    print(np.max(land),np.min(land))

    #currently, the land data is integers coverted to float 32
    if np.max(land) < 0.0001:
        wted_land = 0.0
    elif np.min(land) > 0.9999:
        wted_land = 1.0
    else:
        # not all zero or all 1's -- need to do the integral

        crs_wgs = proj.Proj(init='epsg:4326')  # assuming you're using WGS84 geographic
        cust = proj.Proj(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m")

        #apply the local projection to obtain the mesh points in meters.
        xv, yv = proj.transform(crs_wgs, cust, grid_lons, grid_lats) 

        #convert to km
        yv = yv/1000.0
        xv = xv/1000.0

        dist_from_center = np.sqrt(np.square(xv) + np.square(yv))
        footprint_amplitude = target_gain(dist_from_center,diameter_in_km = diameter_in_km)
        
        #since these are from the EASE 2.0 grid, the cell sizes should be very close to the same,
        # so don't worry about cell size adjustments

        wted_land = np.sum(land*footprint_amplitude)/np.sum(footprint_amplitude)
    
    return wted_land

if __name__ == '__main__':

    from EASE2_grid_land_mask import read_EASE2_grid_land_mask
    from EASE2_grid import read_EASE2_grid_locations

    footprint_diameter_km = 70.0

    weight_map = np.zeros((721,1440))
    weight_map[0:41,:] = 1.0

    # if restart:
    #     nc_file = 'L:/access/land_water/land_fraction_1440_721_30km.nc'
    #     land_mask_DS=xr.open_dataset(nc_file)  
    #     weight_map = land_mask_DS['land_fraction'].values

    lats,lons = read_EASE2_grid_locations(resolution=3,pole='south')
    lats_radians = np.deg2rad(lats)
    lons_radians = np.deg2rad(lons)

    land_present = read_EASE2_grid_land_mask(resolution=3,pole='south',convert_ice=True)

    for latitude0 in np.arange(-80.0,-54.9,0.25):
        ilat = int(np.round((latitude0+90.0)/0.25))
        latitude0_radians = np.deg2rad(latitude0)
        for longitude0 in np.arange(0.0,360.0,0.25):
            ilon = int(np.round(longitude0/0.25))
            longitude0_radians = np.deg2rad(longitude0)

            dlat = lats_radians - latitude0_radians
            dlon = lons_radians - longitude0_radians

            #first guess based on delta_lat,delta_lon
            max_lat_diff = 1.0
            max_lon_diff = max_lat_diff/np.cos(latitude0_radians)
            ok1 = np.all([(np.abs(dlat) < np.deg2rad(max_lat_diff)),np.any([(np.abs(dlon) < np.deg2rad(max_lon_diff)),
                                                                    (np.abs(dlon) > np.deg2rad(360.0-max_lon_diff))],axis=0)],axis=0)
            
            #better guess based on haversine formula
            a = np.square(np.sin(dlat[ok1]/2.0)) + np.cos(latitude0_radians)*np.cos(lats_radians[ok1])*np.square(np.sin(dlon[ok1]/2.0))
            dist_approx = 6371*(2.0*np.arctan2(np.sqrt(a),np.sqrt(1.0-a)))
            ok2 = dist_approx <= 3.0*footprint_diameter_km

            # pass land_present, grid_lat and grid lon arrays selected by both ok1 and ok2

            wt = calc_footprint_LF_from_nsidc_3km(lat = latitude0,
                                                 lon =  longitude0,
                                                 land = land_present[ok1][ok2],
                                                 grid_lats = lats[ok1][ok2],
                                                 grid_lons = lons[ok1][ok2],
                                                 diameter_in_km = footprint_diameter_km)

            weight_map[ilat,ilon] = wt

            print(f'{np.sum(ok1):06d} {np.sum(ok2):05d} {ilat:03d} {ilon:04d} {latitude0:7f}, {longitude0:7f}: {wt:7f}')

        lon_array = np.arange(0.0,360.0,0.25)
        lat_array = np.arange(-90.0,90.1,0.25)
        land_mask_DS = xr.Dataset(
            data_vars = dict(
            land_fraction = (["Latitude","Longitude"],weight_map),
                            ),
                            coords = {
                                "Latitude":lat_array,
                                "Longitude":lon_array
                                }
                            )

        nc_file = f'L:/access/land_water/land_fraction_1440_721_{int(footprint_diameter_km)}km.from_nsidc_3km_mask.south.v3.nc'
        encoding = {"land_fraction":{'zlib' : True, 'complevel': 4 }}
        land_mask_DS.to_netcdf(nc_file,encoding=encoding)   
        print(f'Finished Lat = {latitude0}') 
        print(f'Wrote: {nc_file}')    


