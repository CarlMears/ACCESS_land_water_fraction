
import numpy as np
import os 
import xarray as xr
import matplotlib.pyplot as plt


#from libtiff import TIFF
import pyproj as proj
from pyproj import Geod
import warnings

import sys
sys.path.append("m:/job_access/python/resample_wts/AMSR2")  #contains the target footprint pattern code
sys.path.append("/mnt/ops1p-ren/m/job_access/python/resample_wts/AMSR2")  #contains the target footprint pattern code
from AMSR2_Antenna_Gain import target_gain


warnings.filterwarnings('ignore')
#turn off the warnings for proj



def calc_footprint_LF_from_nsidc_3km(*,lat,lon,land,grid_lats,grid_lons,diameter_in_km):

    print(np.max(land))

    if np.max(land) > 0.001:
        crs_wgs = proj.CRS.from_string("EPSG:4326")
        #crs_wgs = proj.Proj(init='espg:4326')  # assuming you're using WGS84 geographic
        cust = proj.Proj(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m")

        #apply the local projection to obtain the mesh points in meters.
        xv, yv = proj.transform(crs_wgs, cust, grid_lats, grid_lons) 

        #convert to km
        yv = yv/1000.0
        xv = xv/1000.0

        dist_from_center = np.sqrt(np.square(xv) + np.square(yv))

        # g = Geod(ellps='WGS84')
        # n = len(grid_lats)
        # lons = np.full((n),lon)
        # lats = np.full((n),lat)
        # dist_from_center2 = g.inv(lons,lats,grid_lons,grid_lats)[2]/1000.0


        footprint_amplitude = target_gain(dist_from_center,diameter_in_km = diameter_in_km)
        
        #since these are from the EASE 2.0 grid, the cell sizes should be very close to the same,
        # so don't worry about cell size adjustments

        wted_land = np.nansum(land*footprint_amplitude)/np.nansum(footprint_amplitude)
    else:
        wted_land = 0.0

    return wted_land

if __name__ == '__main__':

    from EASE2_grid_land_mask import read_EASE2_grid_land_mask
    from EASE2_grid import read_EASE2_grid_locations

    footprint_diameter_km = 70.0

    weight_map = np.zeros((720,720),dtype=np.float32)

    # if restart:
    #     nc_file = 'L:/access/land_water/land_fraction_1440_721_30km.nc'
    #     land_mask_DS=xr.open_dataset(nc_file)  
    #     weight_map = land_mask_DS['land_fraction'].values

    lats,lons = read_EASE2_grid_locations(resolution=3,pole='south')
    lats_radians = np.deg2rad(lats)
    lons_radians = np.deg2rad(lons)

    land_present = read_EASE2_grid_land_mask(resolution=3,pole='south',convert_ice=True)

    lats_25,lons_25 = read_EASE2_grid_locations(resolution=25,pole='south')

    lats_25 = lats_25.astype(np.float32)
    lons_25 = lons_25.astype(np.float32)

    for ix in range(0,720):
        for iy in range(0,720):
            latitude0 = lats_25[iy,ix]
            longitude0 = lons_25[iy,ix]

            if np.isfinite(latitude0+longitude0):
                latitude0_radians = np.deg2rad(latitude0)
                longitude0_radians = np.deg2rad(longitude0)

                dlat = lats_radians - latitude0_radians
                dlon = lons_radians - longitude0_radians

                #first guess based on delta_lat,delta_lon
                max_lat_diff = 2.0
                max_lon_diff = max_lat_diff/np.cos(latitude0_radians)
                ok1 = np.all([(np.abs(dlat) < np.deg2rad(max_lat_diff)),
                                    np.any([(np.abs(dlon) < np.deg2rad(max_lon_diff)),
                                            (np.abs(dlon) > np.deg2rad(360.0-max_lon_diff))],
                            axis=0)],axis=0)

                ok1 = np.abs(dlat) < np.deg2rad(max_lat_diff)

                ok2 = np.any([(np.abs(dlon[ok1]) < np.deg2rad(max_lon_diff)),
                            (np.abs(dlon[ok1]) > np.deg2rad(360.0-max_lon_diff))],
                            axis=0)

                
                #better guess based on haversine formula
                a = np.square(np.sin((dlat[ok1][ok2])/2.0)) + np.cos(latitude0_radians)*np.cos((lats_radians[ok1][ok2]))*np.square(np.sin((dlon[ok1][ok2])/2.0))
                dist_approx = 6371*(2.0*np.arctan2(np.sqrt(a),np.sqrt(1.0-a)))
                ok3 = dist_approx <= 3.0*footprint_diameter_km

                print(np.sum(ok1),np.sum(ok2),np.sum(ok3),np.max(dist_approx))
                # pass land_present, grid_lat and grid lon arrays selected by both ok1 and ok2

                wt = calc_footprint_LF_from_nsidc_3km(lat = latitude0,
                                                    lon =  longitude0,
                                                    land = land_present[ok1][ok2][ok3],
                                                    grid_lats = lats[ok1][ok2][ok3],
                                                    grid_lons = lons[ok1][ok2][ok3],
                                                    diameter_in_km = footprint_diameter_km)

                weight_map[iy,ix] = wt

                print(f'{np.sum(ok1):06d} {np.sum(ok2):05d} {iy:03d} {ix:03d} {latitude0:7f}, {longitude0:7f}: {wt:7f}')

        ix_array = np.arange(0,720)
        iy_array = np.arange(0,720)
        land_mask_DS = xr.Dataset(
            data_vars = dict(
            land_fraction = (["Y","X"],weight_map),
            latitude = (["Y","X"],lats_25),
            longitude = (["Y","X"],lons_25),
                            ),
                            coords = {
                                "Y":iy_array,
                                "X":ix_array,
                                }
                            )

        if os.name == 'nt':
            nc_file = 'L:/access/land_water/land_fraction_1440_721_70km.from_nsidc_3km_mask.south.ease25.v5.nc'
        elif  os.name == 'posix':
            nc_file = '/mnt/ops1p-ren/l/access/land_water/land_fraction_1440_721_70km.from_nsidc_3km_mask.south.ease25.v5.nc'
        else:
            raise ValueError
        
        encoding = {"land_fraction":{'zlib' : True, 'complevel': 4 }}
        land_mask_DS.to_netcdf(nc_file,encoding=encoding)   
        print(f'Finished Lat = {latitude0}') 
        print(f'Wrote: {nc_file}')    

