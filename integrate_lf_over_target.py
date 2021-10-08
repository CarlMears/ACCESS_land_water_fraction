import sys
sys.path.append("m:/job_access/python/rebinning")
sys.path.append("m:/job_access/python/projections")

import numpy as np 
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from libtiff import TIFF
import pyproj as proj
import warnings
  
from bin_ndarray import bin_ndarray
from gauss_kruger import local_geogauss

warnings.filterwarnings('ignore')

def increment_lon_corner(lon_corner:int):
    x = lon_corner+10
    if x > 170:
        x = x-360
    return x

def decrement_lon_corner(lon_corner:int):
    x = lon_corner-10
    if x < -180:
        x = x+360
    return x

def increment_lat_corner(lat_corner:int):
    x = lat_corner+10
    if x > 90:
        raise ValueError('Cannot increment lat = 80')
    return x

def decrement_lat_corner(lat_corner:int):
    x = lat_corner-10
    if x < -90:
        raise ValueError('Cannot increment lat = 90')
    return x


def find_Hansen_land_mask_filename(*,corner_lat,corner_lon,datapath='L:/access/hansen_land_mask/'):
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

    return f'{datapath}Hansen_GFC2015_datamask_{lat_string}_{lon_string}.regrid2.nc'

if __name__ == '__main__':

    warnings.filterwarnings('ignore')

    weight_map = np.zeros((721,1440))
    for lat_index,latitude0 in enumerate(np.arange(-90.0,90.1,0.25)):
        for lon_index,longitude0 in enumerate(np.arange(0.0,360.0,0.25)):

            if abs(latitude0) > 75.0:
                weight_map[lat_index,lon_index] = np.nan
                continue

            # if abs(latitude0 - 51) > 2.0:
            #     weight_map[lat_index,lon_index] = np.nan
            #     continue

            # if abs(longitude0 - 2) > 1.5:
            #     weight_map[lat_index,lon_index] = np.nan
            #     continue


            #figure out which land fraction tiles to load
            corner_lat = int(10*(np.floor(latitude0/10.0)+1))
            corner_lon = int(10*np.floor(longitude0/10))

            num_lons = 1000
            num_lats = 1000
            lon_scale_factor = 100
            lat_scale_factor = 100
            
            half_numy = 100  #eventually, this will be a function of the footprint size
            half_numx = int(np.ceil(half_numy/np.cos(np.deg2rad(latitude0))))


            ilat0 = np.floor((corner_lat - latitude0)*lat_scale_factor).astype(np.int32)
            ilon0 = np.floor((longitude0 - corner_lon)*lon_scale_factor).astype(np.int32)

            #no create a 2x larger map of land_fract
            
            if ilon0<num_lons/2:  #we'll add maps to the left
                corner_lons = [decrement_lon_corner(corner_lon),corner_lon]
                ilon0 = ilon0 + num_lons
            else:
                corner_lons = [corner_lon,increment_lon_corner(corner_lon)]

            if ilat0>num_lats/2:
                try:
                    corner_lats = [corner_lat,decrement_lat_corner(corner_lat)]
                    #ilat0 = ilat0 + num_lats
                except ValueError:
                    corner_lats = [corner_lat]
            else:
                try:
                    corner_lats = [increment_lat_corner(corner_lat),corner_lat]
                except ValueError:
                    corner_lats = [corner_lat]
                
                ilat0 = ilat0 + num_lats

            land_frac = np.zeros((len(corner_lats)*num_lats,len(corner_lons)*num_lons))
            land_frac_lons = np.zeros(len(corner_lons)*num_lons)
            land_frac_lats = np.zeros(len(corner_lats)*num_lats)
            

            for i in [0,1]:
                for j in [0,1]:
                    corner_lat = corner_lats[j]
                    corner_lon = corner_lons[i]
                    nc_land_fraction_file = find_Hansen_land_mask_filename(corner_lat=corner_lat,corner_lon=corner_lon)
                    try:
                        ds = xr.open_dataset(nc_land_fraction_file)
                        land_frac_temp = ds.land_fraction.values
                        land_frac_lons_temp = ds.Longitude.values
                        land_frac_lats_temp = ds.Latitude.values
                    except FileNotFoundError:

                        land_frac_temp = np.zeros((num_lats,num_lons))
                        land_frac_lons_temp = corner_lon+0.5/lon_scale_factor+np.arange(1000)/lon_scale_factor
                        land_frac_lons_temp[land_frac_lons_temp>180.]=land_frac_lons_temp[land_frac_lons_temp>180.]-360.0
                        land_frac_lats_temp = corner_lat-0.5/lat_scale_factor-np.arange(1000)/lat_scale_factor
                    
                    land_frac[j*1000:(j+1)*1000,i*1000:(i+1)*1000] = land_frac_temp
                    land_frac_lats[j*1000:(j+1)*1000] = land_frac_lats_temp
                    land_frac_lons[i*1000:(i+1)*1000] = land_frac_lons_temp
                    

            

            ilon_extent = [ilon0-half_numx,ilon0+half_numx+1]
            ilat_extent = [ilat0-half_numy,ilat0+half_numy+1]
            
            #print(f'ilat bounds {ilat_extent[0]}, {ilat_extent[1]}')
            #print(f'ilon bounds {ilon_extent[0]}, {ilon_extent[1]}')


            submask = land_frac[ilat0-half_numy:ilat0+half_numy+1,
                                ilon0-half_numx:ilon0+half_numx+1]
            ilat_submask = land_frac_lats[ilat0-half_numy:ilat0+half_numy+1]
            ilon_submask = land_frac_lons[ilon0-half_numx:ilon0+half_numx+1]

            if np.nanmax(submask) > 0.0:
           
                #construct a local projection, centered at the center of the target gaussian

                crs_wgs = proj.Proj(init='epsg:4326')  # assuming you're using WGS84 geographic
                cust = proj.Proj(f"+proj=aeqd +lat_0={latitude0} +lon_0={longitude0} +datum=WGS84 +units=m")

                #print(ilat0,ilon0)

                latv,lonv = np.meshgrid(ilat_submask,ilon_submask,indexing='ij')

                #apply the local projection to obtain the mesh points in meters.
                xv, yv = proj.transform(crs_wgs, cust, lonv, latv)   

                #convert to km
                yv = yv/1000.0
                xv = xv/1000.0

                #target size in km
                sigmax = 30.0
                sigmay = 30.0

                # calculate the footprint
                gauss = np.exp(-((np.square(xv)/(2.0*np.square(sigmax))) + 
                            + (np.square(yv)/(2.0*np.square(sigmay)))))

                # calculate the adjustments for changes in grid size
                cell_width_temp  =  xv[1:2*half_numy,2:2*half_numx+1] - xv[1:2*half_numy,0:2*half_numx-1]
                cell_height_temp =  yv[0:2*half_numy-1,1:2*half_numx] - yv[2:2*half_numy+1,1:2*half_numx]

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

                #do the integration, weighted by the gaussian and the grid cell size
                wted_tot = np.sum(submask*gauss*cell_area)/np.sum(gauss*cell_area)
                wt = gauss*cell_area
                # print(np.max(wt[:,0])/np.max(wt),
                #     np.max(wt[:,2*half_numx])/np.max(wt),
                #     np.max(wt[0,:])/np.max(wt),
                #     np.max(wt[2*half_numy,:])/np.max(wt))
                print(f"lat = {latitude0}, lon = {longitude0}, Gaussing Weighted Land Fraction = {wted_tot}")
                weight_map[lat_index,lon_index] = wted_tot
                #make the figure....
                #fig = plt.imshow(submask*gauss*cell_area)
                #png_file = f'M:/job_access/python/land_water/plots/wted_land_frac_{latitude0:06.2f}_{longitude0:06.2f}.png'
                
                #plt.savefig(png_file)
            else:
                print(f"lat = {latitude0}, lon = {longitude0}, Gaussing Weighted Land Fraction = 0.0")
                weight_map[lat_index,lon_index] = 0.0
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

        nc_file = 'M:/job_access/python/land_water/land_fraction_1440_721_30km.nc'
        encoding = {"land_fraction":{'zlib' : True, 'complevel': 4 }}
        land_mask_DS.to_netcdf(nc_file,encoding=encoding)   
        print(f'Finished Lat = {latitude0}') 
        print(f'Wrote: {nc_file}')    

    print()

