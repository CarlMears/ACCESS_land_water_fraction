
def read_hansen_land_fraction_that_contains(lon = 262.23,lat=52.43,verbose = True):
    import numpy as np
    import xarray as xr
    
    corner_lat = 10 + int(10.*np.floor(lat/10.0))
    corner_lon = int(10.*np.floor(lon/10.0))

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

    nc_land_fraction_file = f'L:/access/hansen_land_mask/Hansen_GFC2015_datamask_{lat_string}_{lon_string}.regrid2.nc'
    if verbose:
        print(nc_land_fraction_file)
  
    ds = xr.open_dataset(nc_land_fraction_file)

    return ds

def read_4_hansen_land_fraction_file_that_contains(lon = 262.23,lat=52.43,verbose = True):
    import numpy as np
    import xarray as xr
    
    corner_lat0 = 10 + int(10.*np.floor(lat/10.0))
    corner_lon0 = int(10.*np.floor(lon/10.0))

    corner_lats = np.zeros((2,2),dtype=np.int32)
    corner_lons = np.zeros((2,2),dtype=np.int32)
    
    print(lat,lon)
    print(corner_lat0,corner_lon0)
    
    if (corner_lat0 - lat) > 5.0:
        corner_lats[0,:] = corner_lat0-10
        corner_lats[1,:] = corner_lat0
    else:
        corner_lats[0,:] = corner_lat0
        corner_lats[1,:] = corner_lat0+10

    if (corner_lon0 - lon) > -5.0:
        corner_lons[:,0] = corner_lon0-10
        corner_lons[:,1] = corner_lon0
    else:
        corner_lons[:,0] = corner_lon0
        corner_lons[:,1] = corner_lon0+10

    print(corner_lats)
    print(corner_lons)

    corner_lons[corner_lons > 180] = corner_lons[corner_lons > 180] - 360

    lf_all = np.full((2000,2000),np.nan,dtype=np.float64)
    lats_all = np.full((2000),np.nan,dtype=np.float64)
    lons_all = np.full((2000),np.nan,dtype=np.float64)

    for ilat in [0,1]:
        for ilon in [0,1]:
            corner_lon = corner_lons[ilat,ilon]
            corner_lat = corner_lats[ilat,ilon]
            if corner_lat > 0:
                lat_string = f"{corner_lat:02d}N"
            else:
                if corner_lat == 0:
                    lat_string = "0"
                else:
                    lat_abs = abs(corner_lat)
                    lat_string = f"{lat_abs:02d}S"

            if corner_lon < 0:
                lon_string = f"{abs(corner_lon):03d}W"
            else:
                lon_string = f"{corner_lon:03d}E"   

            nc_land_fraction_file = f'L:/access/hansen_land_mask/Hansen_GFC2015_datamask_{lat_string}_{lon_string}.regrid2.nc'
            if verbose:
                print(nc_land_fraction_file)
  
            ds = xr.open_dataset(nc_land_fraction_file)
            
            lf_all[1000*ilat:1000*(ilat+1),1000*ilon:1000*(ilon+1)] = np.flipud(ds['land_fraction'])
            lats_all[1000*ilat:1000*(ilat+1)] = ds['Latitude'][::-1]
            lons_all[1000*ilon:1000*(ilon+1)] = ds['Longitude']

            print

    ds_out = xr.Dataset(
                coords = {'Latitude':  (['Latitude'],lats_all),
                            'Longitude': (['Longitude'],lons_all)},
                data_vars = {'land_fraction': (['Latitude','Longitude'],lf_all)}
                )
    return ds_out


def get_hansen_land_fraction_local(lon=262.23,lat=52.43,size_km = 50,verbose = True):
    import sys
    sys.path.append("C:/python_packages_cam/local_earth_grid")
    import numpy as np
    import xarray as xr
    from rss_gridding.local_earth_grid import localearthgrid
    from scipy.interpolate import RectBivariateSpline
    import matplotlib.pyplot as plt

    ds = read_hansen_land_fraction_that_contains(lon = lon,lat=lat,verbose = verbose)

    land_frac = np.flipud(ds.land_fraction.values)
    land_frac_lats = np.flip(ds.Latitude.values)
    land_frac_lons = ds.Longitude.values

    grid = localearthgrid.LocalEarthGrid(center_lon = lon,center_lat = lat,
                          delta_x = 1.,delta_y = 1.,
                          num_x = 1+2*size_km,num_y = 1+2*size_km,
                          name='Land Fraction')

    # do the transformation of the lats and lons to x/y space in the transverse mercator projection

    lons_local_grid = grid.lons
    lats_local_grid = grid.lats

    interp_spline = RectBivariateSpline(land_frac_lats,land_frac_lons, land_frac)

    lats = grid.lats
    lons = grid.lons

    lf_on_local_grid = interp_spline(grid.lats,grid.lons,grid=False)

    grid.setdata(lf_on_local_grid,name='Land Fraction')
    grid.addlayer(name='Tb')
    grid.setdata(150 + lf_on_local_grid*100.0,name='Tb')
    return grid

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    lf_grid = get_hansen_land_fraction_local(lon=238.0,lat=38.,size_km = 100)

    ax = lf_grid.contourplot(name='Land Fraction',vmin=0.0,vmax =1.0,
                    cbticks = [-1.0,-0.5,0.0,0.5,1.0],cmap = None,
                    title = 'Land Fraction',plt_contours=False,
                    xlabel='EW distance (km)',ylabel='NS distance (km)')

    print(f'x extent = {lf_grid.x_extent}')
    print(f'y extent = {lf_grid.y_extent}')

    plt.show()