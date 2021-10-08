
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


def get_hansen_land_fraction_local(lon=262.23,lat=52.43,size_km = 50):
    import sys
    sys.path.append("C:/python_packages_cam/local_earth_grid")
    import numpy as np
    import xarray as xr
    import localearthgrid
    from scipy.interpolate import RectBivariateSpline
    import matplotlib.pyplot as plt

    ds = read_hansen_land_fraction_that_contains(lon = lon,lat=lat)

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