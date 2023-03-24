import os
if os.name == 'nt':
    grid_root = 'l:/access/nsidc_grids/'
elif os.name == 'posix':
    grid_root = '/mnt/ops1p-ren/l/access/nsidc_grids/'
else:
    raise ValueError

def read_EASE2_grid_locations(*,resolution,pole,grid_root = grid_root,verbose=False):

    import numpy as np 

    if pole.lower() in ['south','s']:
        pole_str = 's'
    elif pole.lower() in ['north','n']:
        pole_str  = 'N'
    else:
        raise ValueError(f'Pole {pole} not understood')

    if (resolution - 3.0) < 0.0001:
        res_str = '03km'
        size = 6000
    elif (resolution - 6.25) < 0.0001:
        res_str = '6.25km'
        size = 2880
    elif (resolution - 9.0) < 0.0001:
        res_str = '09km'
        size = 2000
    elif (resolution - 12.5) < 0.0001:
        res_str = '12.5km'  
        size = 1440
    elif (resolution - 25.0) < 0.0001:
        res_str = '25km'  
        size = 720
    elif (resolution - 36.0) < 0.0001:
        res_str = '36km'
        size = 500


    lon_file = f'{grid_root}EASE2_{pole_str}{res_str}.lons.{size}x{size}x1.double'
    lat_file = f'{grid_root}EASE2_{pole_str}{res_str}.lats.{size}x{size}x1.double'

    lons = np.fromfile(lon_file)
    lats = np.fromfile(lat_file)

    lons[lons < -998.] = np.nan
    lats[lats < -998.] = np.nan
    lats = np.reshape(lats,(size,size))
    lons = np.reshape(lons,(size,size))
    
    if verbose:
        print(lons.shape)
        print(lats.shape)
        print(np.nanmin(lats),np.nanmax(lats))
        print(np.nanmin(lons),np.nanmax(lons))
    
    return lats,lons

if __name__ == '__main__':

    lats,lons = read_EASE2_grid_locations(resolution=3,pole='north')
    print(lons.shape)   
    print(lats.shape)