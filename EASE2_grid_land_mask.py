
def read_EASE2_grid_land_mask(*,resolution,pole,land_mask_root = 'L:/access/nsidc_land_masks/',verbose=False,convert_ice=True):

    import numpy as np 

    if pole.lower() in ['south','s']:
        pole_str_long = 'south_pole'
        pole_str = 'S'
    elif pole.lower() in ['north','n']:
        pole_str_long = 'north_pole'
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


    land_file = f'{land_mask_root}{pole_str_long}/EASE2_{pole_str}{res_str}.LOCImask_land50_coast0km.{size}x{size}.bin'
    if verbose:
        print(land_file)

    land = np.fromfile(land_file,dtype=np.uint8)
    

    
    land = np.reshape(land,(size,size))
    
    
    if verbose:
        print(land.shape)
        print(np.nanmin(land),np.nanmax(land))

    if convert_ice:
        land_to_return = np.zeros_like(land,dtype=np.float64)
        land_to_return[land < 110] = 1.0
        land_to_return[land == 254] = np.nan
        return land_to_return
    return land

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    land = read_EASE2_grid_land_mask(resolution=3,pole='north')

    plt.imshow(land)
    plt.colorbar()
    plt.show()

    print