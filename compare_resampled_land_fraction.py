
import argparse
import datetime
from pathlib import Path

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from rss_plotting.global_map import plot_global_map

plot_path = Path('L:/access/land_water/plots/compare/')

land_path = land_path = Path("L:/access/land_water")
land_file = land_path / "land_fraction_1440_721_30km.combined_hansen_nsidc.nc"
land_fraction_xr = xr.open_dataset(land_file)
land_fraction_hansen_np = land_fraction_xr["land_fraction"].values



land_file = land_path / "resampled.modislandwater.nc"
land_fraction_xr = xr.open_dataset(land_file)
land_fraction_modis_np = land_fraction_xr["resampled_modis_land_water_mask"].values
fig0,ax0 = plot_global_map(land_fraction_modis_np,vmin=0.0,vmax=1.0,cmap='viridis',plt_colorbar=True,title='Land Fraction, MODIS',plt_coastlines=False)
png_file = plot_path / 'MODIS_LF_raw.png'
fig0.savefig(png_file,bbox_inches='tight')

land_fraction_modis_np[~np.isfinite(land_fraction_modis_np)] = 0.0
land_fraction_modis_np[0:40,:] = 1.0


fig1,ax1 = plot_global_map(land_fraction_hansen_np,vmin=0.0,vmax=1.0,cmap='viridis',plt_colorbar=True,title='Land Fraction, Hansen/NSIDC',plt_coastlines=False)
png_file = plot_path / 'Hansen_LF.png'
fig1.savefig(png_file,bbox_inches='tight')

fig2,ax2 = plot_global_map(land_fraction_modis_np,vmin=0.0,vmax=1.0,cmap='viridis',plt_colorbar=True,title='Land Fraction, MODIS',plt_coastlines=False)
png_file = plot_path / 'MODIS_LF_filled.png'
fig2.savefig(png_file,bbox_inches='tight')


fig3,ax3 = plot_global_map(land_fraction_modis_np-land_fraction_hansen_np,vmin=-1.0,vmax=1.0,cmap='BrBG',plt_colorbar=True,title='Land Fraction, MODIS - Hansen/NSIDC',plt_coastlines=False)
png_file = plot_path / 'LF_diff.png'
fig3.savefig(png_file,bbox_inches='tight')


fig4,ax4 = plot_global_map(land_fraction_modis_np-land_fraction_hansen_np,vmin=-0.1,vmax=0.1,cmap='BrBG',plt_colorbar=True,title='Land Fraction, MODIS - Hansen/NSIDC',plt_coastlines=False)
png_file = plot_path / 'LF_diff_sensitive_cmap.png'
fig4.savefig(png_file,bbox_inches='tight')

fig5,ax5 = plot_global_map(land_fraction_modis_np,vmin=0.0,vmax=0.1,cmap='viridis',plt_colorbar=True,title='Land Fraction, MODIS',plt_coastlines=False)
png_file = plot_path / 'MODIS_LF_filled_sensitive_cmap.png'
fig5.savefig(png_file,bbox_inches='tight')

fig6,ax6 = plot_global_map(land_fraction_hansen_np,vmin=0.0,vmax=0.1,cmap='viridis',plt_colorbar=True,title='Land Fraction, Hansen/NSIDC',plt_coastlines=False)
png_file = plot_path / 'Hansen_sensitive_cmap.png'
fig6.savefig(png_file,bbox_inches='tight')
plt.show()
print
