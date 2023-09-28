## This project contains code used to construct the land fraction layer for the access project.

The output of the project is organized on several Earth-fixed grids:
* 0.25 x 0.25 degree latitude/longitude global grid
* 25 km EASE2 North Pole grid
* 25 km EASE2 South Pole grid

For each grid, there are two possible footprints:
* 30 km diameter circular gaussians
* 70 km diameter circular gaussians

For the six grid/footprint combinations, we need high accuracy land/water fraction values.  We construct these by sampling higher resolution land/water products to the larger footprints.  We use a combinations of land masks from MODIS, and NSIDC.

Ideally, all locations would contain a land fraction obtained by averaging a high-resolution land mask over the gaussian targets used for the brightness temperatures.

The code in this repo combines the Hansen et al. land mask with the 3km land fraction dataset from NSIDC.  The Hansen et al datasets is good over most land areas except Antarctica and the very far north.  I guess this makes sense because it was developed by a forester!  

In these regions, I fill in the data using the NSIDC land masks.

Caveats:
* The accuracy of the Hansen et al. land mask in not known for small islands.
* It might be good to include permanent ice shelves around Antarctica.

Plans:
* Develop a second land fraction dataset based on the MODIS 500m land mask.
