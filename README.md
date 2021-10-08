This project contains code used to construct the land fraction layer for the access project.

Ideally, all locations would contain a land fraction obtained by averaging a high-resolution land mask over the gaussian targets used for the brightness temperatures.

Currently, I am using the Hansen et al. land mask for regions where it is available.  It is good over most land areas except Antarctica and
the very far north.  I guess this makes sense because it was developed by a forester!  

Something needs to be done to fill in the other areas, and to evaluate how well hansen et al does for small islands in, for example, the tropical Pacific.
