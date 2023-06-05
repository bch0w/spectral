# NZATOM r0.1

June 5, 2023 
Bryant Chow *(bhchow@alaska.edu)*  
[doi.org/10.17611/dp/emc.2021.nzatomnnorthvpvs.1](https://doi.org/10.17611/dp/emc.2021.nzatomnnorthvpvs.1)

This document describes an update to the NZATOM model (revision 0.1) that 
attempts to address numerical artefacts present in the original model available
on the IRIS EMC. The following sections describe the issue and our solution.

## Background

The New Zealand Adjoint TOmography Model (NZATOM) is a 3D velocity model 
derived from a full waveform inversion using the numerical solver SPECFEM3D
Cartesian (SPECFEM3D). SPECFEM3D is a spectral-element solver, which solves 
the seismic wave equation on unstructured hexahedral meshes (deformed cubes). 
These *GLL* (Gauss-Lobatto-Legendre) models stretch and skew cubic elements to
conform to interface boundaries (e.g., topography, moho topography), or 
*coarsening* layers, where elements increase in vertical or lateral size. These
coarsening layers are used to decrease element sizes at depth, where seismic 
wavespeeds typically increase, allowing for less computational requirements
for the same numerical resolution. 

During the inversion procedure, two meshes with varying resolution were used
in order to save computational expense. We call these two meshes *coarse* 
and *fine*. Both coarse and fine meshes had different element sizes and
coarsening layer depths.

Due to decreasing resolution with depth, the final NZATOM model published on 
IRIS EMC is split into three overlapping blocks: shallow (-3--8km), 
crust (7--50km), mantle (44--400km)

![gll_model](https://user-images.githubusercontent.com/23055374/243476592-4078655a-1bef-48de-9341-9b6e1de34b37.png)
*Fig. 1: Side-on view (Y-Z plane) of the NZATOM fine-resolution GLL model showing coarsening layers at approximately 25km and 90km depths.
These coarsening layers lead to artefacts from nearest-neighbor interpolation of the original NZATOM model. Depths are shown in units of meters.*

## Problem

> TL;DR: Original model shows obvious numerical artefacts due to the 
  interpolation of an unregular mesh onto a regular grid at various stages 
  of the inversion workflow.

To obtain NZATOM in a file format usable by non-SPECFEM users, we needed to 
interpolate it onto a regular (XYZ) grid. The original method for doing so
was to use a nearest-neighbor algorithm on the original NZATOM GLL mesh to 
determine model values for a given uniform grid based on the nearest GLL point
available. That is, for a given point **r**(x,y,z) in our regular grid, the 
program would find the nearest GLL point and assign that to point **r**.

This method worked nominally for regular elements, but for skewed elements 
related to coarsening layers the nearest neighbor algorithm tended to create
artefacts as certain points **r** began to show preference for certain element
boundaries as they passed through skewed elements.

During the inversion, this approach was also used to transfer our tomography 
model from the *coarse* to the *fine* mesh. Therefore the coarsening layers
of the coarse mesh were imparted onto the fine mesh. Similarly, extracting the
final velocity model imparted the coarsening layers of the *fine* mesh onto 
the published version of NZATOM.

### Artefact Locations

Artefacts correspond to the coarsening layers of both the coarse and fine mesh:

Coarse Mesh Coarsening Layers: ~15--35km, ~75--125km
Fine Mesh Coarsening Layers: ~14--25km, ~65--100km

## Solution

> TL;DR: We interpolate affected model onto a regular mesh, smooth away 
  artefacts, and re-extract onto regular grid.

To remove these artefacts from the underlying model, while avoiding adding any
new additional artefacts, our approach was to start from a regularized GLL grid.
That is, a GLL mesh with no topography, bathymetry, or coarsening layers, that
shares an origin and domain with the original NZATOM GLL mesh.

We then take the NZATOM model present on IRIS EMC, in a regularized XYZ format 
and containing artefacts from the nearest neighbor approach, and interpolate 
this onto our new regularized mesh. The result is now a regular (cube-shaped
elements) GLL mesh that has model values corresponding to the NZATOM model, and
contains numerical artefacts.

We then apply a 2D Gaussian smoothing to the entire model. Smoothing 
coefficients for each smoothed model are provided in the following subsection.
Smoothing coefficient values were determined through trial-and-error and visual 
inspection of the smoothed model, to determine if artefacts were sufficiently,
without removing real small-scale features in the models. 

The smoothing procedure resulted in significant removal of the numerical 
artefacts, while perserving small-scale heterogeneities (~5km width). Although
some artefacts can still be seen in the resulting model, they should not have
significant affect on future applications of the model. 

### Workflow

1. Generate regular GLL models with SPECFEM for: A) shallow, B) crust, C) mantle
2. Interpolate NZATOM XYZ model from IRIS EMC onto meshes from (1)
3. Smooth SPECFEM GLL models created in (2) using smoothing coefficients listed below
4. Use nearest-neighbor to extract XYZ files from smoothed models of (3)
5. Convert .xyz files to IRIS EMC netCDF format


### Smoothing Coefficients

**Shallow** model (-2.25--8km depth):
- Horizontal Half-Width = .5km
- Vertical Half-width = .5km

**Crustal** model (7--50km depth):
- Horizontal Half-Width = 3km
- Vertical Half-width = 1km

**Mantle** model (44--400km depth):
- Horizontal Half-Width = 4km
- Vertical Half-width = 2km

**Note**: All parameters (vp, vs, rho, qp, qs) are smoothed

## Revised Model

The revised model was again derived by taking nearest-neighbor interpolation 
on a regular grid for each of the sub-models. Now that the underlying GLL model
is regular, there are no artefacts that arise from this approach and the 
resulting XYZ models show smoothly varying values and no artefacting. 

![](https://user-images.githubusercontent.com/23055374/243476681-43106c92-1742-4274-ba4b-2b60f9b4a308.png)
*Fig. 2: Artefact removal, slice through Crust tomography model at X~=171km (Y-Z plane) for Vs in km/s.  
Left: Original NZATOM model showing visible artefacts. Right: Revised (r0.1) NZATOM model smoothed with a 2D Gaussian 
w/ a horizontal half-width of 3km, and a vertical half-width of 1km*

![](https://user-images.githubusercontent.com/23055374/243476805-b3c5fe04-0ca2-4bb1-aa22-b4bfe5950441.png)
*Fig 3: Same as Fig. 2 but for Y=5286km (X-Z plane)*

![](https://user-images.githubusercontent.com/23055374/243476867-11fed8bd-b85b-4501-9bd6-bfc245923d64.png)
*Fig. 4: Same as Fig. 2 but for Z=15km BSL (X-Y plane)*

![](https://user-images.githubusercontent.com/23055374/243477033-7363b743-3cc4-4f57-b173-4996926922c2.png)
*Fig. 5: Same as Fig. 2 but for Z=25km BSL (X-Y plane)*


### Amplitudes

The smoothing procedure will have damped the largest amplitudes of the 
underlying model, resulting in slightly lower minimum and maximum amplitudes for 
all parameters.

### Domain Bounds

Due to a new approach in interpolating the underlying model, the published 
revision has slightly different domain bounds as the original model (smaller). 
The grid spacing is kept the same as the original model.

### Topography

The topography on the original NZATOM GLL model was defined by SRTM-30P, with 
no model values given for values in the air (above topography). Because the 
regular model published to IRIS EMC was given as a regular cube, values in the 
air were extrapolated from the nearest point in the model. This newly revised 
model does the same, so Users are cautioned when querying model parameters above 
topography/bathymetry depth values.

