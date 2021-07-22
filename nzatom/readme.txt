NAME:  
	nz_atom_north_chow_etal_2021_vp+vs
TITLE: 
	New Zealand Adjoint TOMography velocity model (North Island)
TYPE:  
	3-D Tomography Earth Model
SUB TYPE: 
	Shear-wave velocity (km/s) and Compressional-wave velocity (km/s)
YEAR: 
	2021
SHORT DESCRIPTION:
	This is a 3D velocity model of the North Island of New Zealand derived using
	earthquake-based adjoint tomography. The starting model is defined as the 
	ray-based NZ-Wide2.2 Velocity Model from Eberhart-Phillips et al. (2021).
 	To derive this velocity model, we iteratively improve fits between earthquake 
	observations from New Zealand-based broadband seismometers and synthetically
	generated waveforms from spectral element simulations. The waveform bandpass
	of interest is 4-30s. 
	This velocity model defines the following quantities: 
	Vp (km/s), Vs (km/s), density (kg/m^3), Qp and Qs, however only Vp and Vs 
	are updated during our inversion. The remaining quantities are defined by 
	the starting velocity model. 
AUTHORS: 
	Bryant Chow, Yoshihiro Kaneko, Carl Tape, Ryan Modrak, John Townend
REFERENCE MODEL:
	NZ-Wide2.2 Velocity Model (https://zenodo.org/record/3779523#.YPi7ClMzbHo)
	Eberhart-Phillips et al. (2021)
DEPTH COVERAGE:
	-3.0 to 400.0 km (below earth surface)
AREA:
	North Island of New Zealand and the Hikurangi subduction zone
	Coordinate system: UTM 60S
		x-axis 17000,635000, y-axis 5286000, 5905000	
	Coordinate system: WGS84
		latitude -42.5,-37.0, longitude 173.0,178.5
DATA SET DESCRIPTION:
	Dataset includes seismic waveforms filtered at 4-30s period.
	Waveforms from 60 earthquakes recorded on up to 88 broadboad three-component 
	seisimic stations, resulting in 1800 unique source-receiver pairs.
	Broadband stations include the following networks:
		-GeoNet Permanent Seismic Network (https://www.geonet.org.nz/)
		-Broadband East Coast Network (http://www.fdsn.org/networks/detail/2P_2017/)
		-Seismic Analysis of the Hikurangi Experiment (https://www.fdsn.org/networks/detail/X2_2009/)
		-Deep Geothermal HADES Seismic Array (http://www.fdsn.org/networks/detail/Z8_2009/)
		-Gisborne-Mahia Seismic Tremor Array (http://www.fdsn.org/networks/detail/ZX_2011/)
DATA FILES DESCRIPTION:
	The velocity model is separated into three ASCII files covering overlapping 
	depth ranges and with differing grid spacing. These files are formatted in
	the SPECFEM3D_Cartesian external tomography file format. 
	
	+ Header: First 4 lines of each file define the following quantities
		x_min y_min z_min x_max y_max z_max
		dx dy dz
		nx ny nz
		vp_min vp_max vs_min vs_max rho_min rho_max
	+ Data: The remainder of the file contains lines with 8 values each defining
		X Y Z Vp Vs rho Qp Qs 
	+ Notes:
		'dx' is the x-axis grid spacing in m defined in the UTM 60S coord. system
		'nx' is the number of points on the x-axis
		'Vp' and 'Vs' are in units of m/s
		'rho' is density in units of kg/m^3
		'shallow' defines depths -3 to 9 km (BSL)
		'crust' defines depths 7 to 47 km (BSL)
		'mantle' defines depths 44 to 400 km (BSL)
CONTACT:
	Bryant Chow (bryant.chow@vuw.ac.nz)

	
	
	
