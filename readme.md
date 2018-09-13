__**Data directory structure (last updated 13.9.18)**__\
(! = files that are not mine, but useful in my work)

|-**simulations [work related to simulation side of tomography problem]**\
|------_build_event_catalog.py_ [generate event catalog for target region and time]\
|------_quality_control.py_ [event quality control (unfinished)]\
|------_obsynth.py_ [observation/synthetic waveform comparisons]\
|------_quick_mseed_comparison.py_ [quick plotter for synthetic-synthetic comparisons]\
|------**mesh_generation_tools [scripts for setting up and running Trelis mesher with Geocubit]**\
|------------_create_geocubit_config.py_ [make the .cfg file that GEOCUBIT needs]\
|------------_mesh_gen_helper.py_ [helper function for designing mesh dimensions]\
|------------_run_geocubit.sh_ [bash run script for Geocubit to call on Trelis on personal machine]\
|------------_compare_meshes.py_ [plot different topography files to compare mesh resolutions, etc.]\
|------------_! make_new_materials_file.m_ [Carl Tapes script, aptly named]\
|------------_! grd2mat.m_ [convert grd files to matlab data structure]\
|------------_! grdread2.m_ [to read in grd files]\
|------**general_tools [misc. tools to assist tomography work]**\
|------------_ascii2mseed.py_ [convert ascii specfem outputs to mseed]\
|------------_availablestations.py_ [find available data in my data directories]\
|------------_estimate_min_res_period.py_ [minimum resolvable period using NGLL differences]\
|------------_generate_CMTSOLUTION.py_ [create input CMTSOLUTION for specfem, for me]\
|------------_generate_CMTSOLUTION_standalone.py_ [as above but for general purpose]\
|------------_make_master_station_list.py_ [create stationxml files to inform pyatoa and specfem]\
|------**event_catalog [earthquake event catalog for tomography problem]**\
|------------_tomCat_ [event catalog, pickle]\
|------------_tomCat.csv_ [csv of above]\
|------------_errorCat_ [potential events that did not process]\
|------------_errorCat.csv_ [csv of above]\
|-**modules [general helper functions used within this repository]**\
|------**one_off [single use functions to complete various tasks]**\
|------------...\
|------_getdata.py_ [data fetching functions and personal path finding functions]\
|------_mapmod.py_ [mapping functions]\
|------_plotmod.py_ [plotting helper functions]\
|------_procmod.py_ [data processing functions]\
|------_synmod.py_ [synthetic data generation and processing functions]\
|-**observations [related to work concerning data analysis of seismic observations]**\
|------_create_ppsd.py_ [generate ppsd noise analysis plots]\
|------_earthquake_durations.py_ [quantify and plot earthquake duration]\
|------_ppsd_plot.py_ [plot outputs of createppsd]\
|------_waveform_by_event.py_ [waveform fetching and plotting]\
|-**tremor [tremor detection and plotting scripts]**\
|------_pyfreqscan.py_ [main processing script for detecting tremors in data]\
|------_magnify_tremors.py_ [plot zoomed in sections of tremor windows]\
|------_telesearch.py_ [look for teleseismic events temporally near tremor detect]\
|------_plotutils.py_ [folder specific plotting utilities]\
|------_utils.py_ [folder specific helper functions]\
