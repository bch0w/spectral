Data directory structure (21.8.18)

|-**vizmis [misfit visualization and adjoint source creation]
|------_adjointBuilder.py_ [pyflex/pyadjoint run script to create adjoint sources]\
|------_corkBoard.py_ [defines Cork class to parse through PyASDF files]
|------**viztools [waveform plotters, map makers and data visualization]
|------------_dataDepicter.py_ [class for creating statistical plots]
|------------_mapMaker.py_ [plot source receiver information on maps]
|------------_windowMaker.py_ [visualize waveform data etc.]
|------**tests [example run scripts and unit testing]
|------------_func_test.py_ [all testing functions so far]
|------------**test_data [example data for unit testing etc.]
|-**kupe [work related to simulation side of tomography problem]**\
|------_catbuild.py_ [generate tomCat for target region and time]\
|------_comparemeshes.py_ [map and waveform plotting for mesh testing]\
|------_eventQC.py_ [event quality control (unfinished)]\
|------_obsynth.py_ [observation/synthetic waveform comparisons]\
|------**meshGen [scripts for setting up and running Trelis mesher]**\
|------------_createConfig.py_ [make the .cfg file that GEOCUBIT needs]\
|------------_make_new_materials_file.m_ [Carl Tapes script, aptly named]\
|------------_meshGenTools.py_ [helper function for designing mesh dimensions]\
|------------_runMeshNewZealand.sh_ [run script for Trelis on personal machine]\
|------**tomCat [event catalog for tomography problem]**\
|------------_tomCat_ [event catalog, pickle]\
|------------_tomCat.csv_ [csv of above]\
|------------_errorCat_ [potential events that did not process]\
|------------_errorCat.csv_ [csv of above]\
|------**tools [misc. tools to assist tomography work]**\
|------------_ascii2mseed.py_ [convert ascii specfem outputs to mseed]\
|------------_availablestations.py_ [find available data in my data directories]\
|------------_estimateMRP.py_ [minimum resolvable period using NGLL differences]\
|------------_generate_CMTSOLUTION.py_ [create input CMTSOLUTION for specfem, for me]\
|------------_generate_CMTSOLUTION_standalone.py_ [as above but for general purpose]\
|------------_grd2mat.m_ [convert grd files to matlab data structure, not mine]\
|------------_grdread2.m_ [to read in grd files, not mine]\
|-**modules [general helper functions used within this repository]**\
|------**one_off [single use functions to complete various tasks]**\
|------------...\
|------_getdata.py_ [data fetching functions and personal path finding functions]\
|------_mapmod.py_ [mapping functions]\
|------_plotmod.py_ [plotting helper functions]\
|------_procmod.py_ [data processing functions]\
|------_synmod.py_ [synthetic data generation and processing functions]\
|-**spectral [frequency domain work: ppsd noise analysis and spetrograms]**\
|------_createppsd.py_ [generate ppsd noise analysis plots]\
|------_quakeDuration.py_ [quantify and plot earthquake duration]\
|------_ppsdplot.py_ [plot outputs of createppsd]\
|-**tremor [tremor detection and plotting scripts]**\
|------_pyfreqscan.py_ [main processing script for detecting tremors in data]\
|------_magnifytremors.py_ [plot zoomed in sections of tremor windows]\
|------_telesearch.py_ [look for teleseismic events temporally near tremor detect]\
|------_plotutils.py_ [folder specific plotting utilities]\
|------_utils.py_ [folder specific helper functions]\
