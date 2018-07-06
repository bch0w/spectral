Data directory structure (6.7.18)\

|-kupe [work related to tomography problem: meshing, simulations etc.]\
|------meshGen [scripts for setting up and running Trelis mesher]\
|------------createConfig.py [make the .cfg file that GEOCUBIT needs]\
|------------make_new_materials_file.m [Carl Tapes script, aptly named]\
|------------meshGenTools.py [helper function for designing mesh dimensions]\
|------------runMeshNewZealand.sh [run script for Trelis on personal machine]\
|------tomCat [event catalog for tomography problem]\
|------------tomCat [event catalog, pickle]\
|------------tomCat.csv [csv of above]\
|------------errorCat [potential events that did not process]\
|------------errorCat.csv [csv of above]\
|------tools [misc. tools to assist tomography work]\
|------------ascii2mseed.py [convert ascii specfem outputs to mseed]\
|------------availablestations.py [find available data in my data directories]\
|------------generate_CMTSOLUTION.py [create input CMTSOLUTION for specfem, for me]\
|------------generate_CMTSOLUTION_standalone.py [as above but for general purpose]\
|------------grd2mat.m [convert grd files to matlab data structure]\
|------------grdread2.m [to read in grd files]\
|------adjointbuild.py [pyflex/pyadjoint run script to create adjoint sources]\
|------catbuild.py [generate tomCat for target region and time]\
|------comparemeshes.py [map and waveform plotting for mesh testing]\
|------eventQC.py [event quality control (unfinished)]\
|------obsynth.py [observation/synthetic waveform comparisons]\
|-modules [general helper functions used within this repository]\
|------one_off [single use functions to complete various tasks]\
|------------...\
|------getdata.py [data fetching functions and personal path finding functions]\
|------mapmod.py [mapping functions]\
|------plotmod.py [plotting helper functions]\
|------procmod.py [data processing functions]\
|------synmod.py [synthetic data generation and processing functions]\
|-spectral [frequency domain work: ppsd noise analysis and spetrograms]\
|------duration [long duration resonance work]\
|------------...\
|------createppsd.py [generate ppsd noise analysis plots]\
|------eq_duration.py [quantify and plot earthquake duration]\
|------ppsdplot.py [plot outputs of createppsd]\
|------station_plot.py [deprecated mapper]\
|------waveform_by_event.py [deprecated waveform plotter]\
|-tremor [tremor detection and plotting scripts]\
|------pyfreqscan.py [main processing script for detecting tremors in data]\
|------magnifytremors.py [plot zoomed in sections of tremor windows]\
|------telesearch.py [look for teleseismic events temporally near tremor detect]\
|------plotutils.py [folder specific plotting utilities]\
|------utils.py [folder specific helper functions]\
