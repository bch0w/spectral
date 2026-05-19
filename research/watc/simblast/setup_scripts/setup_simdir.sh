# Make SimBlast Directory Files
# SHELOB
# BIN="/home/bhchow/REPOS/SPECFEM/shelob/specfem3d/bin_26-04-02_c172806_t10_gll5/"
# SHADOWFAX
# BIN="/home/bhchow/REPOS/SPECFEM/shadowfax/specfem3d/bin_26-04-07_c172806_t10_gll5/"
# CHINOOK
BIN="/import/home/bhchow/REPOS/specfem3d/bin_26-04-16_21e40fd_t10_gll5"

if [ ! -f ./bin ]; then
	ln -s ${BIN} ./bin
fi
if [ ! -f ./runscripts ]; then
    # SHADOWFAX/SHELOB
	#  ln -s /home/bhchow/REPOS/simutils/cluster/runscripts/shadowfax/specfem3d ./runscripts
    # CHINOOK
    ln -s /import/home/bhchow/REPOS/simutils/cluster/runscripts/chinook/specfem3d ./runscripts
fi

mkdir -p DATA
mkdir -p LOGS
mkdir -p RESULTS
mkdir -p OUTPUT_FILES/DATABASES_MPI

cd DATA
if false; then
    # SHADOWFAX/SHELOB
    cp /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/Par_files/Par_file_SIMBLAST ./Par_file
    ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/CMTSOLUTIONS/paper_events/ .
    ln -s paper_events/CMTSOLUTION_ISO ./CMTSOLUTION
    ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/STATIONS/STATIONS_PAPER_NK_GRID ./STATIONS 
else
    # CHINOOK
    cp /import/home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/Par_files/Par_file_SIMBLAST ./Par_file
    ln -s /import/home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/CMTSOLUTIONS/paper_events/ .
    ln -s paper_events/CMTSOLUTION_TEST ./CMTSOLUTION
    ln -s /import/home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/STATIONS/STATIONS_PAPER_NK_GRID ./STATIONS 

fi
cd ..

