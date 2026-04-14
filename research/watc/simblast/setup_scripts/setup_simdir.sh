# Make SimBlast Directory Files
if false;  # SHELOB
then
	BIN="/home/bhchow/REPOS/SPECFEM/shelob/specfem3d/bin_26-04-02_c172806_t10_gll5/"
else  # SHADOWFAX
	BIN="/home/bhchow/REPOS/SPECFEM/shadowfax/specfem3d/bin_26-04-07_c172806_t10_gll5/"
fi

if [ ! -f ./bin ]; then
	ln -s ${BIN} ./bin
fi
if [ ! -f ./runscripts ]; then
	ln -s /home/bhchow/REPOS/simutils/cluster/runscripts/shadowfax/specfem3d ./runscripts
fi

mkdir -p DATA
mkdir -p LOGS
mkdir -p OUTPUT_FILES/DATABASES_MPI
mkdir -p OUTPUT_FILES/MOVIES

cd DATA
cp /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/Par_files/Par_file_SIMBLAST ./Par_file
ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/CMTSOLUTIONS/paper_events/ .
ln -s paper_events/CMTSOLUTION_ISO ./CMTSOLUTION
ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/STATIONS/STATIONS_PAPER_NK_GRID ./STATIONS 
cd ..

