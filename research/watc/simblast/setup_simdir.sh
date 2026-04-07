# Make SimBlast Directory Files
if [ ! -f ./bin ]; then
	ln -s /home/bhchow/REPOS/specfem3d/bin_26-04-02_c172806_t10_gll5/ ./bin
fi
if [ ! -f ./runscripts ]; then
	ln -s /home/bhchow/REPOS/simutils/cluster/runscripts/shadowfax/specfem3d ./runscripts
fi
ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/bash_scripts/run_all_simblast.sh .

mkdir -p DATA
mkdir -p LOGS
mkdir -p OUTPUT_FILES/DATABASES_MPI
mkdir -p OUTPUT_FILES/MOVIES

cd DATA
cp /home/bhchow/work/simblasts/simulations/Par_file ./Par_file
ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/CMTSOLUTIONS/paper_events/ .
ln -s paper_events/CMTSOLUTION_ISO ./CMTSOLUTION
ln -s /home/bhchow/REPOS/spectral/research/watc/simblast/SPECFEM_DATA/STATIONS/STATIONS_PAPER_NK_GRID ./STATIONS 
cd ..

