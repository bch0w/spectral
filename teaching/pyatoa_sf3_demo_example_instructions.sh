# Pyatoa + SeisFlows3 demo: 
# What it does:
#   Create conda environment, install SeisFlows3, Pyatoa and dependencies. 
#   Download and compile SPECFEM2D and run workflow within the example problem.
# What you need to do:
#   Please change the paths REPO_DIR and WORK_DIR to directories you see fit

# vvv Adjust paths here vvv

# REPO_DIR: Any semi-permanent directory where you would install GitHub 
#   repositories using 'git clone'
# WORK_DIR: Work directory where you can run the example. SeisFlows3 will 
#   generate its own files here so preferably this directory is EMPTY.
REPO_DIR="$HOME/demo_REPOSITORIES" 
WORK_DIR="$HOME/Work/pyatoa_seisflows3_demo" 

# ^^^ Adjust paths here ^^^

# Ensure directories are made
mkdir -p $REPO_DIR
mkdir -p $WORK_DIR

# Exit if any command fails so that we don't install packages into normal 
# environment or try to run an example without all packages installed
set -e  

# Conda create environment, will fail if you don't have Conda installed
conda create -n seisflows python=3.7 -y
source activate seisflows
# conda activate seisflows  # also works but may require `conda init`

# Download SeisFlows3 and Pyatoa from GitHub, install in developer mode
cd $REPO_DIR
git clone --branch devel https://github.com/bch0w/seisflows3.git
cd seisflows3
pip install -e .

cd $REPO_DIR
git clone --branch devel https://github.com/bch0w/pyatoa.git
cd pyatoa
conda install basemap -y  
pip install -e .

# Run the Pyatoa + SeisFlows3 + SPECFEM2D example problem
cd $WORK_DIR
seisflows examples run 2

echo "SeisFlows3 Pyatoa demo execution, results in: $WORK_DIR"
