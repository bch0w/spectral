# Pyatoa + SeisFlows3 demo: Install and run example problem
#
# Please change the paths REPO_DIR and WORK_DIR to directories you see fit

# Adjust paths here
REPO_DIR="$HOME/demo_REPOSITORIES"  # Where you install GitHub repos
WORK_DIR="$HOME/Work/pyatoa_seisflows3_demo"  # Empty directory to run the ex.

# Exit if any command fails so that we don't continue 
set -e  

# Conda create environment, will fail if you don't have Conda installed
conda create -n seisflows python=3.7 -y
source activate seisflows

# Download Pyatoa and SeisFlows3 from GitHub
mkdir $REPO_DIR
cd $REPO_DIR
git clone --branch devel https://github.com/bch0w/seisflows3.git
cd seisflows3
pip install -e .

cd $REPO_DIR
git clone --branch devel https://github.com/bch0w/pyatoa.git
cd pyatoa
pip install -e .

# Run the Pyatoa + SeisFlows3 + SPECFEM2D example problem
mkdir $WORK_DIR
cd $WORK_DIR
seisflows examples run 2
