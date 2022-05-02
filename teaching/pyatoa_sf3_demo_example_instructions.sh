# Pyatoa + SeisFlows3 demo: Install and run examples
# These are example instructions for installing Pyatoa and SeisFlows3, and 
# running an example problem. Please change the paths to where you see fit

# Adjust paths here
REPO_DIR="~/REPOSITORIES_scratch"  # Wherever you install Python repositories
WORK_DIR="~/Work/sf3_pyatoa_demo"  # An empty directory where we can run the ex.

# Conda create environment, will fail if you don't have Conda installed
conda create -n seisflows3 python=3.7
conda activate seisflows3

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
