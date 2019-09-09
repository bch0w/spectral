# Set up directory INPUT_KERNELS/ with symlinks to kernel binary finles
# set up a kernels_list.txt to list all the events placed in INPUT_KERNELS/
# Create or check for directory OUTPUT_SUM/ where summed kernels will be sent

# ALIAS CREATION
RUNFOLDER=`pwd -P`
OUTPUT_FILES=${RUNFOLDER}/OUTPUT_FILES
STORAGE=${RUNFOLDER}/OUTPUT_FILES/STORAGE
echo RUNFOLDER: ${RUNFOLDER}

echo SETTING UP INPUT_KERNELS/ DIRECTORY
if [ ! -d ${RUNFOLDER}/INPUT_KERNELS ]
then
	mkdir ${RUNFOLDER}/INPUT_KERNELS
fi

