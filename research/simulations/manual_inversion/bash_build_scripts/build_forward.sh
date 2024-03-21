#!/bin/bash
# CREATE A NEW SPECFEM RUN BASED ON EVENT ID NUMBER
# SHOULD BE RUN INSIDE THE SPECFEM3DMASTER RUN FOLDER
# EXAMPLE CALL:
# 	./build_forward 2018p130600 -F dry
EVENT_ID=$1

# FLAG FOR change_simulation_type.pl
# -a -- adjoint calculation
# -f -- forward calculation w/ save_forward=.false. (DEFAULT)
# -b -- run both simultaneously
# -F -- forward w/ save_forward=.true.
SIMTYPE="-F"

RUNFOLDER=`pwd -P`
TOMO=/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow
PRIMER=${TOMO}/primer
STORAGE=${RUNFOLDER}/OUTPUT_FILES/${EVENT_ID}
CMTSOLUTION=${PRIMER}/cmtsolution_files/CMTSOLUTION_${EVENT_ID}
TEMPLATE=${PRIMER}/simutils/run_templates/forward_simulation.sh

# ECHO CHECK
echo
echo Event ID: ${EVENT_ID}
echo Run Folder: ${RUNFOLDER}
echo

# CHECK AND EXITS:
# CHECK IF CMTSOLUTION FILE EXISTS
if ! [ -f ${CMTSOLUTION} ]
then
	echo ${CMTSOLUTION} DOES NOT EXIST
	exit
fi
# CHECK IF RUNFOLDER ALREADY EXISTS
if [ -d ${STORAGE}/${EVENT_ID} ]
then
	echo ${EVENT_ID} ALREADY EXISTS IN OUTPUT_FILES
	exit
fi

# CHECK IF OUTPUT_FOLDER EXISTS IN RUNFOLDER
# if [ -d ${RUNFOLDER}/OUTPUT_FILES ]
# then
# 	echo OUTPUT_FILES ALREADY EXISTS IN RUN FOLDER, ATTEMPTING TO MOVE...
# 	#source ${PRIMER}/simutils/output_to_storage.sh
# 	echo
# fi

# IF PASS CHECK-STOPS, CREATE AND RUN
if ! [ -d ${RUNFOLDER}/DATA/tomo_files ]
then
	echo tomo_files IS NOT PRESENT IN DATA
	exit
fi
if ! [ -f ${RUNFOLDER}/DATA/STATIONS ]
then
	echo STATIONS IS NOT PRESENT IN DATA
	exit
fi
if ! [ -d ${RUNFOLDER}/MESH/ ]
then
	echo MESH IS NOT PRESENT
	exit
fi
# IF DATABASES DOESNT EXIST OR IS EMPTY, MAKE IT
if [ ! -d ${RUNFOLDER}/OUTPUT_FILES/DATABASES_MPI ] || [ -z ${RUNFOLDER}/OUTPUT_FILES/DATABASES_MPI ]
then
	echo "DOMAIN DECOMPOSITION AND DATABASE GENERATION REQUIRED"
	exit
fi

echo SYMLINKING CMTSOLUTION
rm ${RUNFOLDER}/DATA/CMTSOLUTION
ln -s ${CMTSOLUTION} ${RUNFOLDER}/DATA/CMTSOLUTION

echo CHANGING SIMULATION TYPE
${RUNFOLDER}/utils/change_simulation_type.pl ${SIMTYPE}
echo

echo CREATING FORWARD RUN SCRIPT: RUNFORWARD_${EVENT_ID}.sh
RF_ID="RUNFORWARD_${EVENT_ID}.sh"
rm ${RUNFOLDER}/${RF_ID}
cp ${PRIMER}/simutils/run_templates/forward_simulation.sh ${RUNFOLDER}/${RF_ID}
SED1="sed -i '3s/.*/#SBATCH --job-name="${EVENT_ID}"_fwd/' ${RUNFOLDER}/${RF_ID}"
SED2="sed -i '19s-.*-cp "${RF_ID}" OUTPUT_FILES/-' ${RUNFOLDER}/${RF_ID}"
SED3="sed -i '71s/.*/ATTENUATION                     = .true./' ${RUNFOLDER}/${RA}"  
eval ${SED1}
eval ${SED2}
eval ${SED3}

echo
echo | grep "SIMULATION_TYPE" ${RUNFOLDER}/DATA/Par_file
echo | grep "SAVE_FORWARD   " ${RUNFOLDER}/DATA/Par_file
echo | grep "NSTEP" ${RUNFOLDER}/DATA/Par_file
echo | grep "DT  " ${RUNFOLDER}/DATA/Par_file
echo | grep "ATTENUATION 	" ${RUNFOLDER}/DATA/Par_file
echo | grep "SAVE_SEISMOGRAMS_DISPLACEMENT" ${RUNFOLDER}/DATA/Par_file
echo | grep "SAVE_SEISMOGRAMS_VELOCITY" ${RUNFOLDER}/DATA/Par_file
echo

