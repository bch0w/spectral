#!/bin/bash
# TO DO: turn off attenuation, change hessian kernel
# CREATE A NEW SPECFEM RUN BASED ON EVENT ID NUMBER
# SHOULD BE RUN INSIDE THE SPECFEM3DMASTER RUN FOLDER
# EXAMPLE CALL: ./build_adjoint 2018p130600
EVENT_ID=$1
if [ -z "$1" ]
then
	echo "EVENT ID REQUIRED"
	exit
fi

# FLAG FOR change_simulation_type.pl
# -a -- adjoint calculation
# -f -- forward calculation w/ save_forward=.false. (DEFAULT)
# -b -- run both simultaneously
# -F -- forward w/ save_forward=.true.
SIMTYPE="-b"

RUNFOLDER=`pwd -P`
OUTPUT_FILES=${RUNFOLDER}/OUTPUT_FILES
TOMO=/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow
PRIMER=${TOMO}/primer
CMTSOLUTION=${PRIMER}/cmtsolution_files/${EVENT_ID}CMTSOLUTION

# ECHO CHECK
echo
echo Event ID: ${EVENT_ID}
echo Run Folder: ${RUNFOLDER}
echo

echo CHANGING SIMULATION TYPE
${RUNFOLDER}/utils/change_simulation_type.pl ${SIMTYPE}
echo

RA="RUNADJOINT_${EVENT_ID}.sh"
echo CREATING ADJOINT RUN SCRIPT: ${RA}
cp ${PRIMER}/simutils/run_templates/adjoint_simulation.sh ${RUNFOLDER}/${RA}
SED1="sed -i '3s/.*/#SBATCH --job-name="${EVENT_ID}"_adj/' ${RUNFOLDER}/${RA}"
SED2="sed -i '71s/.*/ATTENUATION                        = .false./' ${RUNFOLDER}/${RA}"  
eval ${SED1}
eval ${SED2}

echo SETTING CMTSOLUTION
rm ${RUNFOLDER}/DATA/CMTSOLUTION
if [ -f ${OUTPUT_FILES}/CMTSOLUTION ]
then
	rm ${OUTPUT_FILES}/CMTSOLUTION
fi
ln -s ${CMTSOLUTION} ${RUNFOLDER}/DATA/CMTSOLUTION
ln -s ${CMTSOLUTION} ${OUTPUT_FILES}/CMTSOLUTION

echo
echo CHECKING OUTPUT_FILES
if [ -f ${OUTPUT_FILES}/DATABASES_MPI/proc000000_save_forward_arrays.bin ]
then 
	echo proc*_save_forward_arrays.bin exist in output_files, removing
	rm ${OUTPUT_FILES}/DATABASES_MPI/proc*_save_forward_arrays.bin 
	echo symlinking proc*_absorb_field.bin to DATABASES_MPI
	rm ${OUTPUT_FILES}/DATABASES_MPI/proc*_absorb_field.bin
else
	echo symlinking proc*_save_forward_arrays.bin to DATABASES_MPI
	ln -s ${OUTPUT_FILES}/${EVENT_ID}/proc*_save_forward_arrays.bin ${OUTPUT_FILES}/DATABASES_MPI/
	echo symlinking proc*_absorb_field.bin to DATABASES_MPI
	ln -s ${OUTPUT_FILES}/${EVENT_ID}/proc*_absorb_field.bin ${OUTPUT_FILES}/DATABASES_MPI/
fi
echo
echo SYMLINKING SEM FILES
if [ -d ${RUNFOLDER}/SEM ]
then
	echo removing old SEM/ folder 
	rm -r SEM
fi
if [ -d ${RUNFOLDER}/INPUT_SEM/${EVENT_ID} ]
then
	ln -s ${RUNFOLDER}/INPUT_SEM/${EVENT_ID} ${RUNFOLDER}/SEM
	rm ${RUNFOLDER}/DATA/STATIONS_ADJOINT
	cp ${RUNFOLDER}/INPUT_SEM/${EVENT_ID}/STATIONS_ADJOINT ${RUNFOLDER}/DATA/
fi
echo
echo
echo MAKE SURE par_file PARAMETERS ARE SET APPROPRIATELY
echo
echo | grep "SIMULATION_TYPE" ${RUNFOLDER}/DATA/Par_file
echo | grep "NSTEP" ${RUNFOLDER}/DATA/Par_file
echo | grep "DT  " ${RUNFOLDER}/DATA/Par_file
echo | grep "ATTENUATION     " ${RUNFOLDER}/DATA/Par_file
echo | grep "SAVE_SEISMOGRAMS_*" ${RUNFOLDER}/DATA/Par_file
echo | grep "APPROXIMATE_HESS_KL" ${RUNFOLDER}/DATA/Par_file
echo
echo
