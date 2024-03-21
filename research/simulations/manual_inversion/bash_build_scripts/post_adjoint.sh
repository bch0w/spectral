# If OUTPUT_FILES are given in a SPECFEM master folder, check 
# If it is a finished run, move it to the storage folder, if it is not finished
# then don't do anything
# CALL THIS FUNCTION INSIDE YOUR RUNFOLDER
# OR CALL INSIDE RUN SCRIPT

# STORAGE=/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/storage
RUNFOLDER=`pwd -P`
OUTPUT_FILES=${RUNFOLDER}/OUTPUT_FILES

# If given CMTSOLUTION file, get event_id from there, else require user input
if [ -f ${OUTPUT_FILES}/CMTSOLUTION ]
then
	EVENT_ID="$(grep "event name:" ${OUTPUT_FILES}/CMTSOLUTION | cut -c18-)"
else
	EVENT_ID=$1
	if [ -z "$1" ]
	then
		echo "EVENT ID REQUIRED"
		exit
	fi
fi

# CHECK FOR OUTPUT FOLDER
if [ -d ${OUTPUT_FILES} ]
then
	if [ -f ${OUTPUT_FILES}/output_solver.txt ]
	then
		if echo | grep -q "End of the simulation" ${OUTPUT_FILES}/output_solver.txt
		then
			# assumes CMTSOLUTION files are created with the same format always
			STORAGE=${OUTPUT_FILES}/${EVENT_ID}/
			mkdir -vp ${STORAGE}
			# if [ -d ${STORAGE} ]
			# then
			#  	echo OUTPUT_FILE EXISTS IN STORAGE, SORT IT OUT
			# 	exit
			# fi
			echo MOVING CONTENTS OF ${OUTPUT_ID} TO ${STORAGE}
			mv ${OUTPUT_FILES}/timestamp* ${STORAGE}
			mv ${OUTPUT_FILES}/starttimeloop.txt ${STORAGE}
			mv ${OUTPUT_FILES}/sr.vtk ${STORAGE}
			mv ${OUTPUT_FILES}/output_list_sources.txt ${STORAGE}
			mv ${OUTPUT_FILES}/output_list_stations.txt ${STORAGE}
			mv ${OUTPUT_FILES}/output_solver.txt ${STORAGE}/output_solver_adjoint.txt
			mv ${OUTPUT_FILES}/DATABASES_MPI/*kernel.bin ${STORAGE}
			if [ -L ${OUTPUT_FILES}/DATABASES_MPI/proc000000_save_forward_arrays.bin ]
			then
				rm ${OUTPUT_FILES}/DATABASES_MPI/proc*_forward_arrays.bin ${STORAGE}
				rm ${OUTPUT_FILES}/DATABASES_MPI/proc*_absorb_field.bin ${STORAGE}
			else
				mv ${OUTPUT_FILES}/DATABASES_MPI/proc*_forward_arrays.bin ${STORAGE}
				mv ${OUTPUT_FILES}/DATABASES_MPI/proc*_absorb_field.bin ${STORAGE}
			fi
		else
			echo output_solver.txt HAS NO END OF SIMULATION LINE, STILL RUNNING?
			exit
		fi
    else
		echo NO output_solver.txt FILE, DISREGARDING
	fi
fi
