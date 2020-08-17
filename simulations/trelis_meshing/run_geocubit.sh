#!/bin/bash
# SET PATHS ON CARBON
export TRELISHOME=/opt/Trelis-16.1
export CUBITLIB=/opt/Trelis-16.1/bin:opt/Trelis-16.1/structure
export CUBITDIR=/opt/Trelis-16.1
export CUBITHOME=/opt/Trelis-16.1
export LD_LIBRARY_PATH=$CUBITDIR/bin
export PATH=$PATH:$CUBITDIR/bin
export PYTHONPATH=~/.conda/envs/mesher/bin/python:$CUBITDIR/bin:$CUBITDIR/structure

RUNDIR=`pwd -P`
MESH0_MERGE1=$1

# Set the user parameters
FTAG=nz_north_68_3triples
CFG=${RUNDIR}/${FTAG}/${FTAG}.cfg
GEOCUBIT="/seis/prj/fwi/bchow/packages/GEOCUBIT/GEOCUBIT.py"
MAKENEWMAT_DIR="/seis/prj/fwi/bchow/spectral/simulations/meshing/matlab_utils"
MAKENEWMAT=${MAKENEWMAT_DIR}/make_new_materials_file.m

# Get the name of the output directory
OUTPUT=`grep "output_dir" ${CFG} | cut -d"=" -f 2`

# CORES AND PROCESSOR NUMBER FROM CONFIG
nxi=$(awk '{if (/number_processor_xi/) print $2}' FS='=' ${CFG})
neta=$(awk '{if (/number_processor_eta/) print $2}' FS='=' ${CFG})
nproc=$[nxi*neta]
echo nxi=${nxi}, neta=${neta} , nproc=${nproc}
sleep 2s

# RUN MESH
if [ ${MESH0_MERGE1} == 0 ]
then
	echo "running geocubit mesh ${FTAG}"
	cd ${RUNDIR}/${FTAG} 
    echo ${RUNDIR}/${FTAG} 
	
	tstart=0
	tend=$(( $nproc - 1 ))
	for ii in $( seq $tstart $tend);
    do    
        ${GEOCUBIT} --build_volume --mesh --cfg="${CFG}" --id_proc=$ii & 
	done
	wait
	echo "meshing complete"
    cd ${RUNDIR}
    echo ${RUNDIR}

# RUN MERGE
elif [ ${MESH0_MERGE1} == 1 ]
then
	echo "running geocubit merge ${FTAG}"
	cd ${RUNDIR}/${FTAG}/${OUTPUT}

	# MERGE
	${GEOCUBIT} --collect --merge --meshfiles=mesh_vol_*.e \
                                                --cpux="${nxi}" --cpuy="${neta}"
	sleep 2s
	echo "MERGING COMPLETE"
	echo "MOVING FLUFF TO merge_output"
    MERGEOUT=${RUNDIR}/${FTAG}/merge_output
	mkdir -p ${MERGEOUT}
	mv blocks* *.dat *.jou *.cub quality* *.log ${MERGEOUT}

	# GENERATE SPECFEM OUTPUTS
	echo "running geocubit export to specfem3D"
	${GEOCUBIT} --export2SPECFEM3D --meshfiles=TOTALMESH_MERGED.e

	echo "GATHERING EXPORT FILES"
	EXPORT=${RUNDIR}/${FTAG}/export_mesh
    mkdir -p ${EXPORT}
	mv mesh_file materials_file nodes_coords_file free_or_absorbing_surface_file_zmax absorbing_surface_file_* nummaterial_velocity_file ${EXPORT}
    cd ${EXPORT}

	echo "making new materials file for external tomography files"
	cp ${MAKENEWMAT} ${EXPORT}
	matlab -nojvm -nodesktop -r 'try; make_new_materials_file; catch; end; quit;'
    mv ${EXPORT}/materials_file ${EXPORT}/materials_file_old
	mv ${EXPORT}/materials_file_tomo ${EXPORT}/materials_file
	rm ${EXPORT}/make_new_materials_file.m

	echo "FIN"
	
fi
