#!/bin/bash
# SET PATHS
export TRELISHOME=/opt/Trelis-16.1
export CUBITLIB=/opt/Trelis-16.1/bin:opt/Trelis-16.1/structure
export CUBITDIR=/opt/Trelis-16.1
export CUBITHOME=/opt/Trelis-16.1
export LD_LIBRARY_PATH=$CUBITDIR/bin
export PATH=$PATH:$CUBITDIR/bin
export PYTHONPATH=~/.conda/envs/mesher/bin/python:$CUBITDIR/bin:$CUBITDIR/structure

# name of mesh
FTAG=mergeTest

MESH0_MERGE1=$1

# CORES AND PROCESSOR NUMBER FROM CONFIG
nxi=$(awk '{if (/number_processor_xi/) print $2}' FS='=' ./${FTAG}/${FTAG}.cfg)
neta=$(awk '{if (/number_processor_eta/) print $2}' FS='=' ./${FTAG}/${FTAG}.cfg)
nproc=$[nxi*neta]
echo nxi=${nxi}, neta=${neta} , nproc=${nproc}
sleep 2s

# RUN MESH
if [ ${MESH0_MERGE1} == 0 ]
then
	echo RUNNING GEOCUBIT TO MESH ${FTAG}
	cd $FTAG 
	
	tstart=0
	tend=$(( $nproc - 1 ))
	for ii in $( seq $tstart $tend); do
	/seis/prj/fwi/bchow/tomo/GEOCUBIT/GEOCUBIT.py --build_volume --mesh --cfg="$FTAG.cfg" --id_proc=$ii & 
	done
	wait
	echo "MESHING COMPLETE"
	cd ../

# RUN MERGE
elif [ ${MESH0_MERGE1} == 1 ]
then
	echo RUNNING GEOCUBIT TO MERGE ${FTAG}
	cd $FTAG

	# MERGE
	/seis/prj/fwi/bchow/tomo/GEOCUBIT/GEOCUBIT.py --collect --merge --meshfiles=mesh_vol_*.e --cpux=$nxi --cpuy=$neta
	sleep 2s
	echo "MERGING COMPLETE"
	echo "GATHERING OUTPUT FILES"
	mkdir mergeOut
	mv blocks* *.log *.jou *.cub quality_* mergeOut/

	# GENERATE SPECFEM OUTPUTS
	echo "GENERATING SPECFEM OUTPUTS"
	/seis/prj/fwi/bchow/tomo/GEOCUBIT/GEOCUBIT.py --export2SPECFEM3D --meshfiles=TOTALMESH_MERGED.e

	echo "GATHERING EXPORT FILES"
	mkdir export
	mv mesh_file materials_file nodes_coords_file free_or_absorbing_surface_file_zmax absorbing_surface_file_* nummaterial_velocity_file export/
	cd ../
fi
