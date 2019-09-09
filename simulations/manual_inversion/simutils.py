#!/usr/bin/env python3
"""
Simulation Utilities
A collection of utilities for working with Specfem3D Cartesian on an HPC cluster
These scripts and functions glue together the inputs and outputs of various
Fortran codes executed by Specfem3D. Originally a set of bash scripts to move
files around and create run-scripts.

Functions assume that Specfem run folders and files retain a specific format.
We are using Specfem3D Cartesian pulled from github at commit:

6895e2f7f7ad26cd91ff5cdb83bd340fbbe0a46f

Which was committed June 18, 2018

The workflow is as follows:
build_forward > forward simulation > post forward >
build_adjoint > adjoint simulation > post adjoint ...

pre_condition_sum > model_update > post_update

TO DO LIST:
build_forward()
    + add a check for proc**_Database files to see if mesher was run
    + add a check for proc**_external_mesh.bin files to see if generate_db run
    + tomo_files isn't really required past model update 1?
post_forward()
    + if this runs and fails halfway, running it a second time sometimes screws
      things up, need some way to either undo the work or not run if it won't
      complete the run, or dry run?
build_adjoint()
    + SEM folder check?
    + adjoint choice for dynamic filenames?
    + symlinking save_forward_arrays was throwing errors randomly, maybe just
      move the files rather than symlinking?
post_adjoint()
    + should we delete the one off outputs? I think they're the same as forward
pre_condition_sum()
    + kernels_list.txt is listed in constants_tomography.h, dynamically get?
model_update()
    + move mode_update pathnames into dynamic filenames and call from here
post_model_update()
    + check-stop if M00_OUTPUT_FILES exists, otherwise you start making edits to
      the new output_files directory
    + remove the hardcoded mesh_files_m01 at the end of the function
status_check()
    + broken for model updates, kinda
"""
from __future__ import print_function

import os
import sys
import glob


def directories():
    """
    Return a collection of paths that most of the functions will want.
    These are hardcoded and should be standard among runs.
    :return: dict
    """
    # directories with user-information
    home = os.path.expanduser("~")
    primer = os.path.join(home, "tomo", "primer")

    # current specfem3d master folder
    runfolder = os.getcwd()
    data = os.path.join(runfolder, "DATA")
    output_files = os.path.join(runfolder, "OUTPUT_FILES")

    directories_dictionary = {"home": home, "primer": primer, "data": data,
                              "runfolder": runfolder,
                              "output_files": output_files}

    return directories_dictionary


def dynamic_filenames():
    """
    Names of files and folders are set in the Par_file or in constats.h or in
    constants_tomography.h. Rather than hardcoding them in here, dynamically
    fetch them incase the user decides to change them
    :return: dictionary
    """
    drc = directories()

    # grab information from the Par_file
    par_file = os.path.join(drc["data"], "Par_file")
    lines = open(par_file, "r").readlines()
    for line in lines:
        if "TOMOGRAPHY_PATH" in line:  # e.g. DATA/tomo_files/
            tomo_path = line.strip().split()[-1]
            drc["tomo_path"] = tomo_path
        elif "LOCAL_PATH" in line:  # e.g. OUTPUT_FILES/DATABASES_MPI/
            local_path = line.strip().split()[-1]
            drc["local_path"] = local_path
    return drc


def edit_par_file(fid, choice):
    """
    Edits the Par_file in place for forward or adjoint simulations, or update.
    Encompasses the tasks of change_simulation_type.pl, with a few extras.
    Assumes the Par_file will always be the same format (e.g. sim_type   = 1),
    but I don't think Specfem devels will change it during my PhD so okay.
    :type fid: str
    :param fid: path to Par_file
    :type choice: str
    :param choice: forward, adjoint or update
    """
    with open(fid, "r+") as f:
        lines = f.readlines()
        if choice == "forward":
            simulation_type = "1"
            save_forward = ".true."
            attenuation = ".true."
            model = "default"
        elif choice == "adjoint":
            simulation_type = "3"
            save_forward = ".false."
            attenuation = ".false."
            model = "default"
        elif choice == "update":
            simulation_type = "1"
            save_forward = ".false."
            attenuation = ".false."
            model = "gll"
  
        # Change parameters by replacing the old values
        # assuming that the Par_file retains the same format
        for i, line in enumerate(lines):
            if "SIMULATION_TYPE" in line:
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, simulation_type)
                print("\t", lines[i].strip())
            elif "SAVE_FORWARD" in line:
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, save_forward)
                print("\t", lines[i].strip())
            elif "MODEL " in line and not "_MODEL" in line:
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, model)
                print("\t", lines[i].strip())
            elif ("ATTENUATION " in line) and ("_ATTENUATION" not in line):
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, attenuation)
                print("\t", lines[i].strip())
                break

        # delete original file, set cursor to beginning and write in new pars
        f.truncate(0)
        f.seek(0)
        f.writelines(lines)
        print("overwrote Par_file")


def build_forward(event_id):
    """
    Prepare the run folder for a forward simulation for a given event id

    :type: event_id: str
    :param event_id: event id to grab the correct cmtsolution files etc
    """
    drc = dynamic_filenames()

    # check if all the prerequisite files are included
    cmtsolution = os.path.join(
        drc["primer"], "cmtsolution_files", "{}CMTSOLUTION".format(event_id)
    )
    cmtsolution_check = os.path.exists(cmtsolution)
    tomofiles_check = os.path.exists(drc["tomo_path"])
    databases_check = os.path.exists(drc["local_path"])
    stations_check = os.path.exists(os.path.join(drc["data"], "STATIONS"))
    mesh_check = os.path.exists(os.path.join(drc["runfolder"], "MESH"))

    # Check-exits to make sure the run folder is set correctly
    if not databases_check:
        print("{} does not exist".format(local_path))
        sys.exit()
    if not cmtsolution_check:
        print("{}CMTSOLUTION does not exist in primer".format(event_id))
        sys.exit()
    if not tomofiles_check:
        print("{tomofiles} does not exist in {data}".format(
            tomofiles=drc["tomo_path"], data=drc["data"])
        )
        sys.exit()
    if not stations_check:
        print("STATIONS does not exist in DATA/")
        sys.exit()
    if not mesh_check:
        print("MESH/ does not exist in RUNFOLDER")
        sys.exit()

    # set up the run folder
    print("symlinking CMTSOLUTION")
    cmtsolution_destination = os.path.join(drc["data"], "CMTSOLUTION")
    if os.path.exists(cmtsolution_destination):
        os.remove(cmtsolution_destination)
    os.symlink(cmtsolution, cmtsolution_destination)

    print("editing Par_file")
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="forward")

    print("generating forwardrun.sh file")
    fid_in = os.path.join(
        drc["primer"], "simutils", "run_templates", "forward_simulation.sh"
    )
    fid_out = os.path.join(drc["runfolder"], "forwardrun.sh".format(event_id))
    if os.path.exists(fid_out):
        os.remove(fid_out)

    def generate_runforward(fid_in, fid_out, event_id):
        """
        generate a forward run script to be called by sbatch
        """
        lines = open(fid_in, "r").readlines()
        for i, line in enumerate(lines):
            if "${EVENT_ID}" in line:
                lines[i] = line.replace("${EVENT_ID}", event_id)
        with open(fid_out, "w") as f_out:
            f_out.writelines(lines)

    generate_runforward(fid_in, fid_out, event_id)

    print("forward build complete")


def event_id_from_cmt(fid):
    """
    occasionally need to get the event_id from the cmtsolution file
    :type fid: str
    :param fid: path location of the file id
    :rtype: str
    :return: event id from cmtsolution
    """
    lines = open(fid, "r").readlines()
    for line in lines:
        if "event name" in line:
            return line.strip().split()[-1]


def post_forward():
    """
    After a forward simulation, files need to be moved around so that local path
    and output_files can be free for a new forward or adjoint simulation.
    """
    import shutil

    # set some recurring directory names
    drc = dynamic_filenames()
    output = drc["output_files"]
    output_solver = os.path.join(output, "output_solver.txt")
    input_sem = os.path.join(drc["runfolder"], "INPUT_SEM")
    cmtsolution = os.path.join(output, "CMTSOLUTION")
    storage = os.path.join(output, "STORAGE")

    # check the output_solver.txt file to make sure simulation has finished
    if os.path.exists(output_solver):
        text = open(output_solver).read()
        if "End of the simulation" not in text:
            sys.exit("Simulation has not finished, disregarding")
    else:
        sys.exit("No output_solver.txt file, disregarding")

    # get event id from cmtsolution
    if not os.path.exists(cmtsolution):
        sys.exit("CMTSOLUTION doesn't exist in OUTPUT_FILES")
    event_id = event_id_from_cmt(cmtsolution)

    # create storage folder for all the event specific outputs
    event_storage = os.path.join(storage, event_id)
    if not os.path.exists(event_storage):
        if not os.path.exists(storage): 
            os.mkdir(storage)
        print("creating storage directory")
        os.mkdir(event_storage)

    # move files from OUTPUT_FILES/ to STORAGE/event_id
    for quantity in ["*.sem?", "timestamp*"]:
        print("moving {} files".format(quantity))
        for fid in glob.glob(os.path.join(output, quantity)):
            shutil.move(fid, event_storage)

    # the random one-off output files
    print("moving misc. solver outputs")
    for quantity in ["starttimeloop.txt", "sr.vtk", "output_list_sources.txt",
                     "output_list_stations.txt", "Par_file", "CMTSOLUTION",
                     "output_solver.txt"]:
        try:
            shutil.move(os.path.join(output, quantity), event_storage)
        except FileNotFoundError:
            print("{} not found".format(quantity))
            continue
    
    # proc*_save_forward_arrays.bin and absorb_field files
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        i = 0
        print("moving {} files".format(quantity), end="...")
        for fid in glob.glob(os.path.join(drc["local_path"], quantity)):
            if os.path.exists(fid):
                shutil.move(fid, event_storage)
                i += 1
        print("{} files moved".format(i))
    
    # create event folder in INPUT_SEM/ to be populated by adjoint sources
    event_sem = os.path.join(input_sem, event_id)
    if not os.path.exists(event_sem):
        if not os.path.exists(input_sem):
            os.mkdir(input_sem)
        print("creating input_sem event directory")
        os.mkdir(event_sem)

    # remove forwardrun.sh script to signal that we are in a post forward state
    forward_run = os.path.join(drc["runfolder"], "forwardrun.sh")
    if os.path.exists(forward_run):
        print("removing forwardrun.sh")
        os.remove(forward_run) 

    print("post forward complete")


def build_adjoint(event_id):
    """
    Prepare the run folder for an adjoint simulation for a given event id
    :type event_id: str
    :param event_id: event id for creating appropriate files and directories
    """
    # recurring pathnames
    drc = dynamic_filenames()
    storage = os.path.join(drc["output_files"], "STORAGE", event_id)
    adjrun_template = os.path.join(
        drc["primer"], "simutils", "run_templates", "adjoint_simulation.sh")
    adjrun = os.path.join(drc["runfolder"], "adjointrun.sh")
    input_sem = os.path.join(drc["runfolder"], "INPUT_SEM")

    # make sure cmtsolution file is correct, e.g. if another forward run was
    # made between the prerequisite forward for this adjoint...
    print("symlinking CMTSOLUTION")
    cmtsolution = os.path.join(
        drc["primer"], "cmtsolution_files", "{}CMTSOLUTION".format(event_id)
    )
    cmtsolution_destination = os.path.join(drc["data"], "CMTSOLUTION")
    if os.path.exists(cmtsolution_destination):
        os.remove(cmtsolution_destination)
    os.symlink(cmtsolution, cmtsolution_destination)

    print("editing Par_file")
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="adjoint")

    def generate_adjointrun(fid_in, fid_out, event_id):
        """
        generate a forward run script to be called by sbatch
        """
        lines = open(fid_in, "r").readlines()
        for i, line in enumerate(lines):
            if "${EVENT_ID}" in line:
                lines[i] = line.replace("${EVENT_ID}", event_id+"_adj")
        with open(fid_out, "w") as f_out:
            f_out.writelines(lines)

    print("generating adjointrun file")
    if os.path.exists(adjrun):
        os.remove(adjrun)
    generate_adjointrun(adjrun_template, adjrun, event_id)

    print("checking OUTPUT_FILES")
    # check if forward saved files are symlinks and can be removed
    # if so, remove them and symlink from storage for fresh files
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        forward_check = glob.glob(os.path.join(drc["local_path"], quantity))
        if forward_check:
            if os.path.islink(forward_check[0]):
                print("removing symlinks of {}".format(quantity))
                for fid in forward_check:
                    os.remove(fid)
            else:
                sys.exit("{} files are real, please move".format(quantity))

        print("symlinking {} from storage to local path".format(quantity))
        files = glob.glob(os.path.join(storage, quantity))
        for fid in files:
            destination = os.path.join(drc["local_path"], os.path.basename(fid))
            os.symlink(fid, destination)

    # make sure the adjoint source file SEM/ is linked to the correct event
    # if one already exists, delete it as it should just be a symlink
    print("symlinking INPUT_SEM/{} to SEM".format(event_id))
    sem_folder = os.path.join(drc["runfolder"], "SEM")
    event_sem = os.path.join(input_sem, event_id)
    if os.path.exists(sem_folder) or os.path.islink(sem_folder):
        os.remove(sem_folder)
    os.symlink(event_sem, sem_folder)

    # symlink STATIONS_ADJOINT to the data file so Specfem knows where to look
    print("symlinking STATIONS_ADJOINT")
    stations_adjoint = os.path.join(drc["data"], "STATIONS_ADJOINT")
    input_stations_adjoint = os.path.join(event_sem, "STATIONS_ADJOINT")
    if os.path.exists(stations_adjoint) or os.path.islink(stations_adjoint):
        os.remove(stations_adjoint)
    os.symlink(input_stations_adjoint, stations_adjoint)

    print("adjoint build complete")


def post_adjoint():
    """
    After an adjoint simulation, files need to be moved around to get ready
    for more simulations, or for model update proceedings.
    """
    import shutil
    
    print("post adjoint")
    # recurring pathnames
    drc = dynamic_filenames()
    output = drc["output_files"]

    # check the output_solver.txt file to make sure simulation has finished
    if os.path.exists(output_solver):
        text = open(output_solver).read()
        if "End of the simulation" not in text:
            sys.exit("Simulation has not finished, disregarding")
    else:
        sys.exit("No output_solver.txt file, disregarding")

    # get event id from cmtsolution and set storage
    cmtsolution = os.path.join(output, "CMTSOLUTION")
    if not os.path.exists(cmtsolution):
        sys.exit("CMTSOLUTION doesn't exist in OUTPUT_FILES, need for event id")
    event_id = event_id_from_cmt(cmtsolution)
    event_storage = os.path.join(output, "STORAGE", event_id)

    # start moving files from OUTPUT_FILES/ to event storage
    for fid in glob.glob(os.path.join(output, "timestamp*")):
        dst = os.path.join(event_storage, os.path.basename(fid)+"_adjoint")
        shutil.move(fid, dst)

    # move output_solver to storage with new tag
    output_solver = os.path.join(output, "output_solver.txt")
    output_solver_new = os.path.join(event_storage, "output_solver_adjoint.txt")
    if os.path.exists(output_solver):
        print("moving output_solver.txt to storage")
        shutil.move(output_solver, output_solver_new)

    # remove sylinks to the proc*_save_forward_arrays.bin and
    # proc*_absorb_field.bin files that were symlinked by build_adjoint()
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        i = 0
        print("removing {} symlinks".format(quantity), end="...")
        for fid in glob.glob(os.path.join(drc["local_path"], quantity)):
            if os.path.islink(fid):
                os.remove(fid)
                i += 1
        print("{} symlinks removed".format(i))

    # move the output kernels of the adjoint simulation to storage
    print("moving kernels to storage", end="...")
    i = 0
    for fid in glob.glob(os.path.join(drc["local_path"], "*kernel.bin")):
        new_fid = os.path.join(event_storage, os.path.basename(fid))
        shutil.move(fid, new_fid)
        i += 1
    print("{} files moved".format(i))
    
    # remove adjointrun.sh script to signal a post adjoint state
    adjoint_run = os.path.join(drc["runfolder"], "adjointrun.sh")
    if os.path.exists(adjoint_run):
        print("removing adjoint runscript")
        os.remove(adjoint_run) 

    print("post adjoint complete")


def precondition_sum():
    """
    Summing kernels requires a few folders and symlinks to be set up beforehand
    To be run in the master folder
    """
    import shutil

    def symlink_kernels_to_master():
        """
        Because the work is split onto different run folders, we need to
        collect the input kernels for summation afterwards. Assumes that
        slave folders are located one directory up from the master
        """
        # set the master paths
        drc_master = directories()
        input_kernels = os.path.join(drc_master["runfolder"], "INPUT_KERNELS")

        # move into each run folder
        os.chdir("..")
        folders = glob.glob("*")
        folders.sort()
        event_ids = []
        for folder in folders:
            # check to make sure that we are actually looking at a run folder
            folder_check = os.path.exists(os.path.join(folder, "AUTHORS"))
            if os.path.isdir(folder) and folder_check:
                os.chdir(folder)
                # currently in a run folder
                drc = directories()
                storage = os.path.join(drc["output_files"], "STORAGE")

                # assume that only events we want to use are stored in storage
                # test cases and failed runs should be placed in cold storage
                events = glob.glob(os.path.join(storage, "*"))
                for event in events:
                    event_id = os.path.basename(event)
                    kernels = glob.glob(os.path.join(event, "*kernel.bin"))

                    # make directory in the INPUT_KERNELS master IFF there are
                    # kernels in the storage directory, AND doesn't already
                    # exist a folder in the INPUT_KERNELS directory
                    if bool(len(kernels)):
                        input_kernel_event = os.path.join(
                                                        input_kernels, event_id)
                        if not os.path.exists(input_kernel_event):
                            print("symlinking kernels from {dir}/{id}".format(
                                dir=os.path.basename(folder), id=event_id),
                                end="... ")
                            os.mkdir(input_kernel_event)
                            i = 0
                            for kernel in kernels:
                                kernel_sym = os.path.join(
                                    input_kernel_event, os.path.basename(kernel)
                                )
                                os.symlink(kernel, kernel_sym)
                                i += 1
                            print("{} kernels symlinked".format(i))
                            event_ids.append(event_id)
                        else:
                            print("dir exists for {}, skipping...".format(
                                event_id)
                            )
                # return to run folder holding directory
                os.chdir("..")

    # in the master directory, symlink from other run folder
    drc = directories()
    input_kernels = os.path.join(drc["runfolder"], "INPUT_KERNELS")
    kernels_list = os.path.join(drc["runfolder"], "kernels_list.txt")
    output_sum = os.path.join(drc["runfolder"], "OUTPUT_SUM")

    # flush the files from the INPUT_KERNEL directory, or make new
    if not os.path.exists(input_kernels):
        os.mkdir(input_kernels)
    else:
        # flush the existing input kernels,
        # ASSUMING THEY ARE FILLED WITH SYMLINKS.
        # DO NOT PUT REAL FILES IN THIS DIRECTORY
        input_kernel_events = glob.glob(os.path.join(input_kernels, "*"))
        for input_kernel_event in input_kernel_events:
            print("removing {}".format(input_kernel_event))
            shutil.rmtree(input_kernel_event)

    # make a new OUTPUT_SUM/ directory, or flush the existing
    if not os.path.exists(output_sum):
        os.mkdir(output_sum)
    else:
        # flush the existing output_sum
        output_sum_files = glob.glob(os.path.join(output_sum, "*kernel.bin"))
        print("removing files in OUTPUT_SUM/")
        for fid in output_sum_files:
            os.remove(fid) 
    
    print("symlinking kernels from other run folders")
    symlink_kernels_to_master() 

    # writing event numbers into kernels list while symlinking files to
    # INPUT_KERNELS/ directory
    with open(kernels_list, "r+") as f:
        print("writing kernels_list.txt")
        f.truncate(0)
        f.seek(0)
        for event in glob.glob(os.path.join(input_kernels, "*")):
            event_id = os.path.basename(event)
            f.write("{}\n".format(event_id))
    
    print("precondition sum complete")


def model_update():
    """
    Specfem's xmodel_update requires a few directories to be made in the 
    OUTPUT_FILES/ directory, and to have summed and smoothed kernels located
    in some input_kernels directory. Attenuation must also be turned off
    in the Par_file
    """
    drc = directories()
    output_sum = os.path.join(drc["runfolder"], "OUTPUT_SUM")
    model_update = os.path.join(
        drc["runfolder"], "src", "tomography", "model_update.f90")
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="update")

    def strip_markers(string):
        """hacky way to remove the markers from 'text/' 
        """
        return string[1:-2]

    # get model_update pathnames
    lines = open(model_update, "r").readlines()
    for line in lines:
        if ":: INPUT_KERNELS_DIR_NAME" in line:
            drctry = strip_markers(line.strip().split()[-1])
            input_kernels_dir = os.path.join(drc["output_files"], drctry)
        elif ":: LOCAL_PATH_NEW_NAME" in line:
            drctry = strip_markers(line.strip().split()[-1])
            local_path_new = os.path.join(drc["output_files"], drctry)
        elif ":: OUTPUT_STATISTICS_DIR_NAME" in line:
            drctry = strip_markers(line.strip().split()[-1])
            output_statistics_dir = os.path.join(drc["output_files"], drctry)
    
    # make the above named directories if they don't exit
    for drctry in [input_kernels_dir, local_path_new, output_statistics_dir]: 
        if not os.path.exists(drctry):
            os.mkdir(drctry)
   
    # flush symlinks in INPUT_KERNELS_DIR_NAME/
    old_kerns = glob.glob(os.path.join(input_kernels_dir, "*kernel_smooth.bin"))
    for ok in old_kerns:
        if os.path.islink(ok):
            os.remove(ok)   
 
    # symlink kernels into input_kernels_dir
    kernels = glob.glob(os.path.join(output_sum, "*kernel_smooth.bin"))
    if not kernels:
        sys.exit("Kernels must be smoothed before model update")
    for kernel in kernels:
        kernel_new = os.path.join(input_kernels_dir, os.path.basename(kernel))
        os.symlink(kernel, kernel_new)

    print("ready for model_update")
   

def post_model_update():
    """
    After a model update, we need to clean the run folder for the next forward
    simulations, and make sure the old outputs aren't overwritten by the next
    swath of simulations
    """ 
    import shutil
    drc = dynamic_filenames()
    slurm = os.path.join(drc["output_files"], "SLURM")    
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="update")
    
    # move slurm files into a separate folder
    print("moving slurm* files")
    if not os.path.exists(slurm):
        os.mkdir(slurm)
    slurmfiles = glob.glob(os.path.join(drc["runfolder"], "slurm-*.out"))
    for sfile in slurmfiles:
        sfile_new = os.path.join(slurm, os.path.basename(sfile))
        shutil.move(sfile, sfile_new)
   
    # change the name of OUTPUT_FILES, if M00_OUTPUT_FILES doesn't exist
    # if it exists, assume that this has already been run and continue
    old_output_files = os.path.join(drc["runfolder"], "M00_OUTPUT_FILES")
    if not os.path.exists(old_output_files):
        print("moving OUTPUT_FILES/ to M00_OUTPUT_FILES/")
        shutil.move(drc["output_files"], old_output_files)
        
        # set up new OUTPUT_FILES/
        os.mkdir(drc["output_files"])
        surface_h = os.path.join(old_output_files, "surface_from_mesher.h")
        values_h = os.path.join(old_output_files, "values_from_mesher.h")
        for fid in [surface_h, values_h]:    
            shutil.copyfile(fid, os.path.join(
                            drc["output_files"], os.path.basename(fid))
                           )
    else:
        query = input("M00_OUTPUT_FILES exists, "
                      "continue to populate new OUTPUT_FILES? [y/(n)]")
        if query != "y":
            sys.exit()
    
    # create new local_path in output_files
    print("creating new local_path")
    old_local_path = os.path.join(
        old_output_files, os.path.basename(drc["local_path"])
    )
    local_path = drc["local_path"]
    try:
        os.mkdir(local_path)
    except OSError:
        pass
    
    # mv relevant mesh files from old local_path to new local_path
    print("copying *attenuation, *Database, and q* files to new local_path")
    for tag in ["*attenuation.bin", "*Database", "*qkappa.bin", "*qmu.bin"]:
        fids = glob.glob(os.path.join(old_local_path, tag))
        for fid in fids:
            new_fid = os.path.join(local_path, os.path.basename(fid))
            shutil.copyfile(fid, new_fid)
    
    # mv vp_new.bin, vs_new.bin and rho_new.bin files to new local_path
    print("moving and renaming mesh files to new local_path")
    mesh_files = os.path.join(old_output_files, "mesh_files_m01")
    for tag in ["*vp_new.bin", "*vs_new.bin", "*rho_new.bin"]:
        fids = glob.glob(os.path.join(mesh_files, tag))
        for fid in fids:
            # get rid of the _new tag when renaming
            new_tag = os.path.basename(fid).split('_')
            new_tag = "{}_{}.bin".format(new_tag[0], new_tag[1])
            new_fid = os.path.join(local_path, new_tag)
            shutil.copyfile(fid, new_fid)
   
    print("ready for xgenerate_databases") 
      

def status_check():
    """
    For all directories in the current path, check what stage of the run cycle
    you're in by looking at the current state of outputs.

    The state machine is as follows:
    1) pre forward: build_forward has been run and the runfolder is ready for a 
                    a forward simulation
    2) run forward: runfolder is in the middle of a forward run
    3) finished forward: runfolder has finished the forward run but post_forward
                         has not yet been run so outputs are still able to be
                         overwritten by a new simulation
    4) post forward: post_forward has been run and outputs are stored away in
                     the storage folder. a corresponding adjoint run or a new
                     forward run can be pursued
    5) pre adjoint: build_adjoint has been run and the runfolder is ready for
                    an adjoint simulation 
    6) run adjoint: runfolder is in the middle of an adjoint run
    7) finished adjoint: same as finished forward. need to run post_adjoint
    8) post adjoint: post_adjoint has been run, runfolder is ready for new sims
    
    The state machine checks the following parameters:
    1) Par_file: to see what "SIMULATION_TYPE" and "SAVE_FORWARD" vars are
    2) timestamp*: to see if a simulation has been run but not a post script
    3) output_solver.txt: to see if "End of the simulation" has been printed
    4) forwardrun.sh/adjointrun.sh: to see if we are in pre or post simulation
    """ 
    def read_par_file(fid):
        """similar to edit par_file except only read in the choices and return
        a status pertaining to the par file
        :type fid: str
        :param fid: file id of Par_file
        :rtype: str
        :return: forward, adjoint or update
        """
        lines = open(fid, "r").readlines()

        # Change parameters
        for i, line in enumerate(lines):
            if "SIMULATION_TYPE" in line:
                simulation_type = line.strip().split()[-1]
            elif "SAVE_FORWARD" in line:
                save_forward = line.strip().split()[-1]
                break
        if simulation_type == "1":
            if save_forward == ".true.":
                status = "forward"
            elif save_forward == ".false.":
                status = "update"
        elif simulation_type == "3":
            status = "adjoint"
        else:
            status = None

        return status

    # pretty print some dividers and a title for status check
    title= "{:>15} {:>15} {:>18} {:>18}".format(
                                       "FOLDER", "EVENT ID", "STATUS", "QUEUE")
    print("{}\n{}\n{}".format("="*len(title), title, "="*len(title)))
    folders = glob.glob("*")
    folders.sort()

    # check statuses for each run folder in the parent directory
    for folder in folders:
        folder_check = os.path.exists(os.path.join(folder, "AUTHORS"))
        if os.path.isdir(folder) and folder_check:
            os.chdir(folder)

            # set up information to use for checks
            drc = directories()
            event_id = event_id_from_cmt(
                os.path.join(drc["data"], "CMTSOLUTION")
            )
            gen_status = read_par_file(os.path.join(drc["data"], "Par_file"))

            # check files to show where the status is
            output_solver = os.path.join(
                drc["output_files"], "output_solver.txt")
            output_solver_bool = os.path.exists(output_solver)
            timestamps = os.path.join(drc["output_files"], "timestamp*")
            timestamp_bool = bool(len(glob.glob(timestamps)))

            # check if we are in pre, running or post based on output_solver
            # and the existence of timestamp files
            if output_solver_bool and timestamp_bool:
                text = open(output_solver, "r").read()
                # decide between running and finished by checking if sim ended
                if "End of the simulation" in text:
                    status = "finished {}".format(gen_status)
                    queue = "run post_{}".format(gen_status)
                else:
                    status = "running {}".format(gen_status)
                    # check percentage complete
                    lines = text.split("\n")
                    for line in reversed(lines):
                        if "We have done" in line:
                            percentage = line.split()[3]
                            break
                    queue = "{:.2f}%".format(float(percentage))
            else:
                # differentiate pre and post simulation by presence of runscript
                if os.path.exists(os.path.join(
                           drc["runfolder"], "{}run.sh".format(gen_status))):
                    status = "pre {}".format(gen_status)
                    queue = "run {}".format(gen_status)
                else:
                    # change the queue message depending on post fwd, adj or upd
                    if gen_status == "forward":
                        queue = "build adjoint" 
                    elif gen_status == "adjoint":
                        queue = "open"
                    elif gen_status == "update":
                        queue = "-"
                    status = "post {}".format(gen_status)

            # print status per run folder
            print("{rf:>15} {id:>15} {sts:>18} {qu:>18}".format(
                                rf=folder, id=event_id, sts=status, qu=queue))
            os.chdir("..")
    print("="*len(title))


if __name__ == "__main__":
    # simple argument distribution
    available_funcs = ["build_forward", "post_forward", "build_adjoint",
                       "post_adjoint", "precondition_sum", "model_update",
                       "post_model_update", "status_check"]
    try:
        func = sys.argv[1]
    except IndexError:
        sys.exit("argument 1 must be in\n {}".format(available_funcs))

    if func not in available_funcs:
        sys.exit("argument 1 must be in\n {}".format(available_funcs))

    try:
        event_id = sys.argv[2]
    except IndexError:
        if "build" in func:
            sys.exit("argument 2 must be event_id")

    if func == "build_forward":
        build_forward(event_id)
    elif func == "build_adjoint":
        build_adjoint(event_id)
    elif func == "post_forward":
        post_forward()
    elif func == "post_adjoint":
        post_adjoint()
    elif func == "precondition_sum":
        precondition_sum()
    elif func == "model_update":
        model_update()
    elif func == "post_model_update":
        post_model_update()
    elif func == "status_check":
        status_check()

