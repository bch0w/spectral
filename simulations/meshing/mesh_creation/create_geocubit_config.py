"""short script to create GEOCUBIT config files for mesh generation, takes
informatino from topo and moho files and feeds them into a template.
script assumes that file naming remains consistent between the folder,
the topo files, moho files and the output config file
"""
import os
import glob

def set_template():
    """to keep the large format template separate from the function
    """
    config_template = """[cubit.options]
cubit_info                      = on
echo_info                       = on
jou_info                        = off
jer_info                        = off
working_dir                     = WORKINGDIR
output_dir                      = OUTPUTDIR
save_geometry_cubit             = False
save_surface_cubit              = False
export_exodus_mesh              = True
monitored_cpu                   = 0
localdir_is_globaldir           = False
parallel_import                 = False
#scratchdir                     = None

[simulation.cpu_parameters]
number_processor_xi             = 9
number_processor_eta            = 2

[geometry.volumes]
volume_type                     = layercake_volume_ascii_regulargrid_regularmap
longitude_min                   = {LONMIN}
longitude_max                   = {LONMAX}
latitude_min                    = {LATMIN}
latitude_max                    = {LATMAX}
nx                              = {NX}
ny                              = {NY}
unit                            = utm

[geometry.volumes.layercake]
nz                              = 3
bottomflat                      = True
depth_bottom                    = -400000
geometry_format                 = ascii
filename                        = {MOHO},{TOPO},

[meshing]
map_meshing_type                = regularmap
iv_interval                     = 36,2
size                            = 12600
or_mesh_scheme                  = map
ntripl                          = 2
smoothing                       = False
coarsening_top_layer            = False
refinement_depth                = 1,1
    """

    return config_template

def fill_template(MESH_TAG):
    """fill er up
    """
    os.chdir(MESH_TAG)
    cwd = os.getcwd()
    topofile = glob.glob('topo*')[0]
    mohofile = glob.glob('moho*')[0]
    if not topofile:
        import sys
        sys.exit('FOLDER EMPTY!')

    # LAT LON MIN MAX
    with open(topofile,'r') as f:
        lines = f.readlines()
    firstline = lines[0].split()
    lastline = lines[-1].split()

    lonmin,latmin,_ = firstline
    lonmax,latmax,_ = lastline

    # NX AND NY
    meshtagsplit = MESH_TAG.split('_')
    _,_,nx,ny,spacing = meshtagsplit

    # TOPO AND MOHO FILES
    topofullpath = os.path.join(cwd,topofile)
    mohofullpath = os.path.join(cwd,mohofile)


    config_template = set_template()
    config_template = config_template.format(LONMIN=lonmin,
                                             LATMIN=latmin,
                                             LONMAX=lonmax,
                                             LATMAX=latmax,
                                             NX=nx,
                                             NY=ny,
                                             MOHO=mohofullpath,
                                             TOPO=topofullpath
                                             )
    # SAVE IT
    meshfilename = MESH_TAG + ".cfg"
    with open(meshfilename,'w') as f:
        f.write(config_template)


if __name__ == "__main__":
    # MESH_TAG = "nz_SRTM15P_553_622_1000m"
    MESH_TAG = input('Input folder name: ')
    if MESH_TAG[-1] == '/':
        MESH_TAG = MESH_TAG[:-1]
    print("Creating config for",MESH_TAG)
    fill_template(MESH_TAG)
