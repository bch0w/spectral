"""
A script to generate the necessary files for Specfem3D Cartesians internal
mesher, Meshfem3D.
"""
import os
import numpy as np
from pyatoa.utils.operations.source_receiver import lonlat_utm
from pyatoa.utils.operations.calculations import myround


def set_parameters():
    """
    Define the necessary parameters here, these will be accessed by the
    constraint fuctions throughout the script
    :return:
    """
    parameters = {"lat_min": -43.,
                  "lat_max": -37.,
                  "lon_min": 173.,
                  "lon_max": 179.,
                  "utm_projection": -60,
                  "mesh_depth_km": 400.,
                  "ndoublings": 2,
                  "interfaces": [8, 36],
                  "interface_fids": ["interface_crustshallow.dat",
                                     "interface_moho.dat"],
                  "nproc": 160,
                  "shortest_period_s": 10,
                  "vs_min_km_per_s": 1.,
                  }

    assert(parameters["ndoublings"] == len(parameters["interfaces"]))

    # Convert latlon to UTM
    x_min, y_min = lonlat_utm(parameters["lon_min"], parameters["lat_min"],
                              parameters["utm_projection"])
    x_max, y_max = lonlat_utm(parameters["lon_max"], parameters["lat_max"],
                              parameters["utm_projection"])

    # Set the UTM coordinates in the parameters
    parameters["x_min"] = x_min
    parameters["y_min"] = y_min
    parameters["x_max"] = x_max
    parameters["y_max"] = y_max
    parameters["x_length_km"] = (x_max - x_min) * 1E-3
    parameters["y_length_km"] = (y_max - y_min) * 1E-3

    return parameters


def minimum_grid_spacing(slowest_wavespeed, shortest_period):
    """
    Define minimum grid spacing based on slowest wavespeed and shortest period
    :return:
    """
    # Value of 2 comes from two points per wavelength
    min_grid_space = (shortest_period * slowest_wavespeed) / 2
    print(f"\tminimum grid spacing calculated as {min_grid_space}")
    grid_space = myround(min_grid_space, 2, 'down')
    print(f"\t\trounded down to {grid_space}")

    return grid_space


def number_of_processors(nproc, x_length, y_length):
    """
    Set the number of processors base on the ratio of mesh length
    :param nproc:
    :return:
    """
    # Determine the ratio
    ratio = min(x_length, y_length) / max(x_length, y_length)

    # Define which direction is shorter
    if x_length < y_length:
        short_direction = "x"
    else:
        short_direction = "y"

    # Start guessing processor ratios at the square root, until 2 integers found
    guess_a = round(np.sqrt(ratio * nproc))
    while True:
        guess_b = nproc / guess_a
        if guess_b.is_integer():
            # Assign the short direction correctly
            if short_direction == "x":
                nproc_x = min(guess_a, guess_b)
                nproc_y = max(guess_a, guess_b)
            else:
                nproc_x = max(guess_a, guess_b)
                nproc_y = min(guess_a, guess_b)
            print(f"\tNPROC_X = {nproc_x}\n\tNPROC_Y = {nproc_y}")
            return int(nproc_x), int(nproc_y)
        else:
            guess_a += 1


def number_of_elements(nproc_x, nproc_y, x_length, y_length, grid_space,
                       shortest_period_s):
    """
    Define the number of elements in each horizontal direction based on number
    of processors
    :return: 
    """
    # meshfem requires the number of grid points be an integer multiple of 8
    # times the number of processors in a given direction
    nex_x = myround(x_length / grid_space, nproc_x * 8, "near")
    nex_y = myround(y_length / grid_space, nproc_y * 8, "near")

    # ensure that the short direction is maintained
    while (nproc_x < nproc_y) and (nex_x > nex_y):
        nex_y += nproc_y * 8

    # Minimum element number set by Eq. 3.1 in Specfem manual
    assert(nex_x >= 2 * 288 / shortest_period_s)

    dx =  x_length / nex_x
    dy =  y_length / nex_y

    print(f"\tNEX_X = {nex_x}\n\tNEX_Y = {nex_y}\n"
          f"\tdx = {dx:.2f}km\n\tdy = {dy:.2f}km"
          )

    return int(nex_x), int(nex_y)


def vertical_doubling_proportions(depth, ndoublings, interfaces, grid_space):
    """
    Set the number of vertical doubling layers
    :param depth:
    :param ndoublings:
    :param interfaces:
    :return:
    """
    # Start interfaces from the top, include the bottom interface of depth
    all_interfaces = [0] + interfaces + [depth]
    all_interfaces.sort()

    layers = []
    for i in range(ndoublings + 1):
        j = i + 1
        num_layers = ((all_interfaces[j] - all_interfaces[i]) /
                      ((2 ** i) * grid_space))
        layers.append(myround(num_layers, 1, "near"))

    # Get an even number of elements, place new layers on top
    while sum(layers) < myround(sum(layers), 2, "up"):
        layers[0] += 1

    # Start counting layers from bottom
    layers.sort(reverse=True)

    print("\t{nlay} sections; bottom to top {layers} elements per layer".format(
        nlay=len(layers), layers=layers))

    return layers


def write_mesh_par_file(template, lat_min, lat_max, lon_min, lon_max,
                       depth, utm_projection, nex_x, nex_y, nproc_x, nproc_y,
                       ndoublings, layers):
    """
    Write the Mesh_Par_file with the given values and a template script
    :return:
    """
    with open(template, "r") as f:
        lines = f.read()

    lines = lines.format(lat_min=lat_min, lat_max=lat_max, lon_min=lon_min,
                         lon_max=lon_max, depth=depth,
                         utm_projection=utm_projection,
                         nex_xi=nex_x, nex_eta=nex_y, nproc_xi=nproc_x,
                         nproc_eta=nproc_y, ndoublings=ndoublings,
                         nz_doubling_1=sum(layers[:1]),
                         nz_doubling_2=sum(layers[:2]),
                         nmaterials=ndoublings+1, nregions=ndoublings+1,
                         nz_1a=1,
                         nz_1b=sum(layers[:1]),
                         nz_2a=sum(layers[:1]) + 1,
                         nz_2b=sum(layers[:2]),
                         nz_3a=sum(layers[:2]) + 1,
                         nz_3b=sum(layers[:3]),
                         )

    # Put the outputs into a directory for easy transfer
    base = './meshfem3D_files'
    if not os.path.exists(base):
        os.makedirs(base)

    print("\twriting Mesh_Par_file")
    with open(os.path.join(base, "Mesh_Par_file"), "w") as f:
        f.write(lines)


def write_interfaces(template, layers, interfaces, lat_min, lon_min, fids=[]):
    """
    Write the interfaces.dat file as well as the corresponding flat interface
    layers. Topo will need to be written manually
    :param template:
    :param layers:
    :param lat_min:
    :param lon_min:
    :return:
    """
    # Set the names of the interface files
    if not fids:
        fids = []
        for i in range(1, len(layers)):
            fids.append("interface_{}.dat".format(i))

    with open(template, "r") as f:
        lines = f.read()

    # Format the file, interface numbering starts from the bottom, i.e.
    # interface1 is the first from the bottom
    lines = lines.format(number_interfaces=len(layers), lat_min=lat_min,
                         lon_min=lon_min, interface1_fid=fids[-1],
                         interface2_fid=fids[-2], nz_1=layers[0],
                         nz_2=layers[1], nz_3=layers[2]
                         )


    # Put the outputs into a directory for easy transfer
    base = './meshfem3D_files'
    if not os.path.exists(base):
        os.makedirs(base)

    # Write to a new file
    print("\twriting interfaces.dat")
    with open(os.path.join(base, "interfaces.dat"), "w") as f:
        f.write(lines)

    # Write the individual interface files
    for fid, interface in zip(fids, interfaces):
        # Skip top interface (topography)
        print(f"\twriting {fid}")
        with open(os.path.join(base, fid), "w") as f:
            for i in range(5):
                f.write("-{}\n".format(abs(int(interface * 1E3))))


def prepare_meshfem():
    """
    Run all the above functions and output the necessary meshfem3D files
    :return:
    """
    mesh_par_file_template = "./template_Mesh_Par_file"
    interfaces_template = "./template_interfaces.dat"

    print("Preparing Meshfem3D files")
    pars = set_parameters()
    grid_space = minimum_grid_spacing(slowest_wavespeed=pars["vs_min_km_per_s"],
                                      shortest_period=pars["shortest_period_s"]
                                      )
    nproc_x, nproc_y = number_of_processors(nproc=pars["nproc"],
                                            x_length=pars["x_length_km"],
                                            y_length=pars["y_length_km"]
                                            )
    nex_x, nex_y = number_of_elements(
        nproc_x=nproc_x, nproc_y=nproc_y, x_length=pars["x_length_km"],
        y_length=pars["y_length_km"], grid_space=grid_space,
        shortest_period_s=pars["shortest_period_s"]
    )
    layers = vertical_doubling_proportions(depth=pars["mesh_depth_km"],
                                           ndoublings=pars["ndoublings"],
                                           interfaces=pars["interfaces"],
                                           grid_space=grid_space
                                           )
    write_mesh_par_file(template=mesh_par_file_template,
                        lat_min=pars["lat_min"], lat_max=pars["lat_max"],
                        lon_min=pars["lon_min"], lon_max=pars["lon_max"],
                        depth=pars["mesh_depth_km"],
                        utm_projection=pars["utm_projection"],
                        nex_x=nex_x, nex_y=nex_y, nproc_x=nproc_x,
                        nproc_y=nproc_y, ndoublings=pars["ndoublings"],
                        layers=layers
                        )
    write_interfaces(template=interfaces_template, layers=layers,
                     interfaces=pars["interfaces"], lat_min=pars["lat_min"],
                     lon_min=pars["lon_min"], fids=pars["interface_fids"]
                     )


if __name__ == "__main__":
    prepare_meshfem()
