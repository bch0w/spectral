import os
import numpy as np


def xyz_reader(filepath, fidout=''):
    """
    read .xyz file used as input for specfem and return information
    header parsing specifically related to specfem
    data is delineated as :
    x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs
    -origin and endpoints are the coordinates in meters
    -spacing given in units of meters for UTM projections
    -nx, ny, nz = [(end_x - orig_x)/spacing_x] + 1
    :param filepath:
    :return:
    """
    with open(filepath) as f:
        lines = f.readlines()

    # parse the header and make sure the values are floats
    orig_x, orig_y, orig_z, end_x, end_y, end_z = \
        [float(_) for _ in lines[0].strip().split()]
    spacing_x, spacing_y, spacing_z = \
        [float(_) for _ in lines[1].strip().split()]
    nx, ny, nz = [float(_) for _ in lines[2].strip().split()]
    vp_min, vp_max, vs_min, vs_max, rho_min, rho_max = \
        [float(_) for _ in lines[3].strip().split()]

    parsed_header = {"orig_x": orig_x, "orig_y": orig_y, "orig_z": orig_z,
                     "end_x": end_x, "end_y": end_y, "end_z": end_z,
                     "spacing_x": spacing_x, "spacing_y": spacing_y,
                     "spacing_z": spacing_z, "nx": nx, "ny": ny, "nz": nz,
                     "vp_min": vp_min, "vp_max": vp_max, "vs_min": vs_min,
                     "vs_max": vs_max, "rho_min": rho_min, "rho_max": rho_max,
                      }

    if fidout:
        data = np.array([_.strip().split() for _ in lines[4:]])
        import ipdb;ipdb.set_trace()
        np.save(fidout, data)

    return lines[:4], lines[4:], parsed_header


def trim_xyz_file(data, bounds):
    """
    sometime the xyz file is too large, trim it down to new dimensions
    :param header:
    :param data:
    :return:
    """
    # convert data into numpy array for easier working
    import ipdb;ipdb.set_trace()
    data = np.array([_.strip().split() for _ in data])
    for i, bound in enumerate(bounds):
        too_small = np.where(data[:, i] < bound[0])[0]
        too_large = np.where(data[:, i] > bound[1])[0]
    
    import ipdb;ipdb.set_trace()
    to_remove = np.unique(np.concatenate((too_small, too_large), 0))
    new_data = np.delete(data, to_remove)

    # make new header
    new_parsed_header = {"orig_x": data[:, 0].min(), "orig_y": data[:, 1].min(),
                         "orig_z": data[:, 2].min(), "end_x": data[:, 0].max(),
                         "end_y": data[:, 1].max(), "end_z": data[:, 2].max(),
                         "spacing_x": 2000., "spacing_y": 2000.,
                         "spacing_z": 1000.,
                         "nx": len(np.unique(data[:, 0])),
                         "ny": len(np.unique(data[:, 1])),
                         "nz": len(np.unique(data[:, 2])),
                         "vp_min": data[:, 3].min(), "vp_max": data[:, 3].max(),
                         "vs_min": data[:, 4].min(), "vs_max": data[:, 4].max(),
                         "rho_min": data[:, 5].max(),
                         "rho_max": data[:, 5].max(),
                      }
    
    import ipdb;ipdb.set_trace()
    return new_parsed_header, new_data


def write_new_xyz(head, data, fid):
    """
    write out a new xyz file with proper header and data
    """
    with open(fid, "w") as f:
        # write header
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}".format(
            head["orig_x"], head["orig_y"], head["orig_z"],
            head["end_x"], head["end_y"], head["end_z"])
        )
        f.write("{:.1f} {:.1f} {:.1f}".format(
            head["spacing_x"], head["spacing_y"], head["spacing_z"])
                )
        f.write("{:.1f} {:.1f} {:.1f}".format(
            head["nx"], head["ny"], head["nz"])
        )
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}".format(
            head["vp_min"], head["vp_max"], head["vs_min"], head["vs_max"],
            head["rho_min"], head["rho_max"])
        )
        # write data by line
        for line in data:
            x, y, z, vp, vs, rho, qp, qs = line.split()
            f.write(
               "{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}".format(
                   x, y, z, vp, vs, rho, qp, qs)
            )


if __name__ == "__main__":
    path = "/Users/chowbr/Documents/subduction/data/KUPEDATA/tomo_files"
    name_template = "nz_full_eberhart2015_{}.xyz"
    fullpath = os.path.join(path, name_template)

    shallow_header, shallow_data, shallow_parsed = xyz_reader(
        fullpath.format("shallow"), fullpath.format("shallow"))
    # crust_header, crust_data, crust_parsed = xyz_reader(
    #     fullpath.format("crust"))
    # mantle_header, mantle_data, mantle_parsed = xyz_reader(
    #     fullpath.format("mantle"))

    # hardcoded, defined by MESH_SRTM30P_139_162_4000m
    x_bounds = [167500.0, 70800.0]
    y_bounds = [5270000.0, 5915000.0]
    z_bounds = [-4000000, 2100.0]
    bounds = [x_bounds, y_bounds, z_bounds]

    # trim xyz files based on given bounds
    for i, data in enumerate([shallow_data]):  #, crust_data, mantle_data]:
        new_header, new_data = trim_xyz_file(data, bounds[i])

        




