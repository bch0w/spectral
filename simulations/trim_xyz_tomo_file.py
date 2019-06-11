import os
import numpy as np


def xyz_reader(filepath):
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
    
    if os.path.exists(filepath + ".npy") and os.path.exists(filepath + ".npz"):
        header = np.load(filepath + ".npz")
        data = np.load(filepath + ".npy")
        return header, data
    else:
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

            header = {"orig_x": orig_x, "orig_y": orig_y, "orig_z": orig_z,
                      "end_x": end_x, "end_y": end_y, "end_z": end_z,
                      "spacing_x": spacing_x, "spacing_y": spacing_y,
                      "spacing_z": spacing_z, "nx": nx, "ny": ny, "nz": nz,
                      "vp_min": vp_min, "vp_max": vp_max, "vs_min": vs_min,
                      "vs_max": vs_max, "rho_min": rho_min, "rho_max": rho_max,
                     }
            
            data = np.array([_.strip().split() for _ in lines[4:]], dtype=float)
            np.save(filepath, data)  # save data in .npy
            np.savez(filepath, **header)  # save header in .npz

        return header, data


def trim_xyz_file(data, bounds):
    """
    sometime the xyz file is too large, trim it down to new dimensions
    :param header:
    :param data:
    :return:
    """
    # convert data into numpy array for easier working
    for i, bound in enumerate(bounds):
        too_small = np.where(data[:, i] < bound[0])[0]
        too_large = np.where(data[:, i] > bound[1])[0]
    
    to_remove = np.unique(np.concatenate((too_small, too_large), 0))
    new_data = np.delete(data, to_remove, 0)
    
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
    
    return new_parsed_header, new_data


def write_new_xyz(header, data, fidout_template):
    """
    write out a new xyz file with proper header and data
    """
    fid = fidout_template.format(
            x=(header["end_x"] - header["orig_x"])*1E-3,
            y=(header["end_y"] - header["orig_y"])*1E-3,
            )

    with open(fid, "w") as f:
        # write header
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
            header["orig_x"], header["orig_y"], header["orig_z"],
            header["end_x"], header["end_y"], header["end_z"])
        )
        f.write("{:.1f} {:.1f} {:.1f}\n".format(
            header["spacing_x"], header["spacing_y"], header["spacing_z"])
                )
        f.write("{:d} {:d} {:d}\n".format(
            header["nx"], header["ny"], header["nz"])
        )
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
            header["vp_min"], header["vp_max"], header["vs_min"], 
            header["vs_max"], header["rho_min"], header["rho_max"])
        )
        # write data by line
        for line in data:
            x, y, z, vp, vs, rho, qp, qs = line.tolist()
            f.write(
             "{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
                   x, y, z, vp, vs, rho, qp, qs)
            )


if __name__ == "__main__":
    # path = "/Users/chowbr/Documents/subduction/data/KUPEDATA/tomo_files"
    path = "/seis/prj/fwi/bchow/data/KUPEDATA/tomo_files"
    name_template = "nz_full_eberhart2015_{}.xyz"
    fullpath = os.path.join(path, name_template)
   
    # output file names
    out_template = "nz_x{x}_y{y}_eberhart2015_{section}.xyz"
    outpath = os.path.join(path, out_template)
    
    # different names of the tomo files
    sections = ["shallow", "crust", "mantle"]

    # hardcoded, defined by MESH_SRTM30P_139_162_4000m [min, max]
    x_bounds = [167500.0, 708000.0]
    y_bounds = [5270000.0, 5915000.0]
    z_bounds = [-4000000., 2100.0]
    bounds = [x_bounds, y_bounds, z_bounds]

    # trim xyz files based on given bounds
    for i, section in enumerate(sections):  #, crust_data, mantle_data]:
        header, data = xyz_reader(fullpath.format(section))
        new_header, new_data = trim_xyz_file(data, bounds)
        outpath_temp = outpath.format(x="{x:.0f}", y="{y:.0f}", section=section)
        write_new_xyz(new_header, new_data, outpath_temp) 
        




