import os
import numpy as np
from scipy.interpolate import griddata

from checkerboardiphy import xyz_reader
from pyatoa.utils.operations.source_receiver import lonlat_utm


def convert_interp(latlon_fid, xyz_npy_fid, plot=False):
    """
    Matlab spat out some converted lat lon coordinates but they need to be in 
    UTM60 and also need to be interpolated because sampling in Carl's
    custom coordinate system is not the same as sampling in lat/lon or utm
    """
    latlon = np.loadtxt(latlon_fid)
    lat = latlon[:, 0]
    lon = latlon[:, 1]
    x_irreg, y_irreg = lonlat_utm(lon, lat)
    
    # Read in the data file to replace the coordinates
    data = np.load(xyz_npy_fid)

    # These are taken from converting my chosen map coordinates to UTM-60 and
    # putting a small buffer around the dimensions
    # lat_min, lon_min = 172.9998, -42.5007;
    # lat_max, lon_max = 179.5077, -36.9488
    x_min, x_max = 170000., 725000.
    y_min, y_max = 5280000., 5910000.

    # For the coastline plot
    if plot:
        coastline = np.load("./nz_resf_coast_mod_utm60H_xyz.npy")
        coastline = coastline[
            np.where((coastline[:, 0] > x_min) &
                     (coastline[:, 0] < x_max) &
                     (coastline[:, 1] > y_min) &
                     (coastline[:, 1] < y_max)
                     )[0]]

    # Create the regular grid to be interpolated onto
    spacing = 1000.
    x_reg = np.arange(x_min, x_max, spacing)
    y_reg = np.arange(y_min, y_max, spacing)

    x_grid, y_grid = np.meshgrid(x_reg, y_reg)
    x_out = x_grid.flatten()
    y_out = y_grid.flatten()
    
    # Interpolate for each depth slice
    list_interp = []
    data_out = None
    for z_value in np.unique(data[:, 2]):
        data_ind = np.where(data[:, 2] == z_value)[0]

        # Loop over the values to be interpolated
        for val in range(3, 8):
            interp_vals = griddata(
                points=(x_irreg[data_ind], y_irreg[data_ind]),
                values=data[data_ind, val],
                xi=(x_grid, y_grid)
            )
            list_interp.append(interp_vals.flatten())

            if plot:
                plt.contourf(x_grid, y_grid, interp_vals)
                plt.scatter(coastline[:, 0], coastline[:, 1], c='k',
                            marker='.', s=1.5
                            )
                plt.title("{fid} {z}".format(fid=xyz_npy_fid, z=z_value))
                plt.xlabel("Easting (m)")
                plt.ylabel("Northing (m)")
                plt.show()

        # Create the ndarray by creating column vectors and mushing em together
        z_out = np.ones(len(x_grid.flatten())) * z_value
        vp_out = list_interp[0].flatten()
        vs_out = list_interp[1].flatten()
        rho_out = list_interp[2].flatten()
        qp_out = list_interp[3].flatten()
        qs_out = list_interp[4].flatten()

        depth_slice = np.column_stack(
            (x_out, y_out, z_out, vp_out, vs_out, rho_out, qp_out, qs_out)
        )
        if data_out is None:
            data_out = depth_slice
        data_out = np.concatenate((data_out, depth_slice))

    np.save('utm60_{}'.format(xyz_npy_fid), data_out)
    return data_out
convert_interp('./matlab_convert_carl_cust_to_latlon/crust_latlon.txt', 'nz_full_eberhart2015_crust.xyz.npy')


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
    
    # delete in place
    data = np.delete(data, to_remove, 0)
    
    # make new header
    parsed_header = {"orig_x": data[:, 0].min(), "orig_y": data[:, 1].min(),
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

    head_len = parsed_header["nx"] * parsed_header["ny"] * parsed_header["nz"] 
    assert head_len == len(data)

    return parsed_header, data


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
    # path = "/seis/prj/fwi/bchow/data/KUPEDATA/tomo_files"
    path = "./"
    name_template = "nz_north_eberhart2015_{}.xyz"
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
    for i, section in enumerate(sections):
        print("trimming section {}".format(section))
        header, data = xyz_reader(fullpath.format(section))
        new_header, new_data = trim_xyz_file(data, bounds)
        outpath_temp = outpath.format(x="{x:.0f}", y="{y:.0f}", section=section)
        write_new_xyz(new_header, new_data, outpath_temp)
        




