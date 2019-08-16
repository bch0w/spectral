
"""
There are some values required by the internal mesher that control the sizes
of spectral elements in the horizontal and vertical directions. This collection
of functions should help determine those values in a standardized manner
"""
from pyatoa.utils.operations.source_receiver import lonlat_utm


def provide_mesh_par_info():
    """
    Determine the number of vertical layers for a given number of doublings such
    that the top layer contains roughly cube shaped elements
    """
    # User Parameters
    lat_min = -43.
    lat_max = -37.
    lon_min = 173.
    lon_max = 179.
    utm_zone = -60
    shortest_period_s = 10.

    # Find the mesh dimensions in units of km
    x_min, y_min = lonlat_utm(lon_min, lat_min, utm_zone)
    x_max, y_max = lonlat_utm(lon_max, lat_max, utm_zone)
    x_length_km = (x_max - x_min) * 1E-3
    y_length_km = (y_max - y_min) * 1E-3

    ratio_x_to_y = x_length_km / y_length_km
    print(f"The ratio of x to y is {ratio_x_to_y}")

    # Tell the user the smallest value of NEX for shortest period
    # Equation taken from the Specfem3D manual
    nex_min = 288 * 2 / shortest_period_s
    print(f"Smallest number of elements for {shortest_period_s}s is {nex_min}")

    return x_length_km, y_length_km


def vertical_layers():
    """
    Determine the minimum number of vertical layers to keep the target 
    resolution
    """
    # User Parameters
    depth_block_km = 350.
    nex_xi = 72
    nex_eta = 96
    ndoublings = 2
    x_length_km = 504
    y_length_km = 672
    
    # Determine element size in the horizontal directions
    xi_element_length = x_length_km / nex_xi
    eta_element_length = y_length_km / nex_eta

    target_element_size = min(xi_element_length, eta_element_length)
    print(f"Target vertical element length {target_element_size}km")


def create_interface():
    """
    The internal mesher requires files containing interface elevations to 
    produce e.g. topography on the mesh. The files need to be in y.x increasing
    order, this function will take a topography file and parse it into the 
    required file order
    """




def calculate_doubling_layer():
    """
    Just determine the layer number that corresponds to a depth 
    """
    # User parameters
    depth_block_km = 400.
    ndoublings = 2
    number_layers = [7,8,40]
    desired_depth_km = 33.

    assert(ndoublingsi + 1 == len(number_layers))
    
    # Calculate the layer number

    smallest_block_km =




