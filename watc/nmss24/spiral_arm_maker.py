"""
Following Kennet et al. 2015, determine geographic locations for a spiral arm
array for station deployment planning

A:   radial aperture of the array
n_a: number of arms
n_r: number of station rings 
phi: angular span of each arm
k:   arm number
d0:  rotation constant
"""
import numpy as np
import matplotlib.pyplot as plt


def get_initial_angle(rotation_constant=0, arm_number=0, number_spiral_arms=1):
    """Define the initial angle of the kth spiral arm"""
    return rotation_constant + np.rad2deg(2 * np.pi * (arm_number / 
                                                       number_spiral_arms))


def get_ring_radius(radial_aperture_km=10, ring_number=1, 
                    number_station_rings=4):
    """Determine the ring radii for all rings in the array"""
    return radial_aperture_km * (ring_number / number_station_rings)


def get_station_azimuth(initial_angle=0, angular_arm_span_deg=120, 
                        ring_number=1, number_station_rings=4):
    """Calculate the azimuth for the station on the jth ring of the kth arm"""
    return initial_angle + angular_arm_span_deg * (ring_number / 
                                                   number_station_rings)


# Set parameters for spiral design
radial_aperture_km = 20
number_spiral_arms = 3
number_station_rings = 4
angular_arm_span_deg = 120
rotation_constant_deg = 0
include_central_station = False
central_station_location = (0, 0)

f, ax = plt.subplots()

# Calculate parameters for the station on the jth ring and kth arm
if include_central_station:
    start_idx = 0
else:
    start_idx = 1

nsta = 0
for arm_number in range(1, number_spiral_arms + 1):
    print(f"ARM NUMBER {arm_number}")
    initial_angle_deg = get_initial_angle(rotation_constant_deg, arm_number,
                                          number_spiral_arms)
    for ring_number in range(start_idx, number_station_rings + 1):
        ring_radius_km = get_ring_radius(radial_aperture_km, ring_number,
                                         number_station_rings)
        # Plot a circle at the given ring radius
        circ = plt.Circle((0, 0), ring_radius_km, color="k", fill=False, 
                          alpha=0.5, zorder=4)        
        ax.add_patch(circ)
        plt.text(s=f"{ring_radius_km}km", x=0., y=ring_radius_km)

        station_azimuth_deg = get_station_azimuth(initial_angle_deg, 
                                                  angular_arm_span_deg,
                                                  ring_number, 
                                                  number_station_rings) % 360
       
        print(f"(rad, az): {ring_radius_km}, {station_azimuth_deg:.2f}")

        # Convert the radius and azimuth to geographic coordinates
        xdist = ring_radius_km * np.sin(np.deg2rad(station_azimuth_deg))
        ydist = ring_radius_km * np.cos(np.deg2rad(station_azimuth_deg))

        # Calculate distance from central location
        sta_x = xdist + central_station_location[0]
        sta_y = ydist + central_station_location[1]

        plt.scatter(sta_x, sta_y, c="r", marker="v", ec="k", s=60, zorder=5)
        plt.text(s=f"{station_azimuth_deg:.2f}", x=sta_x, y=sta_y)
        nsta += 1

print(f"N_STA={nsta}")
ax.set_aspect("equal")
plt.show()




