"""read in active fault and convert to npz arrays
"""


def offshore_faults():
    with open("./new_zealand_offshore_faults.gmt", 'r') as f:
        lines = f.readlines()

    # 550_641 mesh
    map_corners = [-42.5007, -36.9488, 172.9998, 179.5077]

    faults, lats, lons = [], [], []
    i = 0
    for line in lines[9:]:
        if line == '>\n':
            i += 1
            continue
        elif line[0] == '#':
            continue
            # try:
            #     fault = line.split('"')[1]
            # except IndexError:
            #     fault = line.split('|')[1]
        else:
            lon, lat = line.strip().split(' ')
            lat = float(lat)
            lon = float(lon)
            if (map_corners[0] < lat < map_corners[1]) and \
                    (map_corners[2] < lon < map_corners[3]):
                faults.append(i)
                lats.append(lat)
                lons.append(lon)

    outdict = {"FAULT": faults, "LAT": lats, "LON": lons}
    import numpy as np
    np.savez("north_island_offshore_faults.npz", **outdict)


def onshore_faults():
    import numpy as np

    # 550_641 mesh
    map_corners = [-42.5007, -36.9488, 172.9998, 179.5077]

    with open("./new_zealand_onshore_faults.xy", 'r') as f:
        lines = f.readlines()
    faults, lats, lons = [], [], []
    i = 0
    print(len(lines))
    for line in lines[2:]:
        if line[0] == '>':
            i += 1
            print(i)
            continue
            # fault = line.split('"')[1]
            # faults.append(fault)
        else:
            line = line.split('-')
            lon = float(line[0].strip())
            lat = -1 * float(line[1])
            if (map_corners[0] < lat < map_corners[1]) and \
                    (map_corners[2] < lon < map_corners[3]):
                faults.append(i)
                lats.append(lat)
                lons.append(lon)

    outdict = {"FAULT": faults, "LAT": lats, "LON": lons}
    np.savez("north_island_onshore_faults.npz", **outdict)


