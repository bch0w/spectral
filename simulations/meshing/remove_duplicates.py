"""
My tomo file editing created duplicated entries in the new outputted tomo files
This may have been due to the iterpolation (?), or due to how I was writing
the files out
This quick script will remove those duplicates, check that the new tomo files
are okay, and make depth slice plots for visual inspection
"""
import os
import numpy as np

# qty = 'crust'
for qty in ['mantle', 'crust', 'shallow']:
    tag = f'nz_utm60_2km_1700x_5271y_eberhart15_{qty}'
    fid = f'{tag}.xyz.npy'
    print(qty)

    # Read in the datafile that contains duplicates
    with_dupes = np.load(fid)
    nx_wd = len(np.unique(with_dupes[:, 0]))
    ny_wd = len(np.unique(with_dupes[:, 1]))
    nz_wd = len(np.unique(with_dupes[:, 2]))

    assert(len(with_dupes) != nx_wd * ny_wd * nz_wd)

    # Create a new array without duplicates
    _, nodupes_ind = np.unique(with_dupes[:, :3], axis=0, return_index=True)
    no_dupes = with_dupes[nodupes_ind]
    nx_nd = len(np.unique(no_dupes[:, 0]))
    ny_nd = len(np.unique(no_dupes[:, 1]))
    nz_nd = len(np.unique(no_dupes[:, 2]))

    # Check the duplicates and see if the values are the same
    with_dupes_ind = np.arange(len(with_dupes)) 
    duplicates = np.setdiff1d(with_dupes_ind, nodupes_ind)

    # Loop through the indices that are duplicated
    print(len(duplicates))
    # for dupe in duplicates:
    #     dupe_line = with_dupes[dupe]
    #     duplicate_entries = with_dupes[np.where((with_dupes[:, 0] == dupe_line[0]) & 
    #                                             (with_dupes[:, 1] == dupe_line[1]) &
    #                                             (with_dupes[:, 2] == dupe_line[2]))
    #                                             ]
    #     # Check that the duplicated values are the same, only check vp (ind 3)
    #     assert(len(duplicate_entries) == 2)
    #     assert(duplicate_entries[0][3] == duplicate_entries[1][3])

    assert(len(no_dupes) == nx_nd * ny_nd * nz_nd)


    # rename the file with duplicates
    fid_rename = f'{tag}_with_dupes.xyz.npy'
    os.rename(fid, fid_rename)

    # make sure that worked, and then write the new file
    if not os.path.exists(fid):
        np.save(fid, no_dupes)

