"""
Utility script to convert the Zenodo collection of Liu et al. (2022) 
Empirical Green's Function data (three-station ambient noise cross-correlations)
from the default/native format into a directory structure and formatting that
is suitable for use in SeisFlows. 

This was originally done through a series of functions and scripts over time but
I knew I would forget the overall workflow so it's better to have everything in
one place.

.. warning::

    I never ran this as a script, just copy-pasted into IPython

This script takes into account a few things:

1) If a STATIONS file (SPECFEM format) is provided, it will be used to exclude
   stations NOT included in that file (e.g., to set a simulation domain)
2) SAC stats will be set correctly 
3) Origintime (originally 2019-01-01T00:00:00) can be set manually
4) Source-Receiver plots will be made for each sub-directory

Noteworthy points:
- We are using 'hyp' (hyperbolic) data ONLY because it has higher SNR than 'ell'
"""
import os
import shutil
import matplotlib.pyplot as plt
from glob import glob
from obspy import read, UTCDateTime
from pysep.utils.io import read_stations, write_stations_file


# SET PARAMETERS HERE
INPUT_DIR = "/home/bchow/Work/data/egfs/SAC_I3_stack_4_Zendo/"
OUTPUT_DIR = "/home/bchow/Work/data/egfs/LIU22_EGF"  # !!!
SOURCES_FILE = "SOURCES_FINAL_NALASKA"
STATIONS_FILE = "STATIONS_FINAL_NALASKA"
STATIONS_ALL = "data/ORIGINAL_STATIONS_NALASKA"


def restructure_directory():
    """
    Creates symlinks of each of the default SAC files that is renamed and
    re-organized to fit the directory structure. Does NOT touch the original
    data at all. The resulting directories will be formatted:

    <ELL_OR_HYP>/<SOURCE_STATION>/<ZZ_OR_TT>/<RECEIVER_STATION>

    <RECEIVER_STATION> will be formatted NN.SSS.CXC.SAC
    """
    for stack in glob(os.path.join(INPUT_DIR, "*")):
        # e.g., Lov_I3_hyp_stack
        phase, _, type_, _ = os.path.basename(stack).split("_")

        kernel = {"Lov": "TT", "Ray": "ZZ"}[phase]                
        comp = kernel[-1]                                        
        for src in glob(os.path.join(stack, "*")):
            for fid in glob(os.path.join(src, "*")):                      
                # Separate out the station names from rest of file name       
                _fid = os.path.splitext(os.path.basename(fid))[0]           
                _, net_src, sta_src, net_rcv, sta_rcv = _fid.split("_")       

                # Determine the correct save location                              
                path_out = os.path.join(OUTPUT_DIR, type_.lower(),
                                        f"{net_src}_{sta_src}", kernel)
                if not os.path.exists(path_out):                                    
                    os.makedirs(path_out)                                          

                # Label the station correctly, L because it's 1Hz data, X 
                # because it's a time series derived from observational data
                filename_out = f"{net_rcv}.{sta_rcv}.LX{comp}.SAC"                
                fid_out = os.path.join(path_out, filename_out)                     
                if not os.path.exists(fid_out):                                  
                    os.symlink(fid, fid_out)


def fix_restructured_directory():
    """
    Goes back through restructured directory to fix SAC headers and rename files
    """
    for type_ in glob(os.path.join(INPUT_DIR, "*")):
        for src_sta in glob(os.path.join(type_, "*")):
            for kernel in glob(os.path.join(src_sta, "*")):
                for fid in glob(os.path.join(kernel, "*")):
                    # Separate out the station names from rest of file name       
                    _fid = os.path.splitext(os.path.basename(fid))[0]           
                    net, sta, cha = _fid.split(".")

                    fid_out = fid.replace(".SAC", ".sac")
                    os.rename(fid, fid_out)  # .SAC -> .sac
                    # Correct the SAC header
                    try:
                        st = read(fid_out)
                        if len(st) != 1:
                            print(f"ERROR: {_fid} too many traces")
                            continue
                        if st[0].stats.station == sta:
                            continue

                        for tr in st:
                            tr.stats.network = net
                            tr.stats.station = sta
                            tr.stats.channel = cha
                            tr.stats.sac.evdp = 0.
                        st.write(fid_out, format="SAC")
                        print(f"SUCCESS: {fid}")
                    except Exception as e:
                        print(f"ERROR: {fid} {e}")
                        continue

                    a=1/0


def restructure_directory_mulligan():
    """
    If you have to do it again, edit SAC headers at the same time

    Creates symlinks of each of the default SAC files that is renamed and
    re-organized to fit the directory structure. Does NOT touch the original
    data at all. The resulting directories will be formatted:

    <ELL_OR_HYP>/<SOURCE_STATION>/<ZZ_OR_TT>/<RECEIVER_STATION>

    <RECEIVER_STATION> will be formatted NN.SSS.CXC.SAC
    """
    for stack in glob(os.path.join(INPUT_DIR, "*")):
        # e.g., Lov_I3_hyp_stack
        phase, _, type_, _ = os.path.basename(stack).split("_")

        kernel = {"Lov": "TT", "Ray": "ZZ"}[phase]                
        comp = kernel[-1]                                        
        for src in glob(os.path.join(stack, "*")):
            for fid in glob(os.path.join(src, "*")):                      
                # Separate out the station names from rest of file name       
                _fid = os.path.splitext(os.path.basename(fid))[0]           
                _, net_src, sta_src, net_rcv, sta_rcv = _fid.split("_")       

                # Determine the correct save location                              
                path_out = os.path.join(OUTPUT_DIR, type_.lower(),
                                        f"{net_src}_{sta_src}", kernel)
                if not os.path.exists(path_out):                                    
                    os.makedirs(path_out)                                          

                # Label the station correctly, L because it's 1Hz data, X 
                # because it's a time series derived from observational data
                filename_out = f"{net_rcv}.{sta_rcv}.LX{comp}.SAC"                
                fid_out = os.path.join(path_out, filename_out)                     
                if not os.path.exists(fid_out):                                  
                    os.symlink(fid, fid_out)
                    # Correct the SAC header
                    try:
                        net, sta, cha, sac = basename.split(".")
                        st = read(fid_out)
                        if len(st) != 1:
                            print(f"ERROR: {basename} too many traces")
                            return
                        for tr in st:
                            # Correct the SAC header
                            tr.stats.network = net_rcv
                            tr.stats.station = sta_rcv
                            tr.stats.channel = f"LX{comp}"
                            tr.stats.sac.evdp = 0.
                        st.write(fid, format="SAC")
                        print(f"SUCCESS: {basename}")
                    except Exception as e:
                        print(f"ERROR: {basename} {e}")
                        return


def count_source_receiver_hits(src_threshold=86, rcv_threshold=5):
    """
    Count the number of data for each source and each receiver station as a way
    of determining which ones to keep and which to kick. Sets a threshold for
    minimum number of measurements a source and receiver station can have.
    """
    inv_srcs = read_stations(SOURCES_FILE)
    inv_rcvs = read_stations(STATIONS_FILE)

    src_count, rcv_count = {}, {}
    for src_dir in glob(os.path.join(OUTPUT_DIR, "*")):
        src_net, src_sta = os.path.basename(src_dir).split("_")
        if not inv_srcs.select(network=src_net, station=src_sta):
            continue

        src_code = f"{src_net}_{src_sta}"
        if src_code not in src_count:
            src_count[src_code] = 0

        for kernel in glob(os.path.join(src_dir, "*")):
            for rcv_fid in glob(os.path.join(kernel, "*")):
                rcv_net, rcv_sta, *_ = os.path.basename(rcv_fid).split(".")
                if not inv_rcvs.select(network=rcv_net, station=rcv_sta):
                    continue

                rcv_code = f"{rcv_net}_{rcv_sta}"
                if rcv_code not in rcv_count:
                    rcv_count[rcv_code] = 1
                else:
                    rcv_count[rcv_code] += 1
                src_count[src_code] += 1

    # Sort dictionary by values
    src_count = {k: v for k, v in sorted(src_count.items(),
                                         key=lambda item: item[1])}

    # Sort dictionary by values
    rcv_count = {k: v for k, v in sorted(rcv_count.items(),
                                         key=lambda item: item[1])}

    # Apply thresholding
    if src_threshold:
        src_count_out = {}
        for key, val in src_count.items():
            if val >= src_threshold:
                src_count_out[key] = val
    else:
        src_count_out = src_count

    if rcv_threshold:
        rcv_count_out = {}
        for key, val in rcv_count.items():
            if val >= rcv_threshold:
                rcv_count_out[key] = val
    else:
        rcv_count_out = rcv_count

    return src_count_out, rcv_count_out


def write_source_station_files_with_thresholded_counts():
    """Write out stations files"""
    src_count, rcv_count = count_source_receiver_hits()

    inv = read_stations(STATIONS_ALL)
    src_inv = read_stations(SOURCES_FILE)
    rcv_inv = read_stations(STATIONS_FILE)
    print(f"SRC: {len(src_inv.get_contents()['stations'])}")
    print(f"RCV: {len(rcv_inv.get_contents()['stations'])}")

    for net in inv:
        for sta in net:
            if f"{net.code}_{sta.code}" in src_count:
                pass
            else:
                src_inv = src_inv.remove(network=net.code, station=sta.code)
            if f"{net.code}_{sta.code}" in rcv_count:
                pass
            else:
                rcv_inv = rcv_inv.remove(network=net.code, station=sta.code)

    print(f"SRC: {len(src_inv.get_contents()['stations'])}")
    print(f"RCV: {len(rcv_inv.get_contents()['stations'])}")

    write_stations_file(src_inv, "SOURCES_FINAL_NALASKA")
    write_stations_file(rcv_inv, "STATIONS_FINAL_NALASKA")


def plot_source_receiver_hits(src_threshold=86, rcv_threshold=5):
    """
    Plot src_count and rcv_count to get a map of coverage at each station loc.
    """
    src_count, rcv_count = count_source_receiver_hits(src_threshold,
                                                      rcv_threshold)

    inv = read_stations(STATIONS_ALL)

    # Plot the number of hits and threshold values to provide a 'streamlined'
    # number of source and receiver stations
    for i, counts in enumerate([src_count, rcv_count]):
        lons, lats, cnts, codes = [], [], [], []
        for code, count in counts.items():
            net, sta = code.split("_")
            inv_ = inv.select(network=net, station=sta)

            lons.append(inv_[0][0].longitude)
            lats.append(inv_[0][0].latitude)
            cnts.append(count)
            codes.append(code)

        f = plt.figure(figsize=(20, 10), dpi=200)
        plt.scatter(lons, lats, c=cnts, marker="o",
                    zorder=6, s=40, ec="k", lw=1)
        for lon, lat, code in zip(lons, lats, codes):
            plt.text(s=f"{code}: {counts[code]}", x=lon, y=lat, size=7,
                     zorder=8)
        plt.colorbar(label="counts")

        # Plot all stations in list just for reference of what has been kicked
        for net in inv:
            for sta in net:
                if f"{net.code}_{sta.code}" in counts:
                    continue
                plt.scatter(sta.longitude, sta.latitude, c="k", marker="x",
                            alpha=0.5, s=75, zorder=5)
                # plt.text(s=f"{net.code}.{sta.code}", x=sta.longitude,
                #          y=sta.latitude, size=7, zorder=8)

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        key = ["source", "receiver"][i]
        plt.title(f"{len(lons)} {key} stations\n"
                  f"N_measurements=[{min(cnts)}, {max(cnts)}];\n"
                  )
        plt.show()
        plt.close()
        # plt.savefig(f"coverage_figures/src_count.png")


def remove_doubled_sources():
    """
    AK and TA stations share locations but have different datasets. Need to
    remove one of them to avoid having two simulations with almost exactly the
    same source characeristics
    """
    src_threshold = 125
    rcv_threshold = 20
    src_count, rcv_count = count_source_receiver_hits(src_threshold,
                                                      rcv_threshold)

    # Remove source doubles by choosing the higher number of measurements
    doubles, checker = {}, {}
    for code, count in src_count.items():
        net, sta = code.split("_")
        if sta not in checker:
            checker[sta] = (code, count)
        else:
            code_check, count_check = checker[sta]
            if count > count_check:
                doubles[sta] = f"{code} {count} > {code_check} {count_check} " \
                               f"[{count - count_check}]"
            else:
                doubles[sta] = f"{code} {count} < {code_check} {count_check} " \
                               f"[{count_check - count}]"

    _keys = list(doubles.keys())
    _keys.sort()
    doubles = {i: doubles[i] for i in _keys}

    with open("check_doubles.txt", "w") as f:
        for key, val in doubles.items():
            f.write(f"    {key}: {val}\n")

            
def plot_measurement_coverage():
    """
    Creates very simple scatter plots showing each source station and the
    corresponding receiver stations which have data. Both to show coverage
    and data availability
    """
    inv_srcs = read_stations(SOURCES_FILE)
    inv_rcvs = read_stations(STATIONS_FILE)

    for src_dir in glob(os.path.join(OUTPUT_DIR, "*")):
        src_net, src_sta = os.path.basename(src_dir).split("_")
        src_inv = inv_srcs.select(network=src_net, station=src_sta)
        if not src_inv:
            print(f"SRC ERROR: {src_net}.{src_sta}")
            continue
        src_lat = src_inv[0][0].latitude
        src_lon = src_inv[0][0].longitude

        z_rcv_lats, z_rcv_lons = [], []
        t_rcv_lats, t_rcv_lons = [], []
        for kernel in glob(os.path.join(src_dir, "*")):
            kernel_name = os.path.basename(kernel)
            for rcv_fid in glob(os.path.join(kernel, "*")):
                rcv_net, rcv_sta, *_ = os.path.basename(rcv_fid).split(".")
                rcv_inv = inv_rcvs.select(network=rcv_net, station=rcv_sta)
                if not rcv_inv:
                    print(f"RCV ERROR: {rcv_net}.{rcv_sta}")
                    continue
                if kernel_name == "ZZ":
                    z_rcv_lats.append(rcv_inv[0][0].latitude)
                    z_rcv_lons.append(rcv_inv[0][0].longitude)
                elif kernel_name == "TT":
                    t_rcv_lats.append(rcv_inv[0][0].latitude)
                    t_rcv_lons.append(rcv_inv[0][0].longitude)

        n = len(z_rcv_lats) + len(t_rcv_lats)
        plt.scatter(src_lon, src_lat, c="r", marker="o", zorder=6)
        plt.scatter(z_rcv_lons, z_rcv_lats, c="b", marker="v", zorder=5)
        plt.scatter(t_rcv_lons, t_rcv_lats, c="g", marker="^", zorder=5)
        for rcv_lon, rcv_lat in zip(z_rcv_lons, z_rcv_lats):
            plt.plot([src_lon, rcv_lon], [src_lat, rcv_lat], c="k",
                     ls="-", alpha=0.25, zorder=4)
        for rcv_lon, rcv_lat in zip(t_rcv_lons, t_rcv_lats):
            plt.plot([src_lon, rcv_lon], [src_lat, rcv_lat], c="k",
                     ls="-", alpha=0.25, zorder=4)

        # Plot all stations for refernce
        for net in inv_rcvs:
            for sta in net:
                plt.scatter(sta.longitude, sta.latitude, c="k", marker="s",
                            zorder=3, alpha=0.25)

        plt.xlim([-168, -140])
        plt.ylim([64.5, 72.])
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.title(f"{src_net}.{src_sta}; {kernel_name}; N={n}")
        plt.savefig(f"coverage_figures/{n}_{src_net}_{src_sta}_{kernel_name}.png")
        plt.close("all")


def _combine_coverage_figures():
    """Combine coverage figures for repeated TA and AK stations to decide
    which ones to keep based visually on coverage"""
    from PIL import Image

    png_files = glob("coverage_figures/*_*_*_*.png")
    stations = [_.split("_")[3] for _ in png_files]
    for sta in stations:
        ak_img = glob(f"coverage_figures/*_AK_{sta}_*.png")
        ta_img = glob(f"coverage_figures/*_TA_{sta}_*.png")
        if ak_img and ta_img:
            images = [ak_img[0], ta_img[0]]
        else:
            continue

        images = [Image.open(_).convert("RGBA") for _ in images]

        widths, heights = zip(*(i.size for i in images))
        total_width = sum(widths)
        max_height = max(heights)

        # Create the new image that will be returned
        im_out = Image.new(mode="RGBA", size=(total_width, max_height))
        x_offset = 0
        for im in images:
            im_out.paste(im=im, box=(x_offset, 0))
            x_offset += im.size[0]

        im_out.save(f"coverage_figures/doubles/{sta}.png")


def plot_receiver_coverage():
    """
    Creates very simple scatter plots showing each receiver station and the
    corresponding source stations which have data.
    """
    inv_srcs = read_stations(SOURCES_FILE)
    inv_rcvs = read_stations(STATIONS_FILE)

    rcv_dict = {}
    for src_dir in glob(os.path.join(OUTPUT_DIR, "*")):
        src_net, src_sta = os.path.basename(src_dir).split("_")
        src_inv = inv_srcs.select(network=src_net, station=src_sta)
        if not src_inv:
            print(f"SRC ERROR: {src_net}.{src_sta}")
            continue
        src_lat = src_inv[0][0].latitude
        src_lon = src_inv[0][0].longitude

        for kernel in glob(os.path.join(src_dir, "*")):
            for rcv_fid in glob(os.path.join(kernel, "*")):
                rcv_net, rcv_sta, *_ = os.path.basename(rcv_fid).split(".")
                rcv_inv = inv_rcvs.select(network=rcv_net, station=rcv_sta)
                if not rcv_inv:
                    print(f"RCV ERROR: {rcv_net}.{rcv_sta}")
                    continue
                if f"{rcv_net}_{rcv_sta}" not in rcv_dict:
                    rcv_dict[f"{rcv_net}_{rcv_sta}"] = {"lats": [src_lat],
                                                        "lons": [src_lon]
                                                        }
                else:
                    rcv_dict[f"{rcv_net}_{rcv_sta}"]["lats"].append(src_lat)
                    rcv_dict[f"{rcv_net}_{rcv_sta}"]["lons"].append(src_lon)

    for rcv, vals in rcv_dict.items():
        rcv_net, rcv_sta = rcv.split("_")
        rcv_inv = inv_rcvs.select(network=rcv_net, station=rcv_sta)
        rcv_lat = rcv_inv[0][0].latitude
        rcv_lon = rcv_inv[0][0].longitude

        plt.scatter(rcv_lon, rcv_lat, c="r", marker="v", zorder=6)
        plt.scatter(vals["lons"], vals["lats"], c="g", marker="o", zorder=5)

        for src_lon, src_lat in zip(vals["lons"], vals["lats"]):
            plt.plot([src_lon, rcv_lon], [src_lat, rcv_lat], c="k",
                     ls="-", alpha=0.25, zorder=4)

        plt.xlim([-168, -140])
        plt.ylim([64.5, 72.])
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        n = len(vals['lons'])
        plt.title(f"{rcv_net}.{rcv_sta}; N={n}")
        plt.savefig(f"{n}_{rcv_net}_{rcv_sta}_rcv.png")
        plt.close("all")


def copy_trimmed_directory():
    """
    Creates copies of each of the default SAC files but only for the chosen
    SOURCE and STATION files and type (ELL or HYP). Files are copied because
    their data and metadata will be modified later and we don't want to affect
    the original dataset

    <ELL_OR_HYP>/<SOURCE_STATION>/<ZZ_OR_TT>/<RECEIVER_STATION>

    <RECEIVER_STATION> will be formatted NN.SSS.CXC.SAC
    """
    # OVERWRITE OUTPUT_DIR
    # OUTPUT_DIR = "/home/bchow/Work/data/egfs/NAKANAT_EGF"
    inv_srcs = read_stations(SOURCES_FILE)
    inv_rcvs = read_stations(STATIONS_FILE)

    for stack in glob(os.path.join(INPUT_DIR, "*")):
        # e.g., Lov_I3_hyp_stack
        phase, _, type_, _ = os.path.basename(stack).split("_")
        # We use hyp data because it has higher SNR (from paper)
        if type_.lower() == "ell":
            continue

        kernel = {"Lov": "TT", "Ray": "ZZ"}[phase]
        comp = kernel[-1]
        for src in glob(os.path.join(stack, "*")):
            # Source station needs to be in the simulation domain
            net_src, sta_src = os.path.basename(src).split("_")
            if not bool(inv_srcs.select(network=net_src, station=sta_src)):
                continue

            for fid in glob(os.path.join(src, "*")):
                # Separate out the station names from rest of file name
                _fid = os.path.splitext(os.path.basename(fid))[0]
                _, net_src, sta_src, net_rcv, sta_rcv = _fid.split("_")

                # Rcv station need to be in the simulation domain
                if not bool(inv_rcvs.select(network=net_rcv, station=sta_rcv)):
                    continue
                # Determine the correct save location
                path_out = os.path.join(OUTPUT_DIR, f"{net_src}_{sta_src}",
                                        kernel)
                if not os.path.exists(path_out):
                    os.makedirs(path_out)

                # Label the station correctly, L because it's 1Hz data, X
                # because it's a time series derived from observational data
                filename_out = f"{net_rcv}.{sta_rcv}.LX{comp}.SAC"
                fid_out = os.path.join(path_out, filename_out)
                if not os.path.exists(fid_out):
                    shutil.copy(fid, fid_out)
                    _overwrite_metadata_and_trim_length(fid_out)


def _overwrite_metadata_and_trim_length(
        fid, new_starttime="2000-01-01T00:00:00", new_length_s=60*20):
    """
    Set a new starttime and fix some wonky SAC headers that are present in the
    default data. Also trim the waveforms to a specified time duration
    """
    new_starttime = UTCDateTime(new_starttime)
    try:
        basename = os.path.basename(fid)
        net, sta, cha, sac = basename.split(".")
        st = read(fid)
        if len(st) != 1:
            print(f"ERROR: {basename} too many traces")
            return
        for tr in st:
            # Set the correct stats attribute
            tr.stats.network = net
            tr.stats.station = sta
            tr.stats.channel = cha
            # Set the correct SAC header
            tr.stats.starttime = new_starttime
            tr.stats.sac.nzyear = new_starttime.year
            tr.stats.sac.evdp = 0.
            # Trim data
            tr.trim(new_starttime, new_starttime + new_length_s)
        st.write(fid, format="SAC")
        print(f"SUCCESS: {basename}")
    except Exception as e:
        print(f"ERROR: {basename} {e}")
        return


def reoverwrite_data(new_starttime="2000-01-01T00:00:00",
                     new_length_s=60*20):
    """Second pass because the first run didn't catch all the data"""
    new_starttime = UTCDateTime(new_starttime)

    for src_sta in glob("*"):
        for kernel in glob(os.path.join(src_sta, "*")):
            for fid in glob(os.path.join(kernel, "*")):
                try:
                    basename = os.path.basename(fid)
                    net, sta, cha, sac = basename.split(".")
                    st = read(fid)
                    if len(st) != 1:
                        print(f"ERROR: {basename} too many traces")
                        return
                    if st[0].stats.station == sta:
                        continue

                    for tr in st:
                        # Set the correct stats attribute
                        tr.stats.network = net
                        tr.stats.station = sta
                        tr.stats.channel = cha
                        # Set the correct SAC header
                        tr.stats.starttime = new_starttime
                        tr.stats.sac.nzyear = new_starttime.year
                        tr.stats.sac.evdp = 0.
                        # Trim data
                        tr.trim(new_starttime, new_starttime + new_length_s)
                    st.write(fid, format="SAC")
                    print(f"SUCCESS: {basename}")
                except Exception as e:
                    print(f"ERROR: {basename} {e}")
                    continue
