"""
Plot a moment tensor focal mechanism (beachball) with PyGMT.

Source can be given either as a CMTSOLUTION file (as output by SPECFEM /
Global CMT) or as the six moment tensor components directly on the command
line.

Usage
-----
# From a CMTSOLUTION file
python plot_beachball_pygmt.py --cmtsolution ./CMTSOLUTION_NK6 \
    --save beachball.png

# From explicit moment tensor components (Harvard/GCMT r,t,p convention, in
# dyne-cm - the convention used by CMTSOLUTION/GCMT and GMT's `meca`)
python plot_beachball_pygmt.py --mrr 7.14e23 --mtt 5.17e23 --mpp 4.27e23 \
    --mrt 1.30e22 --mrp -2.59e23 --mtp 4.51e22 --save beachball.png
"""
import sys
from argparse import ArgumentParser

import numpy as np
import pygmt
from obspy import read_events


def mantissa_exponent(mrr, mtt, mpp, mrt, mrp, mtp):
    """
    GMT's `meca` expects moment tensor components as a mantissa (~1-10)
    plus a single shared power-of-ten exponent, e.g., 7.14e+23 -> mantissa
    7.14, exponent 23. Passing raw, un-normalized values (as returned by
    e.g., ObsPy or listed directly in a CMTSOLUTION file) breaks GMT's
    internal size scaling, so all tensor components need to be normalized
    by the same exponent before being handed to `meca`

    :rtype: tuple
    :return: (mrr, mtt, mpp, mrt, mrp, mtp, exponent) with the 6 tensor
        components rescaled to mantissa form
    """
    vals = [mrr, mtt, mpp, mrt, mrp, mtp]
    max_abs = max(abs(val) for val in vals)
    if max_abs == 0:
        return (*vals, 0)

    exponent = int(np.floor(np.log10(max_abs)))
    scaled = [val / 10 ** exponent for val in vals]

    return (*scaled, exponent)


def read_cmtsolution(path_cmtsolution):
    """
    Read a CMTSOLUTION file and return the 6 moment tensor components
    (r,t,p convention, in dyne-cm) plus the event depth in km and a label

    :type path_cmtsolution: str
    :param path_cmtsolution: path to a CMTSOLUTION (or QuakeML-readable)
        source file
    :rtype: tuple
    :return: (mrr, mtt, mpp, mrt, mrp, mtp, depth_km, label)
    """
    event = read_events(path_cmtsolution)[0]
    mt = event.focal_mechanisms[0].moment_tensor.tensor
    origin = event.preferred_origin() or event.origins[0]

    label = event.resource_id.id.split("/")[-1]

    # ObsPy stores the tensor in SI units (N-m); GMT's `meca` expects the
    # GCMT/CMTSOLUTION convention (dyne-cm), so convert back: 1 N-m = 1E7 dyne-cm
    nm_to_dynecm = 1E7

    return (mt["m_rr"] * nm_to_dynecm, mt["m_tt"] * nm_to_dynecm,
            mt["m_pp"] * nm_to_dynecm, mt["m_rt"] * nm_to_dynecm,
            mt["m_rp"] * nm_to_dynecm, mt["m_tp"] * nm_to_dynecm,
            origin.depth / 1E3, label)


def plot_beachball(mrr, mtt, mpp, mrt, mrp, mtp, depth_km=0,
                    label=None, compression_fill="red",
                    extension_fill="white", scale="4c", region=None,
                    projection="X6c", save="beachball.png", dpi=300,
                    show=False):
    """
    Plot a single moment tensor focal mechanism as a PyGMT beachball

    :type mrr/mtt/mpp/mrt/mrp/mtp: float
    :param mrr/mtt/mpp/mrt/mrp/mtp: moment tensor components in the
        Harvard/GCMT (r, t, p) convention, in dyne-cm. GMT's `meca` derives
        the symbol's physical size from the scalar moment assuming dyne-cm
        units, so values in other units (e.g., ObsPy's N-m) must be
        converted to dyne-cm first, see `read_cmtsolution`. Values are
        auto-scaled to the mantissa+exponent form GMT expects internally,
        see `mantissa_exponent`
    :type depth_km: float
    :param depth_km: source depth in km, used only for the `meca` spec
    :type label: str
    :param label: optional text label plotted above the beachball
    :type compression_fill: str
    :param compression_fill: fill color for compressional (P) quadrants
    :type extension_fill: str
    :param extension_fill: fill color for extensional (T) quadrants
    :type scale: str
    :param scale: GMT-style size of the beachball, e.g. '4c'
    :type region: list
    :param region: GMT region [xmin, xmax, ymin, ymax], defaults to a small
        box centered on the beachball
    :type projection: str
    :param projection: GMT projection string
    :type save: str
    :param save: output filename for the figure
    :type dpi: int
    :param dpi: resolution of the saved figure
    :type show: bool
    :param show: open the figure in an external viewer after plotting
    """
    region = region or [-1, 1, -1, 1]
    longitude, latitude = 0, 0

    mrr, mtt, mpp, mrt, mrp, mtp, exponent = \
        mantissa_exponent(mrr, mtt, mpp, mrt, mrp, mtp)
    tensor_dict = {"mrr": mrr, "mtt": mtt, "mff": mpp, "mrt": mrt,
                    "mrf": mrp, "mtf": mtp, "exponent": exponent}

    fig = pygmt.Figure()
    fig.basemap(region=region, projection=projection, frame=0)
    fig.meca(spec=tensor_dict, scale=scale, longitude=longitude,
              latitude=latitude, depth=depth_km,
              compression_fill=compression_fill,
              extension_fill=extension_fill, pen="0.5p,black,solid")

    if label:
        fig.text(x=longitude, y=region[-1] * 0.9, text=label, justify="CB",
                  font="12p")

    fig.savefig(save, dpi=dpi)
    if show:
        fig.show()

    return fig


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cmtsolution", default=None,
                         help="path to a CMTSOLUTION file to read the "
                              "moment tensor from")
    parser.add_argument("--mrr", type=float, default=None)
    parser.add_argument("--mtt", type=float, default=None)
    parser.add_argument("--mpp", type=float, default=None)
    parser.add_argument("--mrt", type=float, default=None)
    parser.add_argument("--mrp", type=float, default=None)
    parser.add_argument("--mtp", type=float, default=None)
    parser.add_argument("--depth_km", type=float, default=0,
                         help="source depth in km (ignored if "
                              "--cmtsolution is given)")
    parser.add_argument("--label", default=None,
                         help="text label to plot above the beachball "
                              "(default: CMTSOLUTION event name, if given)")
    parser.add_argument("--compression_fill", default="red",
                         help="fill color for compressional quadrants "
                              "(default: 'red')")
    parser.add_argument("--extension_fill", default="white",
                         help="fill color for extensional quadrants "
                              "(default: 'white')")
    parser.add_argument("--scale", default="4c",
                         help="GMT-style beachball size (default: '4c')")
    parser.add_argument("--save", default="beachball.png",
                         help="output figure filename (default: "
                              "'beachball.png')")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--show", action="store_true", default=False,
                         help="open the figure after plotting")
    args = parser.parse_args()

    if args.cmtsolution:
        mrr, mtt, mpp, mrt, mrp, mtp, depth_km, label = \
            read_cmtsolution(args.cmtsolution)
        label = args.label or label
    else:
        required = [args.mrr, args.mtt, args.mpp, args.mrt, args.mrp,
                    args.mtp]
        if any(val is None for val in required):
            sys.exit("must provide either --cmtsolution or all six of "
                      "--mrr --mtt --mpp --mrt --mrp --mtp")
        mrr, mtt, mpp, mrt, mrp, mtp = required
        depth_km = args.depth_km
        label = args.label

    plot_beachball(mrr=mrr, mtt=mtt, mpp=mpp, mrt=mrt, mrp=mrp, mtp=mtp,
                    depth_km=depth_km, label=label,
                    compression_fill=args.compression_fill,
                    extension_fill=args.extension_fill, scale=args.scale,
                    save=args.save, dpi=args.dpi, show=args.show)
