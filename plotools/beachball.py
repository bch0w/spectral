"""
ObsPy's Beachball plotter is wrong. Use PyGMT to quickly plot a beachball given
SDR or MT components. For rapid visualization and input into figures

.. rubric::

    python beachball.py <strike> <dip> <rake> 

    OR

    python beachball.py <Mrr> <Mtt> <Mpp> <Mrt> <Mrp> <Mtp>
"""
import sys
import pygmt

vals = [float(_) for _ in sys.argv[1:]]

if len(vals) == 3:
    print("strike dip rake")
    convention = "aki"
    magnitude = 8
    spec = {"strike": vals[0], "dip": vals[1], "rake": vals[2], 
            "magnitude": magnitude}
elif len(vals) == 6:
    print("mrr, mtt, mff, mrt, mrf, mtf")
    convention = "mt"
    exponent = 23
    spec = {"mrr":vals[0], "mtt":vals[1], "mff":vals[2], 
            "mrt":vals[3], "mrf":vals[4], "mtf":vals[5],
            "exponent":exponent
            }

fig = pygmt.Figure()
fig.basemap(region=[-1, 1, -1, 1], projection="X10c/4c", frame=["a"])
fig.meca(spec=spec, convention=convention, scale="1c", longitude=0, latitude=0, 
         depth=0, compression_fill="red", extension_fill="cornsilk",
         pen="1.5p,black,solid")
fig.show()
