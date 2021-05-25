"""
Use PyGMT to plot slip rate deficit
"""
from pygmt import Figure, makecpt

# Set up colormap for plate coupling
cpt = "/Users/Chow/Documents/academic/vuw/data/carto/plate_coupling/tmp.cpt"
makecpt(cmap=cpt, series=[0, 1, 0.01])

# Set up standard map look
fig = Figure()
fig.coast(region=[173, 178.5, -42.5, -37.],
          projection="U60S/12c",
          land="whitesmoke",
          water="white",
          shorelines=["1/1.5p,black", "2/1p,black"],
          frame="afg"
          )

# Plot the plate coupling / slip rate deficit
srd = "/Users/Chow/Documents/academic/vuw/data/carto/plate_coupling/hik.gmt"
fig.plot(data=srd, 
         cmap=True, 
         close=True, 
         color="+z",
         )
fig.colorbar(frame="n")

# Plot the SSEs
sse = "/Users/Chow/Documents/academic/vuw/data/carto/sse_contours/filter_SSE2002_2014.grd"
fig.grdcontour(sse,
               pen=["1"],
               # interval=100,
               limit=[99, 550],
               )

# Plot the locations of the seamounts
fig.plot(x=176.609, y=-40.300, style="c0.3c", color="yellow", pen="2.,black")
fig.plot(x=177.903, y=-39.134, style="c0.3c", color="yellow", pen="2.,black")
fig.plot(x=177.108, y=-40.515, style="c0.3c", color="cyan", pen="2.,black")

fig.savefig("geodetic_setting.png")
fig.show(method="external")
