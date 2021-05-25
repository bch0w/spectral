"""
Generate map of New Zealnd North using PyGMT
"""
import pygmt


# Initiate the map
fig = pygmt.Figure()
fig.coast(region=[173, 178.5, -42.5, -37.],
          projection="U60S/12c",
          land="whitesmoke",
          water="white",
          shorelines=["1/1.5p,black", "2/1p,black"],
          frame="afg"
          )
fig.show(method="external")


# Plot bathymetry from .grd file
grd_file = "/Users/Chow/Documents/academic/vuw/data/carto/topography/dem/bathymetryNZ.grd"
fig = pygmt.Figure()
fig.grdimage(grid=grd_file)

