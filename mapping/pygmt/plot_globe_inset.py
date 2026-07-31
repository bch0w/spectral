"""
Inset map showing a globe with a given lat/lon location to highlight a region
of interest from a global perspective
"""
import numpy as np
import pygmt

GULKANA_LAT = 63.27
GULKANA_LON = -145.45
FBX_LAT = 64.859
FBX_LON = -147.849
ARCTIC_CIRCLE_LAT = 66.5607  # degrees North

# Create a new instance or object of the pygmt.Figure() class
fig = pygmt.Figure()
pygmt.config(MAP_FRAME_PEN="1.5p,black")
fig.coast(
    # Orthographic projection (G) with projection center at 0° East and
    # 15° North and a width of 12 centimeters
    projection="G-145/63/12c",
    region="g",  # global
    land="seashell",  # Color land masses in "gray"
    water="lightblue",  # Color water masses in "lightblue"
    # Add coastlines with a 0.25-point thick pen in "gray50"
    shorelines="1/.8p,black",
    borders="1/0.25p,black",
)

# Draw the black border/frame around the globe
fig.basemap(frame=True)

# Build the Arctic Circle as a line of constant latitude around the globe
circle_lons = np.arange(-180, 180.1, 1)
circle_lats = np.full_like(circle_lons, ARCTIC_CIRCLE_LAT)

fig.plot(
    x=circle_lons,
    y=circle_lats,
    pen="1p,purple,--",  # dashed black line, 1 point thick
)

fig.plot(x=GULKANA_LON, y=GULKANA_LAT, style="d0.2c", fill="red", pen="black")
# fig.text(
#     x=GULKANA_LON,
#     y=GULKANA_LAT,
#     text="Gulkana Glacier",
#     font="10p,Helvetica-Bold,black",
#     justify="LM",   # Left-Middle: text starts to the right of the point
#     offset="0.25c/0c",  # nudge the text away from the marker
# )
fig.plot(x=FBX_LON, y=FBX_LAT, style="c0.2c", fill="green", pen="black")
# fig.text(
#     x=FBX_LON,
#     y=FBX_LAT,
#     text="Fairbanks",
#     font="10p,Helvetica-Bold,black",
#     justify="LM",   # Left-Middle: text starts to the right of the point
#     offset="0.25c/0c",  # nudge the text away from the marker
# )
fig.savefig("globalinset.png",transparent=True)
fig.show()

