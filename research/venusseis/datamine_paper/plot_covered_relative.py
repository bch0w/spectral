#!/usr/bin/env python
"""
Plot noise at KAED and RND relative to each site's Covered (C) station.

These two sites are left out of the permanent-station comparison
(plot_permanent_relative.py) because their permanent station is not a
usable baseline: HV.KAED.EHZ appears broken and has no N/E channels, and
AK.RND shows a high-frequency offset shared by all XX stations. Here the
Uncovered and Windscreen stations are compared with the Covered station
hour by hour instead, using the same method and layout; the permanent
station is left out.

All options of plot_permanent_relative.py are available; only the defaults
differ (``--reference C --omit permanent --sites KAED,RND``). Pass
``--omit ''`` to draw the permanent station as a curve again.

Example::

    python plot_covered_relative.py
    python plot_covered_relative.py --sites RND --components Z,H
"""
from plot_permanent_relative import main


if __name__ == "__main__":
    main(description=__doc__, reference="C", omit="permanent",
         sites="KAED,RND", exclude_sites="",
         outfile="FIGURES/covered_relative.png",
         panel_height=2.4 )
