# AKATOM 
General project notes for orthern Alaska and northeastern Alaska

## REGION (NALASKA)
LAT: 61.8, 71.4 (~10)
LON: -172, -166 (~8)
LON: -141, -135 (~6)

LAT_CENTER: 67.25 (old: 68.25)
LON_CENTER: -154.0
LAT_WIDTH (ETA): 10
LON_WIDTH (XI): 13

TOP: 957 km
BOT: 1334 km
LEFT: 836 km
RIGHT: 836 km

MESH BOTTOM: 600 km

## REGION (NEALASKA)
LAT: 66, 71 (5)
LON: -151, -140 (11)

LAT_CENTER: 68.5 
LON_CENTER: -145.5

TOP: 360 km (adjusted to 442)
BOT: 442 km
LEFT: 551 km
RIGHT: 551 km 

MESH BOTTOM: 400 km

## CHINOOK
- T1SMALL 28 CORES/NODE
- NPROC_XI * NPROC_ETA = 8 * 7 = 56 CORES (2 nodes)
- NEX_XI = 128 (1000 / 128 = 7.8 km / element)
- NEX_ETA = 112 (836 / 112 = 7.4 km / element)

## RELEVANT
https://earthquake.alaska.edu/finding-faults-northeast-alaska
https://earthquake.alaska.edu/m65-kaktovik-earthquake-largest-ever-north-slope

## ANAT
- Point forces at the surface (Wang et al. 2018, 2019)
- 3-station EGF are valid for 8-50s


