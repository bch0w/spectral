"""the figure outputs of adjointBuilder don't play nice together, so they
are outputted as two separate files. this script will combine them into
a single figure
"""
import os
import sys
import glob
from PIL import Image

# filter bounds for naming convention
t0,t1 = 10,30

filepath = '/Users/chowbr/Documents/subduction/plots/adjtomo/2018p130600'
wavpath = os.path.join(filepath,'*_{t0}_{t1}*.png'.format(t0=t0,t1=t1))

for wav in glob.glob(wavpath):
    stationname = os.path.basename(wav).split('_')[0]
    mapfile = os.path.join(filepath,'{}_map.png'.format(stationname))
    if os.path.exists(mapfile):
        images = map(Image.open,[wav,mapfile])
        widths,heights = zip(*(i.size for i in images))
        
        total_width = sum(widths)
        max_height = max(heights)
        
        new_image = Image.new('RGB',(total_width,max_height),
                                                        color=(255,255,255,0))
        
        x_offset = 0
        images = map(Image.open,[wav,mapfile])
        for im in images:
            new_image.paste(im,(x_offset,0))
            x_offset += im.size[0]
        new_image.save(os.path.join(filepath,'{}_composite.png'.format(
                                                                stationname)))
        
    