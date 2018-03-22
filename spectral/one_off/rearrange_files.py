"""
renaming and reorganizing mseed files from RDF deployment from output of
Apollo to the SEED data structure convention
"""
import os
from os.path import join as join
import glob

cwd = "/seis/prj/fwi/bchow/RDF/"
os.chdir(cwd)

SEED_name_template = "{net}.{sta}.{loc}.{cha}.{year}.{jday}"
files_to_rename = glob.glob(join(cwd,"*.ms"))
for old_filename in files_to_rename:
    sta,year,jday,loc_cha,net,fmt = os.path.basename(old_filename).split('.')
    prepath = os.path.dirname(old_filename)
    loc,cha = loc_cha.split('-')
    new_SEED_name = SEED_name_template.format(net=net,
                                            sta=sta,
                                            loc=loc,
                                            cha=cha,
                                            year=year,
                                            jday=jday)
    new_SEED_path = join(prepath,year,net,sta,cha,new_SEED_name)
    os.rename(src=old_filename,dst=new_SEED_path)
