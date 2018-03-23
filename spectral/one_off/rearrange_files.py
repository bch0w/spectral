"""
renaming and reorganizing mseed files from RDF deployment from output of
Apollo to the SEED data structure convention
"""
import os
from os.path import join as join
import glob

def check_files():
    """grab filenames to check if conversion went okay
    """
    filecheck = glob.glob("*.ms")
    jday_list,sta_list = [],[]
    for fid in filecheck:
        sta,year,jday,loc_cha,net,fmt = os.path.basename(fid).split('.')
        sta_list.append(sta)
        jday_list.append(jday)

    new_list = zip(sta_list,jday_list)
    new_list = sorted(new_list, key=lambda x: x[1])
    return new_list

def fixo():
    """some datafiles overlap between two deployments, combine into single
    """
    from os.path import join
    import os
    import glob
    for c in ['HHZ','HHE','HHN']:
        overlap_list = ["RD07.2017.264.*{}*".format(c),
                        "RD08.2017.264.*{}*".format(c),
                        "RD11.2017.264.*{}*".format(c),
                        "RD13.2017.264.*{}*".format(c),
                        "RD15.2017.263.*{}*".format(c),
                        "RD16.2017.263.*{}*".format(c)]
        list1 = '/seis/prj/fwi/yoshi/RDF_Array/Sep2017_Nov2017/DATA_ALL'
        list2 = '/seis/prj/fwi/yoshi/RDF_Array/July2017_Sep2017/DATA_ALL'
        SEED_name_template = "{net}.{sta}.{loc}.{cha}.{year}.{jday}"
        for entry in overlap_list:
            file1 = glob.glob(join(list1,entry))[0]
            file2 = glob.glob(join(list2,entry))[0]
            st = read(file1) + read(file2)
            prepath = '/seis/prj/fwi/bchow/RDF'
            sta,year,jday,loc_cha,net,fmt = os.path.basename(file1).split('.')
            loc,cha = loc_cha.split('-')
            new_SEED_name = SEED_name_template.format(net=net,
                                                sta=sta,
                                                loc=loc,
                                                cha=cha,
                                                year=year,
                                                jday=jday)
            # del_SEED_path = join(prepath,year,net,sta,cha,new_SEED_name+"_new")
            new_SEED_path = join(prepath,year,net,sta,cha,new_SEED_name)

            st.write(new_SEED_path,format='MSEED')
            # os.remove(del_SEED_path)


if __name__ == "__main__":
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
