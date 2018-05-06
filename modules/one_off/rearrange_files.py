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

def VUW_fixo():
    """
    RD01.2017.319.10-HHE.XX.ms RD08.2017.264.10-HHE.XX.ms RD14.2017.321.10-HHE.XX.ms
    RD01.2017.319.10-HHN.XX.ms RD08.2017.264.10-HHN.XX.ms RD14.2017.321.10-HHN.XX.ms
    RD01.2017.319.10-HHZ.XX.ms RD08.2017.264.10-HHZ.XX.ms RD14.2017.321.10-HHZ.XX.ms
    RD03.2017.322.10-HHE.XX.ms RD09.2017.320.10-HHE.XX.ms RD15.2017.263.10-HHE.XX.ms
    RD03.2017.322.10-HHN.XX.ms RD09.2017.320.10-HHN.XX.ms RD15.2017.263.10-HHN.XX.ms
    RD03.2017.322.10-HHZ.XX.ms RD09.2017.320.10-HHZ.XX.ms RD15.2017.263.10-HHZ.XX.ms
    RD04.2017.322.10-HHE.XX.ms RD10.2017.321.10-HHE.XX.ms RD16.2017.263.10-HHE.XX.ms
    RD04.2017.322.10-HHN.XX.ms RD10.2017.321.10-HHN.XX.ms RD16.2017.263.10-HHN.XX.ms
    RD04.2017.322.10-HHZ.XX.ms RD10.2017.321.10-HHZ.XX.ms RD16.2017.263.10-HHZ.XX.ms
    RD05.2017.323.10-HHE.XX.ms RD11.2017.264.10-HHE.XX.ms RD17.2017.322.10-HHE.XX.ms
    RD05.2017.323.10-HHN.XX.ms RD11.2017.264.10-HHN.XX.ms RD17.2017.322.10-HHN.XX.ms
    RD05.2017.323.10-HHZ.XX.ms RD11.2017.264.10-HHZ.XX.ms RD17.2017.322.10-HHZ.XX.ms
    RD06.2017.320.10-HHE.XX.ms RD12.2017.320.10-HHE.XX.ms RD18.2017.322.10-HHE.XX.ms
    RD06.2017.320.10-HHN.XX.ms RD12.2017.320.10-HHN.XX.ms RD18.2017.322.10-HHN.XX.ms
    RD06.2017.320.10-HHZ.XX.ms RD12.2017.320.10-HHZ.XX.ms RD18.2017.322.10-HHZ.XX.ms
    RD07.2017.264.10-HHE.XX.ms RD13.2017.264.10-HHE.XX.ms 
    RD07.2017.264.10-HHN.XX.ms RD13.2017.264.10-HHN.XX.ms
    RD07.2017.264.10-HHZ.XX.ms RD13.2017.264.10-HHZ.XX.ms
    overlap list
    """
    from os.path import join
    import os
    import glob
    from obspy import read
    list1 = '/Users/chowbr/Documents/subduction/spectral/common/DATA/RDF_Array/Sep2017_Nov2017/DATA_ALL'
    list2 = '/Users/chowbr/Documents/subduction/RDF'
    overlap_list = glob.glob(list1 + '/*.ms')
    SEED_name_template = "{net}.{sta}.{loc}.{cha}.{year}.{jday}"
    for entry in overlap_list:
        file1 = glob.glob(join(list1,entry))[0]
        fid = os.path.basename(file1)
        file2 = glob.glob(join(list2,fid))[0]
        st = read(file1) + read(file2)
        st.write(file2,format='MSEED')
        # os.remove(del_SEED_path)

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
            
def was_in_main():
    # cwd = "/seis/prj/fwi/bchow/RDF/"
    cwd = '/Users/chowbr/Documents/subduction/RDF/'
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
        new_SEED_basepath = join(prepath,year,net,sta,cha)
        new_SEED_path = join(prepath,year,net,sta,cha,new_SEED_name)
        if not os.path.exists(new_SEED_basepath):
            os.mkdir(new_SEED_basepath)
        os.rename(src=old_filename,dst=new_SEED_path)
        
def rename_files_and_folders():
    """adding .D to all files and folders in RDF folder
    """
    cwd = "/seis/prj/fwi/bchow/RDF/"
    folders_ = False
    files_ = True
    
    # cwd = '/Users/chowbr/Documents/subduction/RDF/'
    os.chdir(cwd)
    if folders_:
        folders_to_rename = glob.glob(join(cwd,"*","*","*","HH?"))
        for old_folder in folders_to_rename:
            new_folder = old_folder + '.D'
            os.rename(src=old_folder,dst=new_folder)
    
    if files_: 
        files_to_rename = glob.glob(join(cwd,"*","*","*","*","*"))
        for old_file in files_to_rename:
            new_seed_format = "{N}.{S}.{L}.{C}.D.{Y}.{D}"
            
            prepath = os.path.dirname(old_file)
            net,sta,loc,chan,year,day = os.path.basename(old_file).split('.')
            new_file_base = new_seed_format.format(
                                        N=net,S=sta,L=loc,C=chan,Y=year,D=day)

            new_file = join(prepath,new_file_base)
            
            os.rename(src=old_file,dst=new_file)
            
def rearrange_SAHKE():
    """converting SAHKE data from SAC to mseed and rename/organize files
    """
    import os
    import glob
    import numpy as np
    from obspy import read
    
    sacpath = '/Users/chowbr/Documents/subduction/SAHKE/SAC/'
    sacfolders = glob.glob(sacpath + '*')
    
    name_template = "X2.{s}..{c}.D.{y}.{j}"
    savepath_template = ('/Users/chowbr/Documents/subduction/SAHKE/2010/X2/{s}/'
                         '{c}.D/')
    for fold in sacfolders:
        os.chdir(fold)
        sacfiles = glob.glob('*')
        for fid in sacfiles:
            basename = os.path.basename(fid)
            try:
                net,sta,cha,year,jday,fmt = basename.split('.')
            except ValueError:
                sta,year,jday,comp,fmt = basename.split('.')
                cha = "HH{}".format(comp)
                net = "X2"
            
            newbasename = name_template.format(s=sta,c=cha,y=year,j=jday)
            savepath = savepath_template.format(s=sta,c=cha)
            savename = os.path.join(savepath,newbasename)
            if os.path.exists(savename):
                continue
            try:
                st = read(fid)
            except TypeError:
                print(fid)
                continue
            st[0].data = st[0].data.astype(np.int32)
            st.write(savename,format='MSEED',reclen=512,encoding='STEIM1')
        os.chdir(sacpath)
                                


if __name__ == "__main__":
    rearrange_SAHKE()