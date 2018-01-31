"""24/1/18
Use geonet permanent station data for noise spectra analysis
in accordance with McNamara2004 using the obspy PPSD package
!!! currently only works for year-long spectral analyses, grabs entire year
!!! directory from geonet archive for analysis. which year can be specified
"""
import os
import sys
import glob
import time
from obspy import read
from obspy import read_inventory
from obspy.signal import PPSD
from obspy import UTCDateTime
from obspy.imaging.cm import pqlx
from getdata_fdsn import geonet_data

# manual set time interval
start_date = UTCDateTime('2015-01-01')
end_date = UTCDateTime('2016-01-01')
print("Timeframe set: {} through {}".format(start_date.date,end_date.date))

# set instrument choice on command line
station = sys.argv[1].upper()
if station[-1].upper() == 'Z':
    channel = 'HH' + sys.argv[2].upper()
elif station[-1].upper() == 'S':
    channel = 'BN' + sys.argv[2].upper()
instrument_id = 'NZ.{}.*.{}'.format(station,channel)
print("Analyzing noise spectra for {}".format(instrument_id))

# filepaths of geonet archive
net, sta, loc, cha = instrument_id.split('.')
mseed_GNApath = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(
                                                        year=start_date.year,
                                                        sta=sta,
                                                        ch=cha)
resp_GNApath = '/geonet/seed/RESPONSE/{}.NZ/'.format(sta)

# data files
data_files = glob.glob(mseed_GNApath+'*')
data_files.sort()

# response file
NET,STA,LOC,CHA,SUFX,YEAR,JDAY = os.path.basename(data_files[0]).split('.')
response_filename = "RESP.{net}.{sta}.{loc}.{cha}".format(net=NET,
                                                            sta=STA,
                                                            loc=LOC,
                                                            cha=CHA)
resp_file = os.path.join(resp_GNApath,response_filename)

# notify files found
print("++ data file(s): ", mseed_GNApath)
print("\t {} data files found".format(len(data_files)))
print("++ response file: ", resp_file)

# read in response file, set decimate parameter
inv = read_inventory(resp_file)
decimateby = 5

# initialize PPSD with first datafile
print("1/{} Initializing with data file: ".format(len(data_files)),
                                                    data_files[0],end='... ')
start = time.time()
st = read(data_files[0])
st.decimate(decimateby)
tr = st.select(channel=cha)[0]
ppsd = PPSD(tr.stats, metadata=inv)
ppsd.add(st)
end = time.time()
print("complete ({}s)".format(round(end-start,2)))

# loop over rest of datafiles and add to ppsd
for i,filename in enumerate(data_files[1:]):
    print('{0}/{1} {2}'.format(i+2,len(data_files),filename),end='... ')
    try:
        start = time.time()
        st = read(filename)
        st.decimate(decimateby)
        ppsd.add(st)
        end = time.time()
        print("complete ({}s)".format(round(end-start,2)))
    except Exception as e:
        print(e)
        pass

# set output name for saving .npz and f# save data as .npz array, save plot as .png
# i'm assuming here the year stays the same
year_start,jday_start = os.path.basename(data_files[0]).split('.')[-2:]
year_end,jday_end = os.path.basename(data_files[-1]).split('.')[-2:]

output_filename = '{sta}.{cha}.{year}.{day_start}-{day_end}'.format(sta=sta,
                                                        cha=cha,
                                                        year=year_start,
                                                        day_start=jday_start,
                                                        day_end=jday_end)

npz_filepath = './ppsd_arrays/{}.npz'.format(output_filename)
png_filepath = './figures/ppsd_plots/{}.png'.format(output_filename)

# check npz file exists
if os.path.exists(npz_filepath):
    npz_check = input('npz path exists, overwrite? [y]/n')
    if npz_check == 'n':
        ppsd.save_npz(npz_filepath)
else:
    ppsd.save_npz(npz_filepath)
    print("saved .npz file: ",npz_filepath)

# check png file exists
if os.path.exists(png_filepath):
    png_check = input('png path exists, overwrite? [y]/n')
    if png_check == 'n':
        import pdb;pdb.set_trace()
        sys.exit('Aborted')

ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True)
print("saved .png file: ",png_filepath)
