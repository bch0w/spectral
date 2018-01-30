"""24/1/18
Use geonet permanent station data for noise spectra analysis
in accordance with McNamara2004 using the obspy PPSD package
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

# manual set time interval
start_date = UTCDateTime('2015-01-01')
end_date = UTCDateTime('2016-01-01')
print("Timeframe set: {} through {}".format(start_date.date,end_date.date))

# set instrument choice on command line
station = sys.argv[1]
if station[-1].upper() == 'Z':
    channel = 'HH' + sys.argv[2]
elif station[-1].upper() == 'S':
    channel = 'BN' + sys.argv[2]
instrument_id = 'NZ.{}.*.{}'.format(station,channel)
print("Analyzing noise spectra for {}".format(instrument_id))

# gather all data files and response file
net, sta, loc, cha = instrument_id.split('.')
filepath = os.path.join('.',net,sta,cha)
resppath = os.path.join('.',net,sta,'*.xml')
data_files = glob.glob(filepath+'/*')
if not bool(data_files):
    sys.exit('No data available')
data_files.sort()
resp_file = glob.glob(resppath)[0]

# parse list and kick out time outliers
outliers = []
for filename in data_files:
    year,jday = os.path.basename(filename).split('.')[4:6]
    filedate = UTCDateTime("{}-{}".format(year,jday))
    if not start_date <= filedate <= end_date:
        outliers.append(filename)

for ol in outliers:
    print("Removing {}... outside timeframe".format(ol))
    data_files.remove(ol)

# determine start and end days
net,sta,loc,cha,year,jday_start,ext = os.path.basename(
                                                    data_files[0]).split('.')
jday_end = os.path.basename(data_files[-1]).split('.')[-2]
year_end = os.path.basename(data_files[-1]).split('.')[-3]

print(len(data_files)," data files found")
print("timespan: {year}.{day_start}-{year_end}.{day_end}".format(year=year,
                                                    day_start=jday_start,
                                                    year_end=year_end,
                                                    day_end=jday_end))
print("response file: ", resp_file)

# isolate latest response
inv = read_inventory(resp_file, format="STATIONXML")
resp = inv[0].select(channel=cha)[0][-1].response
decimateby = 5
# poles_and_zeros = resp.get_paz()


# initialize PPSD with first datafile then remove it from list
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

# loop over all datafiles and add to ppsd (takes time)
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

output_filename = '{sta}.{cha}.{year}.{day_start}-{day_end}'.format(sta=sta,
                                                        cha=cha,
                                                        year=year,
                                                        day_start=jday_start,
                                                        day_end=jday_end)

# save data as .npz array, save plot as .png
npz_filepath = './ppsd_arrays/{}.npz'.format(output_filename)
png_filepath = './figures/{}.png'.format(output_filename)

# check npz file exists
if os.path.exists(npz_filepath):
    npz_check = input('npz path exists, overwrite? [y]/n')
    if npz_check == 'n':
        ppsd.save_npz(npz_filepath)
else:
    ppsd.save_npz(npz_filepath)

# check png file exists
if os.path.exists(png_filepath):
    png_check = input('png path exists, overwrite? [y]/n')
    if png_check == 'n':
        import pdb;pdb.set_trace()
        sys.exit('Aborted')

ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True)
