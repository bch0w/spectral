import os
import glob
import matplotlib.pyplot as plt
from obspy import read

files_NEW = glob.glob('OUTPUT_FILES_ORIGINAL/*')
files_ORI_path = ('OUTPUT_FILES_SRTM30P')
for fid_NEW in files_NEW:
	try:
		fid_NEW_base = os.path.basename(fid_NEW)
		fid_ORI = os.path.join(files_ORI_path,fid_NEW_base)
		for f_,c,l in zip([fid_NEW,fid_ORI],['r','k'],['NEW','ORIGINAL']):
			st = read(f_)
			st.filter('bandpass',freqmin=1/30,freqmax=1/6)
			plt.plot(st[0].data,c,label=l)
		plt.title(fid_NEW_base)
		plt.legend()
		plt.show()
	except Exception as e:
		plt.close()
		continue
