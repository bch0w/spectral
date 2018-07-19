import glob
import os
from obspy import read
from scipy.io import savemat


def this_doesnt_work():
	"""copied from obspy documentation by tr.stats.iteritems() doesn't work
	"""
	outpath = '/seis/prj/fwi/bchow/mseeds/ASCII/'
	path = '/seis/prj/fwi/bchow/mseeds/ASCII/*.ascii'
	allfiles = glob.glob(path)
	for fid in allfiles:
		st = read(fid)
		filename = os.path.basename(fid)
		for i, tr in enumerate(st):
			mdict = {k: str(v) for k,v in tr.stats.iteritems()}
			mdict['data'] = tr.data
			outname = os.path.join(outpath,filename,'.mat')
			savemat(outname,mdict)

def stream_to_2column_ascii():
	outpath = '/seis/prj/fwi/bchow/mseeds/ASCII/'
	path = '/seis/prj/fwi/bchow/mseeds/ASCII/*KU15-3*.ascii'
	allfiles = glob.glob(path)
	for fid in allfiles:
		st = read(fid)
		filename = os.path.basename(fid)
		for tr in st:
			filename = filename.split('.')[0]
			newfid = "{}".format(filename)
			fidout = os.path.join(outpath,newfid)
			dt = tr.stats.delta
			import ipdb;ipdb.set_trace()
			with open(fidout,'w') as f:
				for i,data in enumerate(tr.data):
					f.write('{:.2f}\t{}\n'.format(i*dt,data))


if __name__ == "__main__":
	stream_to_2column_ascii()
