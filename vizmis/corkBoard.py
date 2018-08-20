"""adjointBuilder saves window information into pyAsdf data format. Cork class
contains functions to parse through this data format and provide easily visible 
information on best windows, number of windows per station, etc. as well as
printing and sorting functionalities to make it easier to interact with 
all the available data
"""
import os
import sys
import pyasdf
import numpy as np
sys.path.append('../modules')
from getdata import pathnames
from procmod import myround

def distribute_to_corkBoard(event_id):
    """initiate a Cork object with the information given its event id
    """
    mycork = Cork(event_id)
    import ipdb;ipdb.set_trace()
    mycork.aggregate()
    
class Cork:
    """a class used to hold information related to source-receiver pairs, 
    misfit windows and adjoint sources
    """
    
    def __init__(self,event_id=None):
        """fill up the Cork with empties
        """
        self.id = event_id
        if self.id:
            self.fid = event_id + '.h5'
        else:
            self.fid = None
        self.path = None
        self.ds = None
        self.event = None
        self.stations = None
        self.adj_srcs = None
        self.misfit_windows = None
        self.misfit_values = None
        
    # def __len__(self):
                    
    def __str__(self):
        """replace print() output
        """
        template = ("corkBoard for {pat}\n"
                    "\t{nos:<4} stations\n"
                    "\t{now:<4} misfit_windows\n"
                    "\t{noa:<4} adj_srcs\n")
        return template.format(pat=self.path,
                               nos=len(self.stations),
                               now=len(self.misfit_windows),
                               noa=len(self.adj_srcs)
                               )
        
    def _check_availability(self):
        """will determine if the data is available
        """
        # local pathing, to be changed if this is used outside 
        data_location = pathnames()["adjtomodata"] + "PYASDF"

        # if an event id is given, give filepath
        if not self.fid:
            import glob
            available = glob.glob(os.path.join(data_location,"*.h5"))
            print("No event id given, available event ids: ")
            for tag in available:
                print(os.path.basename(tag).split('.')[0])
            event_id = input("Event ID: ")
            self.id = event_id
            self.fid = event_id + '.h5'
        
        filepath = os.path.join(data_location,self.fid)
        if os.path.exists(filepath):
            self.path = filepath
        else:
            raise Exception("File does not exist")


    def _read(self):
        """read in datafile
        """
        self._check_availability()
        try:            
            ds = pyasdf.ASDFDataSet(self.path)
        except OSError:
            raise Exception("{} is currently open, please close".format(
                                                                self.filename))
        return ds
    
    def aggregate(self):
        """fills up an initiated cork using all available internal functions
        """
        self.populate()
        self.count_windows()
        self.collect_misfits()
        self.get_srcrcv_information()
        self.collect_misfits()
        
    def populate(self):
        """fills self with all information available in pyasdf dataset
        """
        ds = self._read()
        self.ds = ds
        self.stations = ds.waveforms.list()
        self.aux = ds.auxiliary_data
        self.adj_srcs = ds.auxiliary_data.AdjointSource.list()
        self.misfit_windows = ds.auxiliary_data.MisfitWindows.list()
        
    def count_windows(self,print=False):
        """figure out which stations contain which windows, return a dictionary
        which lists available components. should be run within populate cork
        """
        stations = []
        for win in self.misfit_windows:
            stations.append(
                self.aux.MisfitWindows[win].parameters['channel_id'])
        uniqueid = set(stations)
        counted_windows = {}
        for id in uniqueid:
            counted_windows[id] = stations.count(id)
            
        self.counted_windows = counted_windows
        
    def collect_misfits(self,model="m00"):
        """for each station, collect the misfit value
        """
        misfit_values = {}
        for AS in self.adj_srcs:
            parm = self.aux.AdjointSource[AS].parameters
            channel_id = '{}.{}'.format(parm["station_id"],parm["component"])
            misfit_values[channel_id] = parm["misfit_value"]
        
        self.misfit_values = misfit_values
            
    def get_srcrcv_information(self):
        """calculate source receiver information for each pair
        """
        from obspy.geodetics import gps2dist_azimuth

        event_lat,event_lon = (self.ds.events[0].origins[0].latitude,
                               self.ds.events[0].origins[0].longitude)
        
        srcrcvdict = {}
        for sta in self.stations:
            coordict = self.ds.waveforms[sta].coordinates
            GCDist,Az,BAz = gps2dist_azimuth(lat1=event_lat,lon1=event_lon,
                                             lat2=coordict["latitude"],
                                             lon2=coordict["longitude"])
            srcrcvdict[sta] = {"great_circle_distance":GCDist,
                               "azimuth":Az,
                               "backazimuth":BAz}
        
        self.source_receiver_info = srcrcvdict
        
    def misfit_histogram(self,m0=None,m_a=None,binsize=0.1):
        """plot a histogram of misfits for models m0 and m_a (if available)
        default will be initial model compared to final model
        """
        if self.misfit_values:
            raise Exception("No misfit values available")
        
        # only import plotting functions if necessary
        import matplotlib.pyplot as plt
        
        # gather misfit values
        # np.digitize returns the bin number that values of x fall into
        misfits = np.fromiter(self.misfit_values.values(),dtype="float")
        maxmisfit = myround(misfits.max(),base=1,choice="up")

        n,bins,patches = plt.hist(x=misfits,
                                  bins=len(np.arange(0,maxmisfit,binsize)),
                                  range=(0,maxmisfit))
    
    # def disperse(self):
    #     """data vomit all available data in the Cork object for easy viewing
    #     """
    #     print("corkBoard for {id}".format(self.path)
    #     print("{s:^20}{w:^20}{m:^20}".format(
    #                                      s="STATION",w="WINDOWS",m="MISFIT"))
    #     for sta in self.stations:
    # 
    #         print("{s:^20}")

if __name__ == "__main__":
    distribute_to_corkBoard("2014p240655")
            
        
    

        
        

        
        
        
        
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
            
        
    
