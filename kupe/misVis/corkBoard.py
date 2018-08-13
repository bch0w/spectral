"""adjointBuilder saves window information into pyAsdf data format. this class
contains functions to parse through this data format and provide easily visible 
information on best windows, number of windows per station, etc.
"""
import os
import sys
import pyasdf
sys.path.append('../../modules')
from getdata import pathnames

def populate(event_id):
    """initiate a corkBoard object with the information given its event id
    """
    

class Cork:
    """a class used to hold information related to source-receiver pairs
    """
    # local pathing, to be changed if this is used outside my own directories
    data_location = pathnames["kupedata"] + "PYASDF"
    
    def __init__(self,event_id):
        self.id = event_id
        self.fid = event_id + '.h5'
        self.stations = []
        
        
    def _check_availability(self):
        """will determine if the data is available
        """
        filename = os.path.join(data_location,self.fid)
        if os.path.exists(filename):
            self.filename = filename
        else:
            raise Exception("File does not exist")
    
    def _read(self):
        """read in datafile
        """
        self._check_availability()            
        ds = pyasdf.ASDFDataSet(self.filename)
        
        return ds
    
    def _populate_cork(self):
        """fills self with all information available in pyasdf dataset
        """
        ds = self._read()
        self.event = ds.events[0]
        
        
    def __str__():
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
            
        
    