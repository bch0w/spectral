import os
from obspy import UTCDateTime


def read_write_around_chiapas(sta, plusminus=2):
    """
    read in raw seismic data and response, remove response and save for later
    """
    pathdir = ("/Users/chowbr/Documents/subduction/mseeds/BEACON/2017/XX/"
               "{sta}/HH{comp}.D/XX.{sta}.10.HH{comp}.D.2017.{jday}")
    chiapas = UTCDateTime("2017-09-08T00:00:00")
    inv = read_inventory("/Users/chowbr/Documents/subduction/mseeds/"
                         "BEACON/DATALESS/XX.RDF.DATALESS")
    for comp in ["N", "E"]:
        for jday in range(chiapas.julday-plusminus, chiapas.julday+plusminus):
            fid_in = pathdir.format(sta=sta, comp=comp, jday=jday)
            st = read(fid_in)
            st.attach_response(inv)
            st.remove_response(output="VEL", water_level=60, plot=False)
            fid_out = ("/Users/chowbr/Documents/subduction/spectral/tremor/"
                       "gsnz18/mseed_remove_response/{fid_in}")
            st.write(fid_out.format(
                fid_in=os.path.basename(fid_in)), format="mseed")


def download_geonet_data(sta, plusminus=3):
    from obspy.clients.fdsn import Client
    c = Client("GEONET")
    chiapas = UTCDateTime("2017-09-08T00:00:00")
    fid_template = "{id}.D.{year}.{jday}"
    for i in range(-plusminus, plusminus+1, 1):
        starttime = chiapas + i * 24*3600
        endtime = starttime + 24 * 3600
        for comp in ["E", "N"]:
            try:
                st = c.get_waveforms(network="NZ", station=sta, location="*",
                                     channel="?H{}".format(comp),
                                     starttime=starttime,
                                     endtime=endtime, attach_response=True
                                     )
                fid_out = ("/Users/chowbr/Documents/subduction/mseeds/GEONET/"
                           + fid_template.format(id=st[0].get_id(),
                                                 year=starttime.year,
                                                 jday=starttime.julday
                                                 )
                           )
                if os.path.exists(fid_out):
                    continue
                st.remove_response(output="VEL", water_level=60, plot=False)
                st.write(fid_out, format="MSEED")
                print("wrote"+fid_template.format(
                    id=st[0].get_id(), year=starttime.year, jday=starttime.julday)
                      )
            except Exception as e:
                continue


if __name__ == "__main__":
    for sta in ["BFZ"]:
        download_geonet_data(sta, plusminus=5)