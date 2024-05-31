"""
Gather metadata for IMS stations for reference, maps and simulation work
"""
from obspy import Inventory
from obspy.clients.fdsn import Client


c = Client("IRIS")

# https://watc.alaska.edu/operations-and-maintenance (last access 5-22-24)
watc_seismic_stations = [
        "IU.ANMO",  # FACT, NM
        "US.ELK",   # AS113: Elko, NV
        "IU.GUMO",  # Guam, Mariana Islands,
        "II.KDAK",  # Kodiak, AK
        "US.NEW",   # Newport, WA
        "IM.NVAR",  # Mina, NV (12-station array)
        "II.PFO",   # Pinyon Flats, CA
        "IU.PMSA",  # Palmer Station, Antarctica
        "IU.QSPA",  # AS114: South Pole Station, Antarctica
        "IM.SHEM",  # Shemya Island, AK (2x stations)
        "IU.SJG",   # AS116: San Juan, Puerto Rico
        "IM.TKL",   # Tuckaleechee Caverns, TN
        "BK.YBH",   # Yreka Blue Horn Mine, CA
        ]

# IMS Primary Seismic
ims_primary_seismic = [
        "AU.MAW",   # Mawson, Antarctica
        "AU.STKA",   # Stephens Creek, NSW
        "CM.ROSC",   # El Rosal, Cundinamarca, Colombia
        "CN.SCHQ",   # Schefferville, QC, CA
        "CN.ULM",   # Lac Du Bonnet, MB, CA
        "G.PPT",   # Pamatai - Papeete - Tahiti island - French Polynesia, France
        "GE.KMBO",   # IRIS/GEOFON Station Kilima Mbogo, Kenya
        "GT.BDFB",   # Brasilia, Distrito Federal, Brazil
        "GT.BOSA",   # Boshof, South Africa
        "GT.CPUP",   # Villa Florida, Paraguay
        "GT.DBIC",   # Dimbokro, Cote d'Ivoire
        "GT.DBIC",   # Dimbokro, Ivory Coast
        "GT.LPAZ",   # La Paz, Bolivia
        "GT.PLCA",   # Paso Flores, Argentina
        "GT.PLCA",   # Paso Flores, Argentina
        "GT.VNDA",   # Dry Valley, Vanda, Antarctica
        "IM.AKASG",   # Malin Array Beam, Ukraine
        "IM.ARCES",   # ARCESS Array Beam, Norway
        "IM.ASAR",   # Alice Springs Array Beam
        "IM.CMAR",   # Chiang Mai Array Beam, Thailand
        "IM.ESDC",   # Sonseca Array Beam, Spain
        "IM.FINES",   # FINESS Array Beam, Finland
        "IM.GERES",   # GERESS Array Beam, Bayern, Germany
        "IM.GEYT",   # Alibeck Array Beam, Turkmenistan
        "IM.ILAR",   # ILAR Array Beam, Eilson, AK, USA
        "IM.MJAR",   # Matsushiro Array Beam, Nagano, Japan
        "IM.MKAR",   # Makanchi Array Beam, Kazakhstan
        "IM.NVAR",   # Mina Array Beam, Nevada, U.S.A.
        "IM.PDAR",   # Pinedale Array Beam, Wyoming, U.S.A.
        "IM.PETK",   # Petropavlovsk-Kamchatskiy Array Beam, Kamchatskaya Oblast, R
        "IM.SONM",   # Songino Array Beam, Mongolia
        "IM.TORD",   # Torodi Array Beam, Niger
        "IM.TXAR",   # Lajitas Array Beam, Texas, U.S.A.
        "IM.USRK",   # Ussuriysk Array Beam, Primorskiy Kray, Russia
        "IM.VNDA",   # Vanda, Antarctica
        "IM.WRA",   # Warramunga Array Beam
        "IM.YKA",   # Yellowknife Array Beam
        "IM.ZALV",   # Zalesovo Array Beam, Altayskiy Kray, Russia
        "IU.KMBO",   # Kilima Mbogo, Kenya
        "NC.PPT",   # Peach Tree Valley
        ]

ims_auxiliary_seismic = [
        "AK.MLR",   # Meier's Lake Roadhouse, AK, USA
        "AU.FITZ",   # Fitzroy Crossing
        "AU.NWAO",   # Narrogin, Western Australia
        "AZ.PFO",   # Pinyon Flats Observatory, CA, USA
        "BK.YBH",   # Yreka Blue Horn Mine, Yreka, CA, USA
        "CN.DLBC",   # Dease Lake, BC, CA
        "CN.FRB",   # Iqaluit, NU, CA
        "CN.INK",   # Inuvik, NT, CA
        "CN.RES",   # Resolute, NU, CA
        "CN.SADO",   # Sadowa, ON, CA
        "CZ.VRAC",   # Vranov
        "DT.PMG",   # Port Moresby, Papua New Guinea
        "DW.AFI",   # Afiamalu, Samoa
        "DW.LEM",   # Lembang, Indonesia
        "G.ATD",   # Arta Cave - Arta, Republic of Djibouti
        "G.SEY",   # Seymchan, Russia
        "GE.EIL",   # GII/GEOFON Station Eilat, Israel
        "GE.SFJD",   # IRIS/GEOFON Station Sondre Stromfjord, Greenland
        "GE.SNAA",   # AWI/GEOFON Station Sanae, Antarctica
        "GS.NEW",   # Seattle urban array
        "ID.NNA",   # Nana, Peru
        "ID.RPN",   # Easter Island (Rapa Nui), Chile
        "ID.SJG",   # San Juan, Puerto Rico
        "ID.SUR",   # Sutherland, Republic of South Africa
        "II.AAK",   # Ala Archa, Kyrgyzstan
        "II.BORG",   # Borgarfjordur, Asbjarnarstadir, Iceland
        "II.JTS",   # Las Juntas de Abangares, Costa Rica
        "II.KDAK",   # Kodiak Island, Alaska, USA
        "II.KURK",   # Kurchatov, Kazakhstan - ARRAY
        "II.MBAR",   # Mbarara, Uganda
        "II.OBN",   # Obninsk, Russia
        "II.PALK",   # Pallekele, Sri Lanka
        "IM.TKL",   # Tuckaleechee Caverns, TN, USA
        "IU.ANMO",   # Albuquerque, New Mexico, USA
        "IU.GNI",   # Garni, Armenia
        "IU.GUMO",   # Guam, Mariana Islands
        "IU.HNR",   # Honiara, Solomon Islands
        "IU.MSKU",   # Masuku, Gabon
        "IU.PMSA",   # Palmer Station, Antarctica
        "IU.QSPA",   # South Pole Remote Earth Science Observatory (Quiet Zone)
        "IU.RAO",   # Raoul, Kermadec Islands
        "IU.RCBR",   # Riachuelo, Brazil
        "IU.TEIG",   # Tepich, Yucatan, Mexico
        "IU.TSUM",   # Tsumeb, Namibia
        "JP.JKA",   # Kamikawa Asahi
        "JP.JNU",   # Oita Nakatsue
        "JP.JOW",   # Okinawa Kunigami
        "KZ.AKTO",   # Aktyubinsk, Kazakhstan
        "KZ.BVAR",   # Borovoye, Kazakstan - ARRAY
        "MN.IDI",   # Anogia, Greece
        "MN.MDT",   # Midelt, Morocco
        "MN.VAE",   # Valguarnera, Italy
        "NC.HES",   # Elkhorn Slough - ARRAY
        "NP.ELK",   # Elk Grove
        "NZ.RPZ",   # Rata Peaks
        "NZ.URZ",   # Urewera
        "PS.TGY",   # Tagaytay, Phillipines
        "SY.USHA",   # USHA synthetic
        "VE.PCRV",   # Pto. La Cruz - Anzoategui - CTBTO borehole
        "XY.BBTS",   # Black Bottom Tack and Stables - Summerville, SC
        ]

# From Brown et al. 2014
# ims_primary_seismic = [['IM.AKASG', 'IM.ARCES', 'IM.ASAR', 'GT.BDFB', 'GT.BOSA', 
#                         'IM.CMAR', 'GT.CPUP', 'GT.DBIC', 'IM.ESDC', 'IM.FINES', 
#                         'IM.GERES', 'IM.GEYT', 'IM.ILAR', 'GE.KMBO', 'GT.LPAZ', 
#                         'AU.MAW', 'IM.MJAR', 'IM.MKAR', 'IM.NVAR', 'IM.PDAR', 
#                         'IM.PETK', 'GT.PLCA', 'G.PPT', 'CM.ROSC', 'CN.SCHQ', 
#                         'IM.SONM', 'AU.STKA', 'IM.TORD', 'IM.TXAR', 'CN.ULM', 
#                         'IM.USRK', 'IM.WRA', 'IM.YKA', 'IM.ZALV']
# ims_auxiliary_seismic = ['II.AAK', 'DW.AFI', 'KZ.AKTO', 'IU.ANMO', 'G.ATD', 
#                          'XY.BBTS', 'II.BORG', 'KZ.BVAR', 'CN.DLBC', 'GE.EIL', 
#                          'NP.ELK', 'AU.FITZ', 'CN.FRB', 'IU.GNI', 'IU.GUMO', 
#                          'IU.HNR', 'MN.IDI', 'CN.INK', 'JP.JKA', 'JP.JNU', 
#                          'JP.JOW', 'II.JTS', 'II.KDAK', 'II.KURK', 'DW.LEM', 
#                          'II.MBAR', 'MN.MDT', 'AK.MLR', 'IU.MSKU', 'GS.NEW', 
#                          'ID.NNA', 'AU.NWAO', 'II.OBN', 'II.PALK', 'VE.PCRV', 
#                          'AZ.PFO', 'DT.PMG', 'IU.PMSA', 'IU.QSPA', 'IU.RAO', 
#                          'IU.RCBR', 'CN.RES', 'ID.RPN', 'NZ.RPZ', 'CN.SADO', 
#                          'G.SEY', 'GE.SFJD', 'ID.SJG', 'GE.SNAA', 'ID.SUR', 
#                          'IU.TEIG', 'PS.TGY', 'IM.TKL', 'IU.TSUM', 'NZ.URZ', 
#                          'SY.USHA', 'MN.VAE', 'CZ.VRAC', 'BK.YBH']
# Gather into single 
to_plot = watc_seismic_stations + ims_primary_seismic + ims_primary_seismic

# Get WATC Stations first
inv = Inventory()
for code in to_plot:
    net, sta = code.split(".")
    try:        
        inv_single = c.get_stations(network=net, station=sta, 
                                    channel="BH?,HH?,SH?", level="station"
                                    )
    except Exception as e:
        print(f"Error with {code}")
        continue
    inv.extend(inv_single)

inv.write("inv.xml", format="stationxml")
inv.plot()
    
    


