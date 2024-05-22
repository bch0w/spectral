"""
Gather metadata for IMS stations for reference and maps
"""
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


# From Gibbons 2015 SRL (current IMS stations)
# IMS Primary Seismic Arrays
# ILAR, YKA, PDAR, NVAR, TXAR, ARCES, NOA, FINES, GERES, ESDC, TORD,
# AKASG, BRTR, GEYT, ZALV, MKAR, SONM, USRK, PETK, KSRS, MJAR, CMAR, WRA, ASAR
ims_primary_arrays = [
        "IM.ILAR",   # ILAR Array Beam, Eilson, AK, USA
        "IM.YKA",   # Yellowknife Array Beam
        "IM.PDAR",   # Pinedale Array Beam, Wyoming, U.S.A.
        "IM.NVAR",   # Mina Array Beam, Nevada, U.S.A.
        "IM.TXAR",   # Lajitas Array Beam, Texas, U.S.A.
        "IM.ARCES",   # ARCESS Array Beam, Norway
        "IM.FINES",   # FINESS Array Beam, Finland
        "IM.GERES",   # GERESS Array Beam, Bayern, Germany
        "IM.ESDC",   # Sonseca Array Beam, Spain
        "IM.TORD",   # Torodi Array Beam, Niger
        "IM.AKASG",   # Malin Array Beam, Ukraine
        "SY.AKASG",   # AKASG synthetic
        "IM.GEYT",   # Alibeck Array Beam, Turkmenistan
        "IM.ZALV",   # Zalesovo Array Beam, Altayskiy Kray, Russia
        "IM.MKAR",   # Makanchi Array Beam, Kazakhstan
        "KZ.MKAR",   # Makanchi array,MK31, Kazakhstan
        "NM.MKAR",   # Mark Tree,AR(CERI)
        "IM.SONM",   # Songino Array Beam, Mongolia
        "IM.USRK",   # Ussuriysk Array Beam, Primorskiy Kray, Russia
        "IM.PETK",   # Petropavlovsk-Kamchatskiy Array Beam, Kamchatskaya Oblast, R
        "IM.MJAR",   # Matsushiro Array Beam, Nagano, Japan
        "IM.CMAR",   # Chiang Mai Array Beam, Thailand
        "IM.WRA",   # Warramunga Array Beam
        "IM.ASAR",   # Alice Springs Array Beam
        "NM.ASAR",   # Arkansas State,Jonesboro,AR(CERI)
        ]


# IMS Auxiliary Seismic Arrays
# SPITS, EKA, HES,KVAR, MMAI, BVAR,KURK
ims_auxiliary_arrays = [
        "NC.HES",   # Elkhorn Slough
        "UU.HES",   # Hooper Elementary School, UT, USA
        "KZ.BVAR",   # Borovoye, Kazakstan
        "NM.BVAR",   # Bay Village,AR
        "II.KURK",   # Kurchatov, Kazakhstan
        "KZ.KURK",   # Kurchatov, Kazakhstan
        ]

# IMS Primary 3C Stations
# ULM, SCHQ, ROSC, LPAZ,BDFB, CPUP, PLCA, PPT, DBIC, BOSA, KMBO, STKA, MAW, VNDA
# KBZ
ims_primary_seismic = [
        "CN.ULM",   # Lac Du Bonnet, MB, CA
        "CN.SCHQ",   # Schefferville, QC, CA
        "CM.ROSC",   # El Rosal, Cundinamarca, Colombia
        "GT.LPAZ",   # La Paz, Bolivia
        "GT.BDFB",   # Brasilia, Distrito Federal, Brazil
        "GT.CPUP",   # Villa Florida, Paraguay
        "GT.PLCA",   # Paso Flores, Argentina
        "GT.PLCA",   # Paso Flores, Argentina
        "G.PPT",   # Pamatai - Papeete - Tahiti island - French Polynesia, France
        "NC.PPT",   # Peach Tree Valley
        "GT.DBIC",   # Dimbokro, Cote d'Ivoire
        "GT.DBIC",   # Dimbokro, Ivory Coast
        "GT.BOSA",   # Boshof, South Africa
        "GE.KMBO",   # IRIS/GEOFON Station Kilima Mbogo, Kenya
        "IU.KMBO",   # Kilima Mbogo, Kenya
        "AU.STKA",   # Stephens Creek, NSW
        "AU.MAW",   # Mawson, Antarctica
        "GT.VNDA",   # Dry Valley, Vanda, Antarctica
        "IM.VNDA",   # Vanda, Antarctica
        ]


# Get WATC Stations first
for code in watc_seismic_stations:
    net, sta = code.split(".")
    inv = c.get_stations(network=net, station=sta, channel="BH?,HH?", 
                         level="channel")
    


