"""
Gather metadata for IMS stations for reference, maps and simulation work
"""
from obspy import Inventory
from obspy.clients.fdsn import Client
from pysep.utils.io import write_stations_file


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
        "CB.LZH",  # Lanzhou, Gansu Province, China
        "CM.ROSC",   # El Rosal, Cundinamarca, Colombia
        "CN.SCHQ",   # Schefferville, QC, CA
        "CN.ULM",   # Lac Du Bonnet, MB, CA
        "CS.KNG",  # Kinjal, CIS (Caucuses Array, Georgia)
        "G.PPT",   # Pamatai - Papeete - Tahiti island - French Polynesia, France
        "GE.KMBO",   # IRIS/GEOFON Station Kilima Mbogo, Kenya
        "GT.BDFB",   # Brasilia, Distrito Federal, Brazil
        #ERROR "GT.BGCA",   # Bogoin, Central African Republic
        "GT.BOSA",   # Boshof, South Africa
        "GT.CPUP",   # Villa Florida, Paraguay
        "GT.DBIC",   # Dimbokro, Cote d'Ivoire
        "GT.DBIC",   # Dimbokro, Ivory Coast
        "GT.LPAZ",   # La Paz, Bolivia
        "GT.PLCA",   # Paso Flores, Argentina
        "GT.PLCA",   # Paso Flores, Argentina
        "GT.VNDA",   # Dry Valley, Vanda, Antarctica
        "II.NIL",    # Nilore, Pakistan
        "II.NRIL",   # Norilsk, Russia
        # ERROR "IL.IR3",    # Iran LP Array
        "IM.AKASG",   # Malin Array Beam, Ukraine
        "IM.ARCES",   # ARCESS Array Beam, Norway
        "IM.ASAR",   # Alice Springs Array Beam
        "IM.CMAR",   # Chiang Mai Array Beam, Thailand
        "IM.ESDC",   # Sonseca Array Beam, Spain
        "IM.FINES",   # FINESS Array Beam, Finland
        "IM.GERES",   # GERESS Array Beam, Bayern, Germany
        "IM.GEYT",   # Alibeck Array Beam, Turkmenistan
        "IM.ILAR",   # ILAR Array Beam, Eilson, AK, USA
        "IM.KSAR",   # Wonju Array Beam, South Korea
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
        #ERROR "YV.AT19",  # Kursunkaya, Kirikkale, Turkey
        ]

ims_auxiliary_seismic = [
        "AK.MLR",	# Meier's Lake Roadhouse, AK, USA
        "AU.FITZ",	# Fitzroy Crossing
        "AU.NWAO",	# Narrogin, Western Australia
        "AZ.PFO",	# Pinyon Flats Observatory, CA, USA
        "BK.YBH",	# Yreka Blue Horn Mine, Yreka, CA, USA
        "CN.DLBC",	# Dease Lake, BC, CA
        "CN.FRB",	# Iqaluit, NU, CA
        "CN.INK",	# Inuvik, NT, CA
        "CN.RES",	# Resolute, NU, CA
        "CN.SADO",	# Sadowa, ON, CA
        "CZ.VRAC",	# Vranov
        "DT.PMG",	# Port Moresby, Papua New Guinea
        "DW.LEM",	# Lembang, Indonesia
        "G.ATD",	# Arta Cave - Arta, Republic of Djibouti
        "GE.EIL",	# GII/GEOFON Station Eilat, Israel
        "GE.SFJD",	# IRIS/GEOFON Station Sondre Stromfjord, Greenland
        "GE.SNAA",	# AWI/GEOFON Station Sanae, Antarctica
        "GS.NEW",	# Seattle urban array
        "ID.NNA",	# Nana, Peru
        "ID.RPN",	# Easter Island (Rapa Nui), Chile
        "ID.SJG",	# San Juan, Puerto Rico
        "ID.SUR",	# Sutherland, Republic of South Africa
        "II.AAK",	# Ala Archa, Kyrgyzstan
        "II.BORG",	# Borgarfjordur, Asbjarnarstadir, Iceland
        "II.JTS",	# Las Juntas de Abangares, Costa Rica
        "II.KDAK",	# Kodiak Island, Alaska, USA
        "II.KURK",	# Kurchatov, Kazakhstan
        "II.MBAR",	# Mbarara, Uganda
        "II.OBN",	# Obninsk, Russia
        "II.PALK",	# Pallekele, Sri Lanka
        "IM.TKL",	# Tuckaleechee Caverns, TN, USA
        "IU.ANMO",	# Albuquerque, New Mexico, USA
        "IU.GNI",	# Garni, Armenia
        "IU.GUMO",	# Guam, Mariana Islands
        "IU.HNR",	# Honiara, Solomon Islands
        "IU.MSKU",	# Masuku, Gabon
        "IU.PMSA",	# Palmer Station, Antarctica
        "IU.QSPA",	# South Pole Remote Earth Science Observatory (Quiet Zone)
        "IU.RAO",	# Raoul, Kermadec Islands
        "IU.RCBR",	# Riachuelo, Brazil
        "IU.TEIG",	# Tepich, Yucatan, Mexico
        "IU.TSUM",	# Tsumeb, Namibia
        "JP.JKA",	# Kamikawa Asahi
        "JP.JNU",	# Oita Nakatsue
        "JP.JOW",	# Okinawa Kunigami
        "KZ.AKTO",	# Aktyubinsk, Kazakhstan
        "KZ.BVAR",	# Borovoye, Kazakstan
        "MN.IDI",	# Anogia, Greece
        "MN.MDT",	# Midelt, Morocco
        "MN.VAE",	# Valguarnera, Italy
        "NC.HES",	# Elkhorn Slough
        "NP.ELK",	# Elk Grove
        "NZ.RPZ",	# Rata Peaks
        "NZ.URZ",	# Urewera
        "PS.TGY",	# Tagaytay, Phillipines
        "SY.USHA",	# USHA synthetic
        "VE.PCRV",	# Pto. La Cruz - Anzoategui - CTBTO borehole
        "XY.BBTS",	# Black Bottom Tack and Stables - Summerville, SC
        "XA.I01H",	# in San Juan
        "XA.I02H",	# in San Juan
        "XA.I03H",	# in San Juan
        "XA.I04H",	# in San Juan
        "XA.I06H",	# in San Juan
        "XA.I09H",	# in San Juan
        "XA.I11H",	# in San Juan
        "XA.I15H",	# in San Juan
        "XA.I19H",	# in San Juan
        "XA.O01H",	# in San Juan
        "XA.O02H",	# in San Juan
        "XA.O03H",	# in San Juan
        "XA.O09H",	# in San Juan
        "XA.O10H",	# in San Juan
        "YC.JUAN",	# Pie de Palo (San Juan, ARG.)
        "ZL.BOZA",	# Barboza
        "AS.CTAO",	# Charters Towers, Australia
        "AU.CTA",	# Charters Towers, Queensland
        "HG.CTA",	# Charters Towers, Australia
        "IU.CTAO",	# Charters Towers, Australia
        "IU.NWAO",	# Narrogin, Australia
        "SR.NWAO",	# Narrogin, Australia
        "Z6.BRDL",	# Barydhala, Sitakund (BAEC)
        "GT.LBTB",	# Lobatse, Botswana, Africa
        "IU.PTGA",	# Pitinga, Brazil
        "CN.BBB",	# Bella Bella, BC, CA
        "GE.LVC",	# IRIS/GEOFON Station Calama, Chile
        "IU.LVC",	# Limon Verde, Chile
        "XM.CB04",	# CB04
        "XM.CB05",	# CB05
        "XM.CB17",	# CB17
        "CD.BJI",	# Baijiatuan, Beijing, China
        "IC.BJT",	# Baijiatuan, Beijing, China
        "XG.IGCSB",	# Inst. of Geophysics, China Seismological Bureau
        "CD.KMI",	# Kunming, Yunnan Province, China
        "IC.KMI",	# Kunming, Yunnan Province, China
        "CD.SSE",	# Sheshan, Shanghai, China
        "IC.SSE",	# Shanghai, China
        "IC.XAN",	# Xi'an, China
        "IU.RAR",	# Rarotonga, Cook Islands
        "CZ.KRUC",	# Moravsky Krumlov
        "GE.SFJ",	# IRIS/GEOFON Station Sondre Stromfjord, Greenland
        "IU.SFJ",	# Sondre Stromfjord, Greenland
        "IU.SFJD",	# Sondre Stromfjord, Greenland
        "G.AGD",	# Arta Grotte - Arta, Republic of Djibouti
        "MN.KEG",	# Kottamya, Egypt
        "II.MSVF",	# Monasavu, Fiji
        "G.DZM",	# Dzumac - New Caledonia, France
        "G.NOC",	# Noumea - New Caledonia, France
        "G.NOUC",	# Port Laguerre - New Caledonia, France
        "ZI.PAPA",	# PAPA
        "G.KOG",	# Kourou - French Guiana, France
        "G.MPG",	# Montagne des Peres - French Guiana, France
        "GE.FODE",	# GEOFON Station Fodele, Crete, Greece
        "XD.HOT01",	# Reykir, Iceland
        "PS.JAY",	# Jayapura, Indonesia
        "PS.PSI",	# Parapat, Indonesia
        "II.KAPI",	# Kappang, Sulawesi, Indonesia
        "HG.EIL",	# Eilat, Israel
        "GE.MRNI",	# GII Station Mount Meron, Israel
        "JP.JHJ2",	# Hachijojima Island
        "JP.ASAJ",	# Kamikawa Asahi
        "JP.CBIJ",	# Chichijima Island
        "JP.JCJ",	# Chichijima Island
        "PS.OGS",	# Chichijima, Bonin Islands, Japan
        "II.BORK",	# Burabay, Kazakhstan
        "II.BRVK",	# Borovoye, Kazakhstan
        "KZ.BORK",	# Burabay, Kazakhstan
        "KZ.BRV",	# Borovoye, Kazakstan
        "KZ.BRVK",	# Borovoye, Kazakstan
        "KZ.BVA1",	# Borovoye array,BVA1, Kazakhstan
        "KZ.BVA2",	# Borovoye array,BVA2, Kazakhstan
        "KZ.BVA3",	# Borovoye array,BVA3, Kazakhstan
        "KZ.BVA4",	# Borovoye array,BVA4, Kazakhstan
        "KZ.BVA5",	# Borovoye array,BVA5, Kazakhstan
        "KZ.BVA6",	# Borovoye array,BVA6, Kazakhstan
        "KZ.BVA7",	# Borovoye array,BVA7, Kazakhstan
        "KZ.BVA8",	# Borovoye array,BVA8, Kazakhstan
        "KZ.BVA9",	# Borovoye array,BVA9, Kazakhstan
        "KZ.KUR",	# Kurchatov, Kazakstan
        "KZ.KUR01",	# Cross-array, site1, Kurchatov, Kazakstan
        "KZ.KUR02",	# Cross-array, site2, Kurchatov, Kazakstan
        "KZ.KUR03",	# Cross-array, site3, Kurchatov, Kazakstan
        "KZ.KUR04",	# Cross-array, site4, Kurchatov, Kazakstan
        "KZ.KUR05",	# Cross-array, site5, Kurchatov, Kazakstan
        "KZ.KUR06",	# Cross-array, site6, Kurchatov, Kazakstan
        "KZ.KUR07",	# Cross-array, site7, Kurchatov, Kazakstan
        "KZ.KUR08",	# Cross-array, site8, Kurchatov, Kazakstan
        "KZ.KUR09",	# Cross-array, site9, Kurchatov, Kazakstan
        "KZ.KUR10",	# Cross-array, site10, Kurchatov, Kazakstan
        "KZ.KUR11",	# Cross-array, site11, Kurchatov, Kazakstan
        "KZ.KUR12",	# Cross-array, site12, Kurchatov, Kazakstan
        "KZ.KUR13",	# Cross-array, site13, Kurchatov, Kazakstan
        "KZ.KUR14",	# Cross-array, site14, Kurchatov, Kazakstan
        "KZ.KUR15",	# Cross-array, site15, Kurchatov, Kazakstan
        "KZ.KUR16",	# Cross-array, site16, Kurchatov, Kazakstan
        "KZ.KUR17",	# Cross-array, site17, Kurchatov, Kazakstan
        "KZ.KUR18",	# Cross-array, site18, Kurchatov, Kazakstan
        "KZ.KUR19",	# Cross-array, site19, Kurchatov, Kazakstan
        "KZ.KUR20",	# Cross-array, site20, Kurchatov, Kazakstan
        "KZ.KUR21",	# Kurchatov-Cross, Reference point, Kazakstan
        "KZ.KURK",	# Kurchatov, Kazakhstan
        "KZ.AKT",	# Aktyubinsk, Kazakstan
        "KZ.AKTK",	# Aktyubinsk, Kazakstan
        "KN.AAK",	# -
        "KR.FRU",	# Frunze, Kyrgyzstan
        "KR.FRU1",	# Frunze, Kyrgyzstan
        "KR.FRU2",	# Bishkek, Kyrgyzstan
        ]


# Stations not found on IRIS
not_found = {
        "PS27": [60.8, 10.8],  # Hamar, Norway
        "PS42": [35.7, 9.3],   # Kesra, Tunisia
        "PS16": [26.0, 33.5],  # Luxor, Egypt
        "PS38": [23.4, 44.5],  # Haleban, Saudi Arabia
        "PS35": [59.6, 112.6], # Peleduy, Russia
        # AS: DW.AFI, G.SEY, CD.BJI

    }

# Gather into single 
sta_dict = {"watc": watc_seismic_stations,
            "primary": ims_primary_seismic,
            "auxiliary": ims_auxiliary_seismic
            }

inv_write = Inventory()
for name, stations in sta_dict.items():
    inv = Inventory()
    for code in stations:
        net, sta = code.split(".")
        try:        
            inv_single = c.get_stations(network=net, station=sta, 
                                        # channel="BH?,HH?,SH?", 
                                        channel="*", level="station"
                                        )
        except Exception as e:
            print(f"Error with {code}")
            continue
        inv.extend(inv_single)

    inv.plot(outfile=f"figures/stations_{name}.png")
    inv.write(f"data/inv_{name}.xml", format="stationxml")
    inv_write.extend(inv)

write_stations_file(inv_write, fid=f"data/STATIONS", order_by="network")
    
    


