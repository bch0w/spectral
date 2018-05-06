#!/bin/bash
# Playback set up - As of Aug2013-21018
######################

### Paths      
#Path2XML="/home/salichon/seiscomp3VM/Sc3hik/PlayBck/"
Path2Res="/home/salichon/seiscomp3VM/Sc3hik/PlayBck/RESULTS"
### Output directory 
OutputDir="$Path2Res/OutputXML_and_Bulls"

###  streams and  data services 
SDS=sdsarchive
SDSDirectory="/work/seiscomp"
Source=$SDS://$SDSDirectory

ScartFile="Hik_Stat.txt" # To be edited  - see other batch 

### Seedlink combined services slink... with SDS archive -  to de set up if ever required
#service=combined://;
#service=sdsarchive://;
### 

### DB read options 
DBFLAG="-d mysql://sysopGNS:sysop@localhost/seiscomp3"
## DB write options Storage of results
WDBFLAG="-d mysql://sysopGNS:sysop@localhost/seiscomp3"
WSTORAGE=$WDBFLAG

## no DB (full) read option ?  
#XMLFILES="--inventory-db $Path2XML/...xml  --config-db $Path2XML/.....xml " 

## Logging
VERBOS=" --verbosity=4 "

### Options set-up 

## Using Database
FLAGS=" --console=0 --verbosity=0 $DBFLAG"
FLAGS_DeBg="  $VERBOS  --trace $DBFLAG  --debug "

## Don t read database  
#FLAGS_DeBg=" $VERBOS --debug $XMLFILES "

####### SC3 seattle command
Sc3ex="seiscomp exec "

### Channel used (scart usage)
#ChAnnels="(HH|EH)(Z|N|E|1|2)"
ChAnnels="(H|B)(N)(Z)" # only Z are picked there !
#Netwk="NZ"

### Time window set up 
	#echo $tnew1 $tnew2 
    	#echo "   2008-02-14 00:00:00~2008-03-01 00:00:00"
    	#echo "   		t1		t2	"
    	#echo " read  start time t1"

	## Default time windows
	t1="2017-08-11 09:50:00"
	#t1="2017-08-11 00:00:00"
#	read  tnew # exported
	[ -n "$tnew1" ] && t1=$tnew1

	#echo " read end time t2"
	t2="2017-08-11 10:35:00"
	#t2="2017-08-11 00:35:00"

        #echo "Default end  date:  $t2"  >> ResplayLog.txt

        #echo "enter another one if required:"
    	#read  tnew2 # exported
	[ -n "$tnew2" ] && t2=$tnew2

	Time_Stamp1=`date -d "$t1" +%Y%m%d-%H%M%S  `
	Time_Stamp2=`date -d "$t2" +%Y%m%d-%H%M%S  `
###
    #read  name_xml
	name_xml="$Time_Stamp1"_"$Time_Stamp2".xml
    #echo ""
    #echo "Usage: $0  [time window: $t1~$t2] [output-xml: $name_xml]"
    #echo ""
	Event_file="${Time_Stamp1}_${Time_Stamp2}.EvtListe.txt"
    #echo ""
echo ""
echo ""

######################
	echo "Starting scart|scautopick ..."
echo "$Sc3ex scart -dsE --list $ScartFile   $SDSDirectory  | $Sc3ex scautopick --ep --playback -I - $FLAGS_DeBg --logging.file=false > Picks_${Time_Stamp1}_${Time_Stamp2}.xml"
#read
#$Sc3ex scart -dsE --list $ScartFile   $SDSDirectory  | $Sc3ex scautopick --ep --playback -I - $FLAGS_DeBg --logging.file=false > Picks_${Time_Stamp1}_${Time_Stamp2}.xml
$Sc3ex scart -dsE --list $ScartFile   $SDSDirectory  | $Sc3ex scautopick --ep --playback -I - $FLAGS --debug  --logging.file=false > Picks_${Time_Stamp1}_${Time_Stamp2}.xml

$Sc3ex scautoloc --ep Picks_${Time_Stamp1}_${Time_Stamp2}.xml --playback $DBFLAG --debug > Origins_${Time_Stamp1}_${Time_Stamp2}.xml
$Sc3ex scamp   -I $Source --ep Origins_${Time_Stamp1}_${Time_Stamp2}.xml  --playback  $DBFLAG --debug > amps_${Time_Stamp1}_${Time_Stamp2}.xml
$Sc3ex scmag --ep amps_${Time_Stamp1}_${Time_Stamp2}.xml  --playback  $DBFLAG --debug > mags_${Time_Stamp1}_${Time_Stamp2}.xml
$Sc3ex scevent --ep mags_${Time_Stamp1}_${Time_Stamp2}.xml  --playback  $DBFLAG --debug >  events_${Time_Stamp1}_${Time_Stamp2}.xml

#$Sc3ex scevent $XMLFILES $DBFLAG --debug --db-disable --start-stop-msg=1 --auto-shutdown=1 --shutdown-master-module=scmag &

######################

# NB generate mseed sample instead
##	echo "$Sc3ex scart  -vvv --test -n $Netwk -dsE -c"$ChAnnels" -t "$t1"~"$t2" $SDSDirectory "
##	$Sc3ex scart  -vvv  -n $Netwk -dsE -c"$ChAnnels" -t "$t1"~"$t2" $SDSDirectory  > test.mseed 
