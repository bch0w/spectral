#!/bin/bash
# Catalog process management

#nbr of days to process
NDay=1 

#### Check if the process is already running 
if mkdir ~/.Catalogsproc.exclusivelock_Sc3Hik
then
	echo `tail -1 DaysDone.txt`
	firstdate=`tail -1 DaysDone.txt`

	### Define the dates to replay
	# Starting date
	# We are going backward ....
	#tINIT=`date -d "2013-06-01" +%Y-%m-%d ` 
	tINIT=$firstdate 
	
	echo "##"
	echo `more Readme.txt` 
	echo "##"
	echo "Replay"
	echo "starting date $tINIT "
		for (( i=1 ; i<=$NDay ; i++ ))
		do 
		nrday=$i
		tnew1=`date -d "$tINIT -$i days" +%Y-%m-%d` 
		tnew2=`date -d "$tnew1 +1 day +15 minutes" +'%Y-%m-%d %H:%M:%S' ` 
		#tnew2=`date -d "$tnew1  +15 minutes" +'%Y-%m-%d %H:%M:%S' ` 

		Trun=$tnew1 

		tnew1=`date -d "$tINIT -$i days" +'%Y-%m-%d %H:%M:%S'` 
		export tnew1 
		export tnew2 
		export nrday 
		echo $nrday > daynbr.tmp
		## if using scart + a list of station isntead of the whole network:
			Stat_Ref="Hik_Stat.txt_ref"
			Stat_file=" Hik_Stat.txt"
			rm Hik_Stat.txt
			st1=`date -d "$tINIT -$i days" +%Y-%m-%d` 
			st2=`date -d "$tnew1 +1 day " +'%Y-%m-%d'` 
			#sed 's/1999-99-99/"$st1"/' $Stat_Ref > $Stat_file
			sed s/1999-99-99/"$st1"/ $Stat_Ref > $Stat_file
			sed -i s/2999-99-99/"$st2"/ $Stat_file
			#sed -i 's/2999-99-99/"$st2"/' $Stat_file
		#	exit	
		#process to run
		time ./Hik_Stat_TimeWindoProcess.bash
		echo $Trun >> DaysDone.txt
		echo "Days to be treated:" $NDay 
		echo "Days  treated:" $nrday 
		echo  "" 
		done
	echo "Ending Loop  dates $tnew1 to  $tnew2  " 
	echo " ************  " 
	echo "  " 
	rmdir ~/.Catalogsproc.exclusivelock_Sc3Hik
	echo ""
	echo ""
	echo "#### FINISHED SET ####"
	echo "dates done $tnew1 to  $tINIT  " 
	echo ""
	echo ""
	echo "######"
else
	echo ""
	echo "Already running script .... LOCKED safety"
	echo "" 
	echo "days to be treated:" $NDay 
	echo "days  treated:" $nrday 
	echo "Exit script" 
	echo  `date`
	echo " " 
	exit
fi

echo "#### FINISHED LOOP ####"
echo "######"
echo "######"
