#!/bin/bash
# sync common folder and data folder up and down on convolution
# UP means from convolution up to computer
# DOWN means down from computer to convolution
# to run: sh rsync_folders.sh GNS/VIC UP/DOWN DRY/GO 
VICORGNS=$1
UPORDOWN=$2
DRYRUN=$3
VIC='VIC'
GNS='GNS'

# set variables 
if [ $VICORGNS = 'GNS' ]
	then
		COMMONDIR='/seis/prj/fwi/bchow/spectral/common'
		PAPERDIR='/home/bchow/Documents/papers'
		CONVOLUTION='/run/media/bchow/convolution'
elif [ $VICORGNS = $VIC ]
	then
		COMMONDIR='/Users/chowbr/Documents/subduction/spectral/common'
		PAPERDIR='/Users/chowbr/Documents/papers'
		CONVOLUTION='/Volumes/convolution'
else
	echo 'INVALID FIRST COMMANDLINE ARGUMENT'
fi

# check if dry run, default to dry run for safety
if [ $DRYRUN = 'GO' ]
	then
		echo [REALRUN] rsyncing $UPORDOWN, $VICORGNS
		DRY=''
else
	echo [DRYRUN] rsyncing $UPORDOWN, $VICORGNS  
	DRY='n'
fi

# run rsync
SLASH='/'
COMMON='common'
PAPER='papers'
if [ $UPORDOWN = 'UP' ]
	then
	    echo $ rsync -gloptrucv$DRY $CONVOLUTION$SLASH$PAPER $PAPERDIR$SLASH
	    rsync -gloptrucv$DRY $CONVOLUTION$SLASH$PAPER $PAPERDIR$SLASH
		echo $ rsync -gloptrucv$DRY $CONVOLUTION$SLASH$COMMON $COMMONDIR$SLASH
		rsync -gloptrucv$DRY $CONVOLUTION$SLASH$COMMON $COMMONDIR$SLASH
elif [ $UPORDOWN = 'DOWN' ]
	then
		echo $ rsync -gloptrucv$DRY $PAPERDIR $CONVOLUTION$SLASH
		rsync -gloptrucv$DRY $PAPERDIR $CONVOLUTION$SLASH
		echo $ rsync -gloptrucv$DRY $COMMONDIR $CONVOLUTION$SLASH
		rsync -gloptrucv$DRY $COMMONDIR $CONVOLUTION$SLASH
else
	echo 'INVALID SECOND COMMANDLINE ARGUMENT'
fi
