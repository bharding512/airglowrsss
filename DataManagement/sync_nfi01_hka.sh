#!/bin/bash
data='/cygdrive/c/Data'
airglowdata='/rdata/airglow/imaging/cnfi01/hka/'
servername='airglow@remote2.ece.illinois.edu'
FILENAME=/cygdrive/c/Data/foldersizes.txt

echo "You will need to enter the scp password..."

scp $servername:$airglowdata/foldersizes.txt $data/.
if [ -e "$data/foldersizes.txt" ]
then
	while read LINE
	do
		listday=$( echo $LINE | cut -d ' ' -f2 )
		listyear=$( echo $LINE | cut -d ' ' -f1 )
		temp=$( echo $LINE | cut -d ' ' -f3 )
		let listfilesize=$temp
		let index=10#$listday+366*10#$listyear-366*2001
		sizearray[$index]=$listfilesize
	done <$FILENAME

	echo "text file parsed"
	for year in $( find $data/20* -maxdepth 0 -type d )
	do
		justyear=$( echo "$year"  | cut -c18-21)
		
		for day in $( find $year/* -maxdepth 0 -type d )
		do
			justmonth=$( echo "$day" | cut -c23-25 )
			justday=$( echo "$day" | cut -c26-27 )
			doy=$( date +%j -d "$justmonth $justday $justyear" )
			let doy=10#$doy+1
			doy=$( printf "%03d" $doy )
			imsize=$( find $day'/' -name *_*.tif -print0 | du --files0-from=- -hc -b | tail -n1 | cut -f1 )
			temp=$( du -b $day | cut -f1 )
			let filesize=$temp
			let index=10#$doy+366*10#$justyear-366*2001

			if [ -z ${sizearray[$index]} ]
			then
				echo $justyear $doy "doesn't exist"
				scp -r $day/* $servername:/$airglowdata/$justyear/$doy
			elif [ $filesize -gt ${sizearray[$index]} ]
			then
				echo $justyear $doy "needs to send" 
				echo $filesize ${sizearray[$index]}
				scp -r $day/* $servername:/$airglowdata/$justyear/$doy
			elif [ $filesize -eq ${sizearray[$index]} ]
			then
				rm -r $day
				echo $justyear $doy "is fine"
				echo $filesize ${sizearray[$index]} 
			fi
				

		done

	done

else
	echo 'shit'
fi
