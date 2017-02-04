data='/rdata/airglow/fpi/minime06/par'
if [ -e "$data/foldersizes.txt" ]
then
	rm $data/foldersizes.txt
fi

for year in $( find $data/20* -maxdepth 0 -type d )
do
	justyear=$( echo "$year"  | cut -c33-36)
	
	# Skip this year because we had filenaming issues.
	if [$justyear == 2012]
	then
		continue
	fi

	for day in $( find $year/* -maxdepth 0 -type d )
	do
		justday=$( echo "$day" | cut -c42-46 )
		imsize=$( find $day'/' -name *_*.img -print0 | du --files0-from=- -hc -b | tail -n1 | cut -f1 )
		filesize=$( du -b $day | cut -f1 )
		echo $justyear $justday $imsize >> $data/foldersizes.txt
	done

done
