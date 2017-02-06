data='/rdata/airglow/imaging/picasso04/mor'
if [ -e "$data/foldersizes.txt" ]
then
	rm $data/foldersizes.txt
fi

for year in $( find $data/20* -maxdepth 0 -type d )
do
	justyear=$( echo "$year"  | cut -c38-42)

	for day in $( find $year/* -maxdepth 0 -type d )
	do
		justday=$( echo "$day" | cut -c43-45 )
		imsize=$( find $day'/' -name *_*.tif -print0 | du --files0-from=- -hc -b | tail -n1 | cut -f1 )
		filesize=$( du -b $day | cut -f1 )
		echo $justyear $justday $imsize >> $data/foldersizes.txt
	done

done
