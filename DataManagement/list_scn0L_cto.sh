data='/rdata/airglow/gps/scintmonL/cto'
if [ -e "$data/foldersizes.txt" ]
then
	rm $data/foldersizes.txt
fi

for stuff in $( find $data/20*/raw_data/* )
do
	justyear=$( echo "$stuff"  | cut -c34-37 )
	justday=$( echo "$stuff" | cut -c48-59 )
	imsize=$( stat -c '%s' $stuff )
	echo $justyear $justday $imsize >> $data/foldersizes.txt

done
