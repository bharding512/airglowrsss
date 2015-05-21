data='/rdata/airglow/gps/scintmonU/caj'
if [ -e "$data/foldersizes.txt" ]
then
	rm $data/foldersizes.txt
fi

# Look in folder
for stuff in $( find $data/20*/*[.fsl,.nav,.obs] -maxdepth 0 -type f)
do
	justfile=$( echo "$stuff"  | cut -c39-50 )
	imsize=$( stat -c '%s' $stuff )
    echo $justfile  $imsize >> $data/foldersizes.txt

done

# Look in raw_data
for stuff in $( find $data/20*/raw_data/*[.fsl,.nav,.obs] -maxdepth 0 -type f)
do
	#justyear=$( echo "$stuff"  | cut -c34-37 )
	justfile=$( echo "$stuff"  | cut -c48-59 )
	imsize=$( stat -c '%s' $stuff )
	#echo $justyear $justdate  $imsize >> $data/foldersizes.txt
    echo $justfile  $imsize >> $data/foldersizes.txt

done
