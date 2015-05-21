#!/bin/bash

# This is specifically for HKA, (the CASES receiver), if you ever have another system that doesn't have perl,
# I would need to think of a clever way to use hashes in bash... (can you?)
# Modificaiton for Sortinghat v3.0 on 2/28/14 - djf


# Get Filename
path="/mnt/data/scripts/DataManagement/"
name="cas01_hka_"
time="$(date +'%Y%m%d')"
end=".txt"
fn=$path$name$time$end

# Scp Location
scp_rcv="tx@remote2.ece.illinois.edu:/rdata/airglow/rx/tracking/."

# Get most recent file
file=$(find /mnt/data -type f | xargs ls -ltr | tail -n 1 | awk '{ print $9 }')

# Create Checkfile
cat > check.txt << EOF
cas01_hka_createandsendfilesh
0
0
EOF
stat -c %y "$file" >> check.txt
df -B 1073741824 /mnt/data | awk '{ print $4 }' | tail -n 1 >> check.txt
echo "Latest File: $file" >> check.txt
mv check.txt $fn

# Send it over and remove
scp $fn $scp_rcv
rm $fn

