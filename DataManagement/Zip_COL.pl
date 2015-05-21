#!/usr/bin/perl -w

# Script to Zip, Split and Move recent data to Sending Folder
# Daniel Fisher
# 30 August 2011

#use strict;
chdir;
chomp(my $dir_doy = `pwd`);
require "$dir_doy/scripts/DataManagement/Doy2Date.pl";
my ($doy,$year,$yr,$month,$day);
my $dayback = $ARGV[0];  #Days back or year
my $doyd = $ARGV[1];    #Doy 
my $site = 'bog';
my $dir_send = '/home/gps/Sending';

if($dayback < 2000) {
   # Days Back, get yesterdays's date
   my @date = localtime(time-108000 - $dayback*86400); # Shifts 1 day back (Windows reads UTC -> 1day + 4 hrs offset)
   $doy = $date[7]+1;		# Day of year
   $year = $date[5] + 1900;	# 4 digit year
   $yr = $year - 2000;		# 2 digit year
   $month = $date[4]+1;		# Month
   $day = $date[3];		# Day of month
} else {
   # Actual Given Date
   $year = $dayback;
   $doy = $doyd;  
   ($month,$day) = doy2date($year,$doy);
   $yr = $year - 2000;		# 2 digit year
}


# Send SCN O ===============================================
my $instrument = 'scn0O';
my $local_stub = sprintf("/home/gps/cascade");
my $name = sprintf("%02d%02d%02d*",$yr,$month,$day);
my $filename = sprintf("%s\_%s\_%04d%02d%02d.tar.gz",$instrument,$site,$year,$month,$day);
chdir("$local_stub");
# Gun-Tar Folder
`tar czvf $filename $name`;
# Split File & Move it
my $filesize = -s "$local_stub/$filename";
print("$filesize\n");
if($filesize ge 1000) {
   system("split -a 6 -b 5242880 -d $filename $dir_send/$filename");
   print "$filename Ready to Send\n";
}
# Remove Intermediary
`rm -f $local_stub/$filename`;

print "Zip Complete\n";
