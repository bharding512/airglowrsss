#!/usr/bin/perl -w

# Script to rotate the temperature logfile
# Needs 3 LETTER SITE when called and number of days back.
# Written by Jonathan J. Makela on 4 Aug 2011

#use Date::Manip;
## Initialize the Timezone
#Date_Init("TZ=BRT");

# Send in Site Location
my $site = @ARGV[0];
my $dayback = @ARGV[1];

# Days Back, get today's date
my @date = localtime(time - $dayback*86400); # Shifts 1 day back (Windows reads UTC -> 1day + 4 hrs offset)
$year = $date[5] + 1900;	# 4 digit year
$month = $date[4]+1;		# Month
$day = $date[3];		# Day of month


if($site eq 'CAR') {
   $Thermo_dir = "/cygdrive/c/Docume~1/MiniME/Desktop/ThermoHID/";
   $log_file = sprintf("TempHIDData_AEROLUME_%04d_%02d.txt",$year,$month);
   $sending_dir = "/cygdrive/c/Sending/";
} elsif($site eq 'CAJ') {
   $Thermo_dir = "/cygdrive/c/Docume~1/MiniME/Desktop/ThermoHID/";
   $log_file = sprintf("TempHIDData_WINXP-014BA5828_%04d_%02d.txt",$year,$month);
   $sending_dir = "/cygdrive/c/Sending/";
}

my $header_file = "ThermoHID_head.txt";

# Create filename
my $fname = sprintf("TempL_%s_%04d%02d%02d.txt", $site, $year, $month, $day);

# Copy the file over
my $cmd = "cp $Thermo_dir$log_file $sending_dir/$fname";
system($cmd);

## Copy the header file
#$cmd = "cp $Thermo_dir$header_file $Thermo_dir$log_file";
#system($cmd);
