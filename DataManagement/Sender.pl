#!/usr/bin/perl -w

# Script to Send all tar.gz to the Rx Airglow Folder, Needs 3 Letter Site Abbr.
# Daniel Fisher
# 29 August 2011

use File::Basename;

## Code to stop if alread running this script
#chomp(my $result = `pidof -x Sender.pl`);
#my @result = split(/ /,$result);
#my $processes = scalar(@result);
#print "$processes\n";
#
#if($processes gt 1) {
#   exit;
#}

# Get in Proper Sending Directory
my $site = $ARGV[0];
print("-$site-\n");
if($site eq 'CAR') {
   $local = sprintf("/cygdrive/f/Sending");
} elsif($site eq 'NSO') {
   $local = sprintf("/data/Sending");
} elsif($site eq 'CAJ_P') {
   $local = sprintf("/data/Sending");
} elsif($site eq 'CAJ_S') {
   $local = sprintf("/home/scintmon/Sending");
} elsif($site eq 'CAJ_D') {
   $local = sprintf("/home/gps/Sending");
} elsif($site eq 'CTO') {
   $local = sprintf("/home/airglow/Sending");
} elsif($site eq 'SGT_P') {
   $local = sprintf("/home/airglow/Sending");
} elsif($site eq 'SGT_T') {
   $local = sprintf("/home/gps/Sending");
} elsif($site eq 'COL') {
   $local = sprintf("/home/gps/Sending");
} elsif($site eq 'CTO_L') {
   $local = sprintf("/home/gps/Sending");
} elsif($site eq 'UAO_Cloud') {
   $local = sprintf("/cygdrive/d/Sending");
} else {
   $local = sprintf("/cygdrive/c/Sending");
} 

my $dest = sprintf('tx@remote2.ece.illinois.edu:/rdata/airglow/rx');
chdir($local);

# Send All MiniME files
@files = `ls TempL* Cloud* *.tar.gz*`;
foreach $file(@files) {
   chomp($file);
   system('chmod 777',"$file")
   my $result = system('scp', "$file", "$dest");
   if($result eq 0) {
      system("rm -f $file");
      print("$file Sent Successfully\n");
   }
}

print "Complete!\n";

