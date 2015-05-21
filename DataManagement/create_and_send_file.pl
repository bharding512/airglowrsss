#!/usr/bin/perl
use warnings;
use strict;

# version 1.3
# Tim Duly (duly2@illinois.edu)
# Created Aug 2011
# 
# (1.3) updated to send systems disk usage.

my $site_name = $ARGV[0];
# =-=-=-=-==-=-=-=-==-=-=-=-==-=-=-=-=


my $scp_rcv;
# we have to treat trinidad special:
if ($site_name =~ /SGT_PICASSO|SGT_SCINDA/) {
  $scp_rcv = "airglow\@192.168.0.100:~/internet_connectivity_rcv/"; # for Trinidad machines, send to the Windows Internet NetBook
} elsif ($site_name =~ /CAJ_SCINT/){
  $scp_rcv = "picasso\@192.168.0.177:~/internet_connectivity_rcv/"; # for CAJ_SCINT
} else {
  $scp_rcv = "data\@airglow.ece.illinois.edu:/usr/local/share/datamanagement/internet_connectivity_rcv/";
# everyone else send to airglow
}


my %df_path;

$df_path{TIM} = "/dev/sda5";

$df_path{ANN} = "C:";
$df_path{EKU} = "C:";
$df_path{UAO} = "C:";
$df_path{BON} = "";
$df_path{CAJ} = "C:";
$df_path{CAJ_PICASSO} = "/dev/sda1";
$df_path{CAJ_SCINDA} = "/dev/sda1";
$df_path{CAJ_SCINT} = "/dev/hda6";
$df_path{CAR} = "C:";
$df_path{COL_SCINT} = "/dev/hda6";
$df_path{CTO} = "/dev/hda2";
$df_path{CTO_SCINT_L} = "/dev/hda6";
$df_path{CTO_SCINT_M} = "/dev/hda6";
$df_path{HKA_CASI} = "C:/cygwin";
$df_path{HKA_CNFI} = "D:";
$df_path{MRH} = "D:";
$df_path{NSO} = "/dev/sda5";
$df_path{PAR} = "";
$df_path{SCO} = "C:/cygwin/bin";
$df_path{SGT_PICASSO} = "/dev/sda7";
$df_path{SGT_SCINDA} = "/dev/sda1";

my %run_path;

$run_path{TIM} = "/home/duly/Documents/SourceCode/DataManagement";
$run_path{ANN} = "/cygdrive/c/Scripts/DataManagement";
$run_path{EKU} = "/cygdrive/c/Scripts/DataManagement";
$run_path{UAO} = "/cygdrive/c/Scripts/DataManagement";
$run_path{BON} = "/home/airglow/scripts/DataManagement";
$run_path{CAJ} = "/cygdrive/c/Scripts/DataManagement";
$run_path{CAJ_PICASSO} = "/home/picasso/scripts/DataManagement";
$run_path{CAJ_SCINDA} = "/home/gps/scripts/DataManagement";
$run_path{CAJ_SCINT} = "/home/scintmon/scripts/DataManagement";
$run_path{CAR} = "/cygdrive/c/Scripts/DataManagement";
$run_path{COL_SCINT} = "/home/gps/scripts/DataManagement";
$run_path{CTO} = "/home/airglow/scripts/DataManagement";
$run_path{CTO_SCINT_L} = "/home/gps/scripts/DataManagement";
$run_path{CTO_SCINT_M} = "/home/gps/scripts/DataManagement";
$run_path{HKA_CASI} = "/home/AGUser/scripts/DataManagement";
$run_path{HKA_CNFI} = "/home/AGUser/scripts/DataManagement";
$run_path{MRH} = "/home/merihill/scripts/DataManagement";
$run_path{NSO} = "/home/airglow/scripts/DataManagement";
$run_path{PAR} = "/cygdrive/c/Scripts/DataManagement";
$run_path{SCO} = "/home/AGUser/scripts/DataManagement";
$run_path{SGT_PICASSO} = "/home/airglow/scripts/DataManagement";
$run_path{SGT_SCINDA} = "/home/gps/scripts/DataManagement";

my $result = system("touch $run_path{$site_name}/$site_name");# or die "couldnt write file";
#my $result = system("df -h > $site_name");
chomp (my $time = `date`);

if ($result eq 0) {
  print "[$time]: created internet connectivity test file for $site_name.\n"
}


my @dfs = `df`;
my $df_line = "na";
foreach my $df (@dfs) {
  chomp $df;
  if ($df =~ /($df_path{$site_name}.*)/) {
    $df_line = $1;
  }
}

print "captured:\n$df_line\n";

my @disk_cap = split /\s+/, $df_line;
my $disk_capacity = @disk_cap[3];

print "captured:\n$disk_capacity\n";

open FILE,'>',"$run_path{$site_name}/$site_name";
printf FILE $disk_capacity;
close FILE;


my $scp_result = 999; 
my $count = 1;

while ( ($count <= 5) and ($scp_result ne 0)){
  print "Try $count of scp-ing...\n";
  $scp_result = system("scp $run_path{$site_name}/$site_name $scp_rcv");
  $count += 1;
}

print "$scp_result\n";

if ($scp_result eq 0) {
  print "[$time]: file \"$site_name\" successfully scp'd.\n";
  system("rm $run_path{$site_name}/$site_name");
} else {
  print "there was a problem scp'ing at $site_name\n";
}



print "\n";

