#!/usr/bin/perl

# CDAS_monitor.pl - Monitor process queue to make sure that the
# CDAS program is running when it should be running.  This
# script should be placed in crontab to run every 15 minutes.
#
# */15  * * * *   /usr/local/bin/CDAS_monitor.pl
# 
# Written by Jonathan J. Makela (jmakela@illinois.edu) based upon
# a script used to control MythTV (written by David Brodbeck; 
# gull@gull.us)

# If any of the programs in this list are running, we do nothing.
# This should be the program name for the CDAS control program on
# this particular machine.

`PATH=/bin:/usr/bin:/usr/local/bin:/usr/X11R6/bin:/sbin:/usr/sbin:/usr/local/sbin`;
`export PATH`;

@expected_progs = ("cdas");

# Programs we need
$ps = "/bin/ps";
#$cdas = "cd /home/airglow; export LD_LIBRARY_PATH=/usr/local/lib; /usr/bin/xterm -e /home/airglow/run_cdas";
#$cdas = "cd /home/airglow; /usr/bin/xterm -e /home/airglow/run_cdas";
$cdas = "cd /home/airglow; /usr/bin/xterm -e /home/airglow/run_cdas";

# Check to make sure we are at a time when we expect the CDAS
# program to be running
if (cdastime()) {
	# CDAS should be running.  Verify this.
	$pslist = `$ps ax`;
	foreach $prog(@expected_progs) {
		if($pslist =~ $prog) {
			# The program is running, exit
			print "CDAS currently running.  Doing nothing.\n";
			exit;
		}
	}

	# We should only get here if the CDAS program is NOT
	# running.
	#
	# Run CDAS

	print "CDAS not running.  Starting program ($cdas).\n";
	`$cdas`;	
}

# Subroutine to figure out if CDAS should be running or not.  The
# variables ontime and offtime should be the localtime at which
# we need to make sure the system is turned on and off, respectively.
# These should match the times used in crontab to turn the respective
# power switches on/off if using web power switches.
sub cdastime {
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;

    $ontime = 18;
    $offtime = 8;  
 
    return 1 if (($hour >= $ontime || $hour <= $offtime));

    print "Not the correct time to be running CDAS.\n";
    return 0;
}
