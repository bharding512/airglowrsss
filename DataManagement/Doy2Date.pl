#!/usr/bin/perl -w

# Script to Convert Year and Doy to Date
# Returns Month and Day
# Daniel Fisher
# 20 December 2011
sub doy2date {

use strict;

my $year = $_[0];    #year
my $doy = $_[1];     #doy
my (@eom, $month);

if($year%4 == 0) {
   if($year%100 == 0) {
      if($year%400 == 0) {
         @eom = qw( 31 29 31 30 31 30 31 31 30 31 30 31 );
         print "Leap Century\n";
      } else {
         @eom = qw( 31 28 31 30 31 30 31 31 30 31 30 31 );
         print "365 Days in Year\n";
      }
   } else {
      @eom = qw( 31 29 31 30 31 30 31 31 30 31 30 31 );
      print "Leap Year\n";
   }
} else {
   @eom = qw( 31 28 31 30 31 30 31 31 30 31 30 31 );
   print "365 Days in Year\n";
}

for($month=0; $doy>0; $month++) {
   $doy = $doy - $eom[$month];
   if($month == 12) {   #Catch Errors when >366 days are used - out will be wrong
      $doy = -30;
      $month = 0;
   }
}
$doy = $doy + $eom[$month-1];

print("$month-$doy\n");
return($month,$doy);
}
1;
