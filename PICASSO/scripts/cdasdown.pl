#!/usr/bin/perl -w 
# cdasdown script connects to CDAS via port 30000 on localhost and sends 
#  shutdown command (0:)
# Written by Ethan S. Miller (esmiller@illinois.edu) 13 March 2009
#
# Mostly stolen from:
#  http://www.perl.com/doc/FMTEYEWTK/IPC/inet.html
#

require 5.002; 
use strict; 
use Socket; 

print "Shutting down CDAS.\n";

my ($remote,$port, $iaddr, $paddr, $proto, $line); 
$remote = shift || 'localhost'; 
$port = shift || 30000; 

$iaddr = inet_aton($remote) || die "no host: $remote"; 
$paddr = sockaddr_in($port, $iaddr); 
$proto = getprotobyname('tcp'); 
socket(SOCK, PF_INET, SOCK_STREAM, $proto) || die "socket: $!"; 
connect(SOCK, $paddr) || die "connect: $!"; 

print SOCK "0:";

close (SOCK) || die "close: $!"; 
exit; 
