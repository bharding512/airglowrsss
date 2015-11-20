#!/bin/sh

# Turn the filterwheel and camera off
/home/airglow/scripts/lpcperl.pl 192.168.0.100 admin:ionosphere 1off
/home/airglow/scripts/lpcperl.pl 192.168.0.100 admin:ionosphere 5off
