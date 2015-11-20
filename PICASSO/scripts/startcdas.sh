#!/bin/sh

# turn off the power to the filterwheel and camera
/home/airglow/scripts/peripherals_off.sh
sleep 1

# turn on the power to the filterwheel and camera
/home/airglow/scripts/peripherals_on.sh
sleep 10

# run cdas
cd home/airglow; cdas
