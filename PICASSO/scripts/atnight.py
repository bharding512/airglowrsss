#!/usr/bin/env python

"""
Dependencies:
>= python2.7
PyEphem
psutil
sh
"""

import sys
import os
import logging
import time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple
from datetime import datetime, timedelta
from StringIO import StringIO

import ephem
import psutil
import sh

LOGGER = logging.getLogger('atnight')


"""
Default location for CDAS configuration file.
"""
CONFIG_FNAME = os.path.expanduser('~/CDAS.conf')

"""
Default log file.
"""
LOG_FNAME = os.path.expanduser('~/atnight.log')

"""
Name of PICASSO system process.
"""
CDAS_PROCESS_NAME = 'cdas'

"""
Default period (in hours) after sunset to turn on system.
"""
DELTA_SUNRISE = 1.

"""
Default period (in hours) prior to sunset to turn on system.
"""
DELTA_SUNSET = 1.

"""
Peripherals on script.
"""
PERIPHERALS_ON_SCRIPT = os.path.expanduser('~/scripts/peripherals_on.sh')

"""
Peripherals off script.
"""
PERIPHERALS_OFF_SCRIPT = os.path.expanduser('~/scripts/peripherals_off.sh')

"""
Location of cdas program.
"""
CDAS = '/usr/local/bin/cdas'

"""
Location of cdas kill script.
"""
CDAS_DOWN_SCRIPT = os.path.expanduser('~/scripts/cdasdown.pl')


"""
Container for ASI location information.
"""
class SiteInfo(namedtuple('SiteInfo', 'location id lat lon alt')):
    pass


def get_site_info(config_fname):
    """
    Parse *config_name* and return :class:`SiteInfo` containing
    information regarding the site location.
    """
    with open(config_fname) as fid:
        while True:
            try:
                line = fid.next()
            except StopIteration:
                break
            if line.startswith('#Site'):
                location = fid.next().strip()
                site_id = fid.next().strip()
                lat = float(fid.next())
                lon = float(fid.next())
                # convert longitude to range [-180, 180)
                if lon >= 180:
                    lon -= 360
                # convert [km] to [m]
                alt = float(fid.next()) * 1000
                return SiteInfo(location,
                                site_id,
                                lat,
                                lon,
                                alt)
    raise RuntimeError('could not find #Site block in configuration file '
                       + config_fname)


def is_system_on(name=CDAS_PROCESS_NAME):
    """
    Return `True` if the system is currently on. System status is
    determined by the presence of the active process named *name*.
    """
    for process in psutil.process_iter():
        try:
            if process.name() == name:
                return True
        except:
            # zombie processes cause problems --- skip them
            continue
    return False


def is_nighttime(utc_dt,
                 lat,
                 lon,
                 elev,
                 pressure=0,
                 horizon='-0:34',
                 delta_sunrise=timedelta(hours=1),
                 delta_sunset=timedelta(hours=1)):
    """
    Return `True` if :class:`datetime` *utc_dt* at the given latitude
    *lat* [+N deg], longitude *lon* [+E deg], and elevation *elev* [m]
    is between sunset (minus *delta_sunrise*) and sunrise (plus
    *delta_sunset*). The optional parameter *pressure* is the
    atmospheric pressure in [mBar] --- atmospheric refraction is
    ignored when set to 0. The optional parameter *horizon* manually
    sets the horizon to determine if the Sun is visible (use `None` to
    ignore this parameter). The defaults *pressure*=0 and
    *horizon*=-0:34 correspond to the USNO standard of -34 arcminutes
    below the normal horizon (see
    http://rhodesmill.org/pyephem/rise-set.html).
    """
    site = ephem.Observer()
    site.lat = str(lat)
    site.lon = str(lon)
    site.elevation = elev
    site.pressure = pressure
    if horizon is not None:
        site.horizon = horizon
    site.date = utc_dt.strftime('%Y-%m-%d %H:%M:%S')
    sun = ephem.Sun()
    sun.compute(site)
    if sun.alt <= 0:
        # approximately nighttime (Sun is below the horizon)
        return True
    else:
        # approximately daytime
        try:
            sunset = site.next_setting(sun).datetime()
        except ephem.AlwaysUpError:
            return False
        try:
            sunrise = site.previous_rising(sun).datetime()
        except ephem.NeverUpError:
            return True
        return utc_dt >= sunset - delta_sunset or utc_dt <= sunrise + delta_sunrise


def do_system_on(peripherals_on=PERIPHERALS_ON_SCRIPT,
                 peripherals_off=PERIPHERALS_OFF_SCRIPT,
                 cdas=CDAS,
                 config_fname=CONFIG_FNAME,
                 no_wps=False):
    """
    Turn the camera system on (mimics startcdas.sh script).
    """
    LOGGER.info('Turning system ON.')
    if not no_wps:
        # turn off peripherals via the web power switch
        LOGGER.info('Turning off peripherals')
        peripherals_off_cmd = sh.Command(peripherals_off)
        log_string = StringIO()
        peripherals_off_cmd(_out=log_string,
                            _err_to_out=True)
        LOGGER.debug(log_string.getvalue())
        time.sleep(1)
        # turn on peripherals via the web power switch
        LOGGER.info('Turning on peripherals')
        peripherals_on_cmd = sh.Command(peripherals_on)
        peripherals_on_cmd(_out=log_string,
                           _err_to_out=True)
        LOGGER.debug(log_string.getvalue())
        time.sleep(10)
    # run cdas program in the background
    LOGGER.info('Running {} in background'.format(cdas))
    sh.nohup(cdas,
             _bg=True,
             _cwd=os.path.dirname(config_fname),
             _in=os.devnull,
             _out=os.devnull,
             _err=os.devnull)


def do_system_off(peripherals_off=PERIPHERALS_OFF_SCRIPT,
                  cdas_down=CDAS_DOWN_SCRIPT,
                  no_wps=False):
    """
    Turn the camera system off.
    """
    LOGGER.info('Turning system OFF.')
    # cleanly kill the cdas program
    LOGGER.info('Cleanly stopping cdas with {}'.format(cdas_down))
    log_string = StringIO()
    cdas_down_cmd = sh.Command(cdas_down)
    cdas_down_cmd(_out=log_string,
                  _err_to_out=True)
    LOGGER.debug(log_string.getvalue())
    if not no_wps:
   	time.sleep(10 * 60)  # required 10 minute sleep
        # turn off peripherals via the web power switch
        LOGGER.info('Turning off peripherals')
        peripherals_off_cmd = sh.Command(peripherals_off)
        peripherals_off_cmd(_out=log_string,
                            _err_to_out=True)
        LOGGER.debug(log_string.getvalue())
        time.sleep(1)


def atnight(config_fname,
            delta_sunrise=DELTA_SUNRISE,
            delta_sunset=DELTA_SUNSET,
            peripherals_on=PERIPHERALS_ON_SCRIPT,
            peripherals_off=PERIPHERALS_OFF_SCRIPT,
            cdas=CDAS,
            cdas_down=CDAS_DOWN_SCRIPT,
            no_wps=False):
    """
    Given the PICASSO system configuration file *config_fname*
    determine if it is nighttime (sunset - *delta_sunset* [hours] <=
    current time <= sunrise + *delta_sunrise*) or daytime at the
    system location. Then, take action if the system needs to be
    turned on or off. Return `True` if the system is running at the
    completion of the function and `False` if it is off.

    Additional parameters:

    *peripherals_on*: script to turn on power to the system.

    *peripherals_off*: script to turn off power to the system.

    *cdas*: location of the cdas program.

    *cdas_down*: script to shutdown the cdas process.

    *no_wps*: indicate that there is no web power switch (and system
              power need not be turned on/off)
    """
    # get system location from configuration file
    site_info = get_site_info(config_fname)
    LOGGER.debug('Found site location in {}'.format(config_fname))
    LOGGER.debug('Location = {} '.format(site_info.location))
    LOGGER.debug('ID = {} '.format(site_info.id))
    LOGGER.debug('Latitude = {} [deg]'.format(site_info.lat))
    LOGGER.debug('Longitude = {} [deg]'.format(site_info.lon))
    LOGGER.debug('Altitude = {} [m]'.format(site_info.alt))
    # get current UTC time
    utc_dt = datetime.utcnow()
    LOGGER.info('Current time [UTC] = {:%Y-%m-%d %H:%M:%S}'.format(utc_dt))
    # determine if it is nighttime
    nighttime = is_nighttime(utc_dt,
                             site_info.lat,
                             site_info.lon,
                             site_info.alt,
                             delta_sunrise=timedelta(hours=delta_sunrise),
                             delta_sunset=timedelta(hours=delta_sunset))
    nighttime_str = 'YES' if nighttime else 'NO'
    LOGGER.info('Nighttime (sunset - {} [hours] <= current time <= sunrise + {} [hours])? {}'.format(delta_sunset,
                                                                                                     delta_sunrise,
                                                                                                     nighttime_str))
    # determine if the software system is running
    system_on = is_system_on()
    LOGGER.info('CDAS running? {}'.format('YES' if system_on else 'NO'))
    # take action depending on daytime/nighttime and system running/not running
    if nighttime and not system_on:
        # turn on the camera
        do_system_on(peripherals_on=peripherals_on,
                     peripherals_off=peripherals_off,
                     cdas=cdas,
                     config_fname=config_fname,
                     no_wps=no_wps)
        system_on = True
    elif not nighttime and system_on:
        # turn camera off
        do_system_off(peripherals_off=peripherals_off,
                      cdas_down=cdas_down,
                      no_wps=no_wps)
        system_on = False
    else:
        LOGGER.info('No action required')
    return system_on


def main(args):
    parser = ArgumentParser('Turn on/off PICASSO system based on location and time.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--config',
                        '-c',
                        type=str,
                        default=CONFIG_FNAME,
                        help='CDAS configuration file containing site latitude, longitude, and altitude.')
    parser.add_argument('--delta_sunrise',
                        type=float,
                        default=DELTA_SUNRISE,
                        help='Period (in hours) after sunset to turn on system.')
    parser.add_argument('--delta_sunset',
                        type=float,
                        default=DELTA_SUNSET,
                        help='Period (in hours) prior to sunset to turn on system.')
    parser.add_argument('--peripherals-on',
                        type=str,
                        default=PERIPHERALS_ON_SCRIPT,
                        help='Script to turn on web power switch outlets.')
    parser.add_argument('--peripherals-off',
                        type=str,
                        default=PERIPHERALS_OFF_SCRIPT,
                        help='Script to turn on web power switch outlets.')
    parser.add_argument('--cdas',
                        type=str,
                        default=CDAS,
                        help='Location of cdas program.')
    parser.add_argument('--cdas-down',
                        type=str,
                        default=CDAS_DOWN_SCRIPT,
                        help='Location of script to turn off cdas.')
    parser.add_argument('--no-wps',
                        action='store_true',
                        help='Do not send on/off to web power switch (i.e., do not run peripherals on/off scripts).')
    test_group = parser.add_mutually_exclusive_group()
    test_group.add_argument('--test-on',
                            action='store_true',
                            help='Test the power on sequence and exit.')
    test_group.add_argument('--test-off',
                            action='store_true',
                            help='Test the power off sequence and exit.')
    parser.add_argument('--log',
                        type=str,
                        default=LOG_FNAME,
                        help='Log file location.')
    args = parser.parse_args(args)

    logging.basicConfig(filename=args.log,
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    # squelch sh module output
    logging.getLogger('sh').setLevel(logging.WARNING)


    if args.test_on:
        # test power on sequence and exit
        do_system_on(peripherals_on=args.peripherals_on,
                     peripherals_off=args.peripherals_off,
                     cdas=args.cdas,
                     config_fname=args.config,
                     no_wps=args.no_wps)
        sys.exit(0)

    if args.test_off:
        # test power on sequence and exit
        do_system_off(peripherals_off=args.peripherals_off,
                      cdas_down=args.cdas_down,
                      no_wps=args.no_wps)
        sys.exit(0)

    atnight(args.config,
            delta_sunrise=args.delta_sunrise,
            delta_sunset=args.delta_sunset,
            peripherals_on=args.peripherals_on,
            peripherals_off=args.peripherals_off,
            cdas=args.cdas,
            cdas_down=args.cdas_down,
            no_wps=args.no_wps)


if __name__ == '__main__':
    main(sys.argv[1:])
