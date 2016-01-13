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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple
from datetime import datetime, timedelta

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


def system_on():
    """
    Turn the camera system on.
    """
    pass


def system_off():
    """
    Turn the camera system off.
    """
    pass


def atnight(config_fname,
            delta_sunrise=DELTA_SUNRISE,
            delta_sunset=DELTA_SUNSET):
    """
    Given the PICASSO system configuration file *config_fname*
    determine if it is nighttime (sunset - *delta_sunset* [hours] <=
    current time <= sunrise + *delta_sunrise*) or daytime at the
    system location. Then, take action if the system needs to be
    turned on or off. Return `True` if the system is running at the
    completion of the function and `False` if it is off.
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
        LOGGER.info('Turning system ON')
        system_on()
        system_on = True
    elif not nighttime and system_on:
        # turn camera off
        LOGGER.info('Turning system OFF')
        system_off()
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
    parser.add_argument('--cdas',
                        type=str,
                        default=CDAS_PROCESS_NAME,
                        help='Process name for PICASSO system (used to determine if the system is on or off)')
    parser.add_argument('--log',
                        type=str,
                        default=LOG_FNAME,
                        help='Log file location.')
    args = parser.parse_args(args)

    logging.basicConfig(filename=args.log,
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    atnight(args.config,
            delta_sunrise=args.delta_sunrise,
            delta_sunset=args.delta_sunset)


if __name__ == '__main__':
    main(sys.argv[1:])
