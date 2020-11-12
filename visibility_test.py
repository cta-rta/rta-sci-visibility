# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import argparse
import warnings
import yaml
import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from lib.visibility_tools import Visibility, complete_irf_name


# parse command line inputs
parser = argparse.ArgumentParser(description='Simple test for source visibility.')
parser.add_argument('config', metavar='cf', type=str, help='configuration yaml file')
# configuration file
cf = parser.parse_args().config
# load params configuration from cf
with open('config_visibility.yaml') as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)

filename = cfg['path']['filename']
site = cfg['setup']['site']
twilight = cfg['setup']['twilight']
sites = cfg['sites_list']
thresholds = cfg['setup']['thresholds']
zenith_angles = cfg['setup']['zenith_angles']
use_visibility_table = cfg['use_visibility_table']
if not use_visibility_table:
    total_points = cfg['total_points']
window_points = cfg['window_points']

# initialise site coordinates
site_coords = EarthLocation.of_site(sites[site])
# load template
with fits.open(filename) as hdul:
    hdr = hdul[0].header
    # source coordinates
    source_radec = SkyCoord(ra=hdr['RA'] * u.deg, dec=hdr['DEC'] * u.deg, frame='icrs')
    # source trigger time
    t_trigger = Time(hdr['GRBJD'] * u.day, format='jd')
    if use_visibility_table:
        # using visibility table
        visibility_table = Table.read(hdul, hdu=1)
        # time windows [THIS SHOULD BE CHANGED ACCORDING TO THE TEMPLATE FORMAT]
        t = Time(visibility_table['True'].data, format='jd')
        t_true = {'North': t[0:2], 'South': t[2:4]}
        # minimum altitude
        min_altitude = visibility_table.meta['MIN_ALT'] * u.deg
        # replace minimum altitude within thresholds
        thresholds[0] = min_altitude.value
    else:
        # otherwise use event total duration
        times = np.array(hdul['TIMES (AFTERGLOW)'].data.tolist())
        afterglow_duration = Time(((times[-1] + times[1]) / 2)[0] / 86400, format='jd')


# ------------------------- EXAMPLE 1 (COMPACT) :: USING VISIBILITY TABLE -------------------------- !!!

if use_visibility_table:
    print('Example of use: visibility table from template')
    # set start time and duration
    t_start = t_true[site][0]
    duration = t_true[site][1] - t_true[site][0]
    # ignore warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        # initialise
        visibility = Visibility()
        # short-cut: IRFs and relative time intervals
        irfs = visibility.associate_irf_one_night(source_radec, t_start, duration, sites[site], window_points, thresholds, zenith_angles)
        # complete IRFs
        irfs['irf'] = complete_irf_name(irfs['irf'], site, '0.5h')
        print('IRFs', irfs)
        del visibility

# ------------------------------ EXAMPLE 2 (EXPLICIT) :: W/O USING VISIBILITY TABLE -------------------- !!!
else:
    print('\nExample of use: event full duration')
    # set start time and duration
    t_start = t_trigger
    duration = afterglow_duration
    # ignore warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        # initialise
        visibility = Visibility()
        # visibility points in JD and AltAz
        visibility.visibility_points(t_start, duration, total_points)
        visibility.visibility_altaz(source_radec, sites[site])
        # find nights
        nights = visibility.get_nighttime(twilight=twilight)
        print('nights:', nights)
        del visibility
        # within each night find IRFs
        for i in range(len(nights['start'])):
            print('Night', i+1, 'of', len(nights['start']))
            t_start = Time(nights['start'][i], format='jd')
            duration = Time(nights['stop'][i] - nights['start'][i], format='jd')
            # initialise
            visibility = Visibility()
            # visibility points in JD and AltAz
            visibility.visibility_points(t_start, duration, window_points)
            visibility.visibility_altaz(source_radec, sites[site])
            # IRFs and relative time intervals
            irfs = visibility.associate_irf_zenith_angle(thresholds, zenith_angles)
            # complete IRFs
            irfs['irf'] = complete_irf_name(irfs['irf'], site, '0.5h')
            print('IRFs', irfs)
            del visibility

# ------------------------------ EXAMPLE 3 (IMPLICIT) :: PROLONGED EVENT DURATION -------------------- !!!

t_start = t_trigger
# The template duration was too short to test >1 windows of observability
duration = Time(3, format='jd')
# ignore warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore')
    print('\nEvent full duration')
    # initialise
    visibility = Visibility()
    # visibility points in JD and AltAz
    visibility.visibility_points(t_start, duration, total_points)
    visibility.visibility_altaz(source_radec, sites[site])
    # find nights
    nights = visibility.get_nighttime(twilight=twilight)
    print('nights:', nights)
    del visibility
    # within each night find IRFs
    for i in range(len(nights['start'])):
        print('Night', i+1, 'of', len(nights['start']))
        t_start = Time(nights['start'][i], format='jd')
        duration = Time(nights['stop'][i] - nights['start'][i], format='jd')
        # initialise
        visibility = Visibility()
        # short-cut: IRFs and relative time intervals
        irfs = visibility.associate_irf_one_night(source_radec, t_start, duration, sites[site], window_points, thresholds, zenith_angles)
        # complete IRFs
        irfs['irf'] = complete_irf_name(irfs['irf'], site, '0.5h')
        print('IRFs', irfs)
        del visibility