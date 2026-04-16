import healpy as hp
import csv
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpixel_functions
from sys import argv

def read_event_list(event_file_path):
    data = []
    with open(event_file_path, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            s = SkyCoord(row[1], row[2], frame='icrs', unit=(u.hourangle, u.deg))
            data.append(row + [s])
    data = np.array(data)

    event_table = Table(
                        [Column(data=data[:,0], name='Event'),
                        Column(data=data[:,1], name='RA'),
                        Column(data=data[:,2], name='Dec'),
                        Column(data=data[:,3], name='Coordinate')]
                        )

    return event_table


def read_event_catalog(file_path):
    column_types = {
        'Name': 'str',
        'RA': 'str',
        'Dec': 'str',
        't0[HJD]': 'float',
        'tE[days]': 'float',
        'u0': 'float',
        'Ibase': 'float',
        'Ibase_err': 'float'
    }
    df = pd.read_table(file_path, sep='\s+', dtype=column_types)

    return df

def output_survey_fields_csv(survey_field_centres, output_file_prefix):
    output_file_name = output_file_prefix+'_survey_fields.csv'
    with open(output_file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for entry in survey_field_centres:
            writer.writerow([entry[0], entry[1]])

if len(argv) == 1:
    event_file_path = input('Please enter the path to the catalog of events: ')
    output_file_prefix = input('Please enter the prefix for the output files: ')
    nyears_survey = input('Please enter the number of years of survey: ')
else:
    event_file_path = argv[1]
    output_file_prefix = argv[2]
    nyears_survey = argv[3]
nyears_survey = float(nyears_survey)

NSIDE = 64
NPIX = hp.nside2npix(NSIDE)

# Read the table of events
event_table = read_event_catalog(event_file_path)

# Convert the list of event coordinates to HEALpixel indices,
# and from this extract the unique HEALpixel IDs where the events
# occurred.  Two calls to this function are used here.  The first
# looks at the sky regions the survey is pointing at.  The second
# counts how many events occur within each specific HEALpixel.
hp_index = []
event_hpindices = []
event_map = np.zeros(NPIX)

event_coords = SkyCoord(event_table['RA'], event_table['Dec'], frame='icrs', unit=(u.hourangle, u.deg))

for coord in event_coords:
    survey_pixels = healpixel_functions.skycoord_to_HPindex(coord, NSIDE)
    hp_index += survey_pixels.tolist()
    event_map[survey_pixels[0]] += 1.0
map_pixels = list(set(hp_index))

# Plot the map of survey fields
map = np.zeros(NPIX)
map[map_pixels] = 1.0
fig = plt.figure(1,(10,10))
hp.mollview(map)
hp.graticule()
plt.tight_layout()
plt.savefig(output_file_prefix+'_survey_fields.png')
plt.close(1)

# Plot the map of event rate
fig = plt.figure(1,(10,10))
hp.mollview(event_map/nyears_survey)
hp.graticule()
plt.tight_layout()
plt.savefig(output_file_prefix+'_event_rate_map.png')
plt.close(1)

# Output the list of survey fields
survey_field_centres = []
for pixel in map_pixels:
    field_centre = healpixel_functions.HPindex_to_skycoord(pixel, NSIDE)
    (ra,dec) = field_centre.to_string('hmsdms', sep=':')[0].split()
    survey_field_centres.append([ra, dec])
output_survey_fields_csv(survey_field_centres, output_file_prefix)

# Output a FITS table of the event rate data:
event_rate = event_map/nyears_survey
hp_fields_ra = []
hp_fields_dec = []
for ihp in range(0,len(event_map),1):
    field_centre = healpixel_functions.HPindex_to_skycoord(ihp, NSIDE)
    (ra,dec) = field_centre.to_string('hmsdms', sep=':')[0].split()
    hp_fields_ra.append(ra)
    hp_fields_dec.append(dec)

hdr = fits.Header()
hdr['NSIDE'] = NSIDE
hdr['NPIX'] = NPIX
hdr['ORDER'] = 'RING'
hdr['MAPTITLE'] = 'Event rate'
phdu = fits.PrimaryHDU(header=hdr)

c1 = fits.Column(name='hp_index', array=range(1,len(map_pixels),1), format='I8')
c2 = fits.Column(name='RA', array=hp_fields_ra, format='12A')
c3 = fits.Column(name='Dec', array=hp_fields_dec, format='12A')
c4 = fits.Column(name='gamma', array=event_rate, format='E')
hdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4])

hdul = fits.HDUList([phdu,hdu])
hdul.writeto(output_file_prefix+'_event_rate.fits', overwrite=True)
