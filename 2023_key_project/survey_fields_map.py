import csv
import healpy as hp
import healpixel_functions
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt

NSIDE = 64
NPIX = hp.nside2npix(NSIDE)

def read_survey_field_list(file_path):
  fields = []
  with open(file_path, 'r', newline='') as csvfile:
      reader = csv.reader(csvfile, delimiter=',', quotechar='|')
      for row in reader:
          s = SkyCoord(row[0], row[1], frame='icrs', unit=(u.hourangle, u.deg))
          fields.append(row + [s])

  return fields

# Survey definitions.  Parameters are:
# Field locations, radius of field [deg]
survey_regions = {
    'ogle': ['ogle/ogle_survey_fields.csv', np.sqrt(1.4/np.pi)],
    'kmtnet': ['kmtnet/kmtnet_survey_fields.csv', np.sqrt(4.0/np.pi)]
}

map = np.zeros(NPIX)

for survey, survey_params in survey_regions.items():
    field_centres = read_survey_field_list(survey_params[0])
    #print(field_centres)
    hp_index = []
    for coord in field_centres:
        pixels = healpixel_functions.skycoord_to_HPindex(coord[2], NSIDE,
                                                         radius=survey_params[1])
        hp_index += pixels.tolist()
    map_pixels = list(set(hp_index))

    map[map_pixels] += 1.0

fig = plt.figure(1,(10,10))
hp.mollview(map)
hp.graticule()
plt.tight_layout()
plt.savefig('mulens_survey_fields.png')
plt.close(1)
