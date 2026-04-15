import csv
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, Column
import numpy as np

def read_event_list(event_file_path, coords='sexigesimal', ref_frame='icrs'):
    data = []
    with open(event_file_path, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            if coords == 'decimal_deg':
                s = SkyCoord(row[1], row[2], frame=ref_frame, unit=(u.deg, u.deg))
            else:
                s = SkyCoord(row[1], row[2], frame=ref_frame, unit=(u.hourangle, u.deg))
            data.append(row + [s])
    data = np.array(data)

    event_table = Table(
                        [Column(data=data[:,0], name='Event'),
                        Column(data=data[:,1], name='RA'),
                        Column(data=data[:,2], name='Dec'),
                        Column(data=data[:,3], name='Coordinate')]
                        )

    return event_table
