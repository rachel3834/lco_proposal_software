import numpy as np
import h5py
from os import path
from astropy.time import Time
from scipy import interpolate

def output_visibility_data(NPIX, dates, visibility_table, file_path):

    jd = np.array([(ts.decimalyear-int(ts.decimalyear)) for ts in dates])
    site_codes = {'LSC': 0, 'CPT': 1, 'COJ': 2, 'ELP': 3, 'TFN': 4, 'OGG': 5}
    sites = np.array([site_codes[x] for x in observatories.keys()])

    with h5py.File(file_path, "w") as f:

        dset1 = f.create_dataset('healpix',
                                (NPIX),
                                dtype='int',
                                data=np.arange(1,NPIX+1,1))

        dset2 = f.create_dataset('decimalyear',
                                (len(dates)),
                                dtype='float64',
                                data=jd)

        dset3 = f.create_dataset('observatories',
                            (len(sites)),
                            dtype='int',
                            data=sites)

        dset4 = f.create_dataset('visibility_data',
                            visibility_table.shape,
                            dtype='float64',
                            data=visibility_table)

        dset5 = f.create_dataset('total_hrs_visible',
                            total_hrs_visible.shape,
                            dtype='float64',
                            data=total_hrs_visible)
        f.close()

def read_visibility_data(file_path):

    # The list of observatory sites is encoded as an integer array,
    # because HD5 doesn't handle string arrays well.
    site_codes = {0: 'LSC', 1: 'CPT', 2: 'COJ', 3: 'ELP', 4: 'TFN', 5: 'OGG'}

    if not path.isfile(file_path):
        raise IOError('Cannot find visibility data file '+file_path)

    dataset_list = ['healpix', 'decimalyear', 'observatories',
                    'visibility_data', 'total_hrs_visible']

    f = h5py.File(file_path, "r")

    visibility_data = {}
    for dataset in dataset_list:
        dset = f[dataset]
        if dataset == 'observatories':
            codes = np.array(dset[:])
            visibility_data[dataset] = np.array([site_codes[x] for x in codes])
        else:
            visibility_data[dataset] = np.array(dset[:])

    return visibility_data

def calc_network_hrs_visible(visibility_data, ipix):
    """Function to calculate the total number of hours that a given HEALpixel
    can be observed for, combining data from different sites.

    This is calculated separately for telescopes in the northern and
    southern hemispheres, since the sites can observe sequentially.

    total_hrs_visible = t_north + t_south,

    where t_north = t_TFN + t_ELP
          t_south = t_COJ + t_CPT + t_LSC

    Note the input HEALpixel identifier should be the array index, NOT
    the HEALpixel number
    """

    # Set up storage array, and identify the array indices of the telescopes in
    # the northern and southern hemispheres.
    data = np.zeros((len(visibility_data['dates']),2))

    observatories = list(visibility_data['observatories'])
    southern_ring = [observatories.index('LSC'),
                     observatories.index('CPT'),
                     observatories.index('COJ')]
    northern_ring = [observatories.index('ELP'),
                     observatories.index('TFN')]

    # Extract the visibility data for the requested pixel for
    # all sites as a 2D array:
    hrs_visible = visibility_data['hrs_visible'][ipix,:,:]

    for idate in range(0,len(visibility_data['dates']),1):
        t_north = hrs_visible[idate, northern_ring].sum()
        t_south = hrs_visible[idate, southern_ring].sum()

        data[idate,0] = visibility_data['dates'][idate]
        data[idate,1] = min((t_north + t_south),24.0)

    return data

def interpolate_visibility(hp_visibility, dates):
    """Function to create an interpolate function object for an array
    containing the total_hrs_visible for a given HEALpixel, as a function
    of an array of JD dates."""

    return interpolate.interl1d(hp_visibility, dates)
