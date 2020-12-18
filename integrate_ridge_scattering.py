#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 4 10:37:04 2018

@author: juanenciso
"""

import argparse
from datetime import datetime
import numpy as np
import fabio
from scipy.signal import argrelmax
import peakutils
from misc_funcs import simplify_header
from misc_funcs import configure_azimuthal_integrator as cf_az_int

# find_ridge_scattering(fabio.obj, np.array) -> (np.int64, np.int64)
def find_ridge_scattering(par_fabio_obj, par_mask=None):
    """
    Finds lower and upper bounds for the scattering of the ridge, we
    assume the largest q possible will be due to scattering of the dorsal scales
    which are likely to have the smallest ridge spacing

    As not all experiments will yield the same levels of intensity, it
    should be changed accordingly using

    Parameters:
    -----------

    - A fabio object containing a 2D scattering pattern
    - A np.array of shape equal to the shape of the fabio object data with
    masked pixels as 1 and normal pixels as 0 (default = None)

    Returns:
    --------

    (lower, upper), a tuple of np.in64 integers corresponding to the lower and
    upper bounds for the azimuthal angle of ridge scattering.

    This function is expected to work if par_fabio_obj contains a well formed
    scattering pattern from a butterfly scale (orientated sample)
    """
    az_int = cf_az_int(par_fabio_obj)
    qscat, _ = az_int.integrate1d(par_fabio_obj.data, NPOINTS,
                                  azimuth_range=(0, WINDOW),
                                  radial_range=(L_Q, U_Q),
                                  mask=par_mask)
    lb_array = np.arange(PF_START, PF_STOP - WINDOW + 1, STEP)
    ub_array = lb_array + WINDOW
    indices = np.arange(0, lb_array.size)
    q_tv = np.full_like(lb_array, 0.1, dtype=float)
    i_tv = np.full_like(lb_array, 10.0, dtype=float)
    for low, upp, index in zip(lb_array, ub_array, indices):
        _, intens = az_int.integrate1d(par_fabio_obj.data, NPOINTS,
                                       azimuth_range=(low, upp),
                                       radial_range=(L_Q, U_Q),
                                       mask=par_mask)
        q_tv[index] = qscat[np.argmax(intens)]
        i_tv[index] = intens.max()
    pkindexes = peakutils.indexes(i_tv, thres=0.25, min_dist=3)
    q_top_idx = pkindexes[np.argmax(q_tv[pkindexes])]
    return (lb_array[q_top_idx], ub_array[q_top_idx])


# valid_frame(np.array, np.array, float) -> bool
def valid_frame(par_q, par_i, par_thresh_intensity=1000.0):
    """
    This function validates that the frame contains useful scattering data

    Parameters:
    -----------

    par_q np.array with values of a scattering vector
    par_i np.array with values of scattered intensity
    par_thresh_intensity float containing the threshold of minimum intensity
    (default = 1000.0)

    Returns:
    --------

    True if all conditions for a 'good' frame are satisfied, False otherwise
    """
    # ======================================================
    # 1. Highest intensity needs to be at least 3000.0
    # 2. There should be at least two q peaks identified by
    #    argrelmax(intensity, order=ORDER_AMAX)
    # 3. The error between predicted and estimated Bragg peaks
    #    should not be larger than 1e-04
    # ======================================================

    max_values = argrelmax(par_i, order=ORDER_AMAX)[0]

    # Verify condition 1
    if par_i.max() <= par_thresh_intensity:
        print("1. Intensity is too low!")
        return False
    # Verify condition 2
    if max_values.size < 2:
        print("2. Not enough peaks!")
        return False
    # Remove spurious peaks
    if par_i[max_values][0] < par_i[max_values][1]:
        max_values = np.delete(max_values, 0)
    predicted_bragg = par_q[np.argmax(par_i)] * np.array([1., 2.])
    estimated_bragg = par_q[max_values[:2]]
    err_bragg = np.array([1e-04, 1e-03])
    # Verify condition 3
    if not np.all(np.abs(predicted_bragg - estimated_bragg) < err_bragg):
        print("3. No Bragg scattering!")
        return False

    return True


# integrate_angles(int, int, fabio.obj) -> (np.array, np.array)
def integrate_angles(par_l, par_u, par_fabio_obj, par_mask=None):
    """
    Given lower and upper azimuthal angles, integrate between these and
    return q (scat) and i (intens) vectors

    Parameters:
    -----------
    par_l lower bound for azimuthal angle
    par_u upper bound for azimuthal angle
    par_fabio_obj fabio object containing a 2D scattering pattern

    Returns:
    --------
    (scat, intens) a tuple of arrays
    scat contains the scattering vector, intens contains scattered intensity
    """
    az_int = cf_az_int(par_fabio_obj)
    scat, intens = az_int.integrate1d(par_fabio_obj.data, NPOINTS,
                                      azimuth_range=(par_l, par_u),
                                      radial_range=(L_Q, U_Q),
                                      mask=par_mask)
    return (scat, intens)


def main():
    """
    Main execution block
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_edf", type=str,
                        help=("Path to the sample edf file to be analysed."))
    parser.add_argument("empty_edf", type=str,
                        help=("Path to the edf file of the empty beam used "
                              "to apply intensity corrections."))
    parser.add_argument("mask_file", type=str,
                        help=("Path to the text file containing the masked"
                              " pixels of the beamstop. The value of the "
                              "un-masked pixels should be 1 and the value of"
                              "the masked pixels should be 0, as in a Foxtrot"
                              " mask file."))
    args = parser.parse_args()

    mask_array = np.loadtxt(args.mask_file, delimiter=";") * -1 + 1
    empty_edf = fabio.open(args.empty_edf)
    sample_edf = fabio.open(args.sample_edf)
    ind_name = simplify_header(sample_edf)['IndName']
    frame_name = args.sample_edf.split("/")[-1].rstrip("_raw.edf")
    l_az, u_az = find_ridge_scattering(sample_edf, par_mask=mask_array)
    # Find the ridge bits
    r_q, r_i = integrate_angles(l_az, u_az, sample_edf, par_mask=mask_array)
    # Validate frame
    if valid_frame(r_q, r_i):
        # If valid, get cross rib intensity
        # In theory angle should be 90  but in test 85 yields a clearer signal
        cr_q, cr_i = integrate_angles(l_az + 85, u_az + 85, sample_edf,
                                      par_mask=mask_array)
        # Get intensity for empty_beam at specific angles
        _, emp_r_i = integrate_angles(l_az, u_az, empty_edf)
        _, emp_cr_i = integrate_angles(l_az + 85, u_az + 85, empty_edf,
                                       par_mask=mask_array)
        # Apply background correction
        br_r_i = r_i - emp_r_i
        br_cr_i = cr_i - emp_cr_i
        # Dump to file
        out_name_r = "{}.{}.RDV.txt".format(ind_name, frame_name)
        np.savetxt(out_name_r, np.transpose((r_q, br_r_i)),
                   delimiter='\t',
                   fmt='%.8f')
        out_name_cr = "{}.{}.CRV.txt".format(ind_name, frame_name)
        np.savetxt(out_name_cr, np.transpose((cr_q, br_cr_i)),
                   delimiter='\t',
                   fmt='%.8f')
    else:
        print("{} FAILED!".format(args.sample_edf))
        nv_file = "{}.{}.NonValid.txt".format(ind_name, frame_name)
        with open(nv_file, 'w') as non_valid:
            nowtime = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
            non_valid.write("{}\t{}\n".format(nowtime,
                                              args.sample_edf.split("/")[-1]))

# =============================================================================
# Constants defining the boundaries of search range
# WARNING: These settings work ok for experiments with a camera
# to detector distance of ~31m. They need to be revised/changed for other
# experimental settings.
# =============================================================================

# CONSTANTS
PF_START = -9
PF_STOP = 171
WINDOW = 5
STEP = 1
L_Q = 0.003
U_Q = 0.03
NPOINTS = 1000
ORDER_AMAX = 100

if __name__ == "__main__":
    main()
