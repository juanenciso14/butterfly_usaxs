#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 10:08:39 2018

@author: juanenciso
"""

import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import pyFAI

# Get filenames
def get_filenames(par_pattern, par_path=None):
    """
    get_filenames(str) -> list[str]

    Parameters:
    -----------
    An extension so that the program can look for files with such extension
    and make a list of those files

    Returns:
    --------
    A list with the names of the files containing the data that is going to be
    in the plot

    Example:
    --------
    get_filenames('.txt') -> ['lol.txt', 'haha.txt']
    """
    if par_path is None:
        the_path = os.getcwd()
    else:
        the_path = par_path
    return list(filter(lambda x: par_pattern in x, os.listdir(the_path)))


def configure_azimuthal_integrator(par_edf_hanldle):
    """
    A wrapper around a pyFAI.AzimuthalIntegrator constructor and its setFit2D
    method that returns a properly configured integrator for a specific edf
    handle. The function assumes that the detector has no tilt and that the
    beam center is equal to the point of normal incidence.

    configure_azimuthal_integrator(edf_handle) -> AzimuthalIntegrator
    """
    # Get header attributes and convert to numeric
    header = par_edf_hanldle.header
    centre1px = float(header['Center_1'])
    centre2px = float(header['Center_2'])
    pixel1m = float(header['PSize_1'])
    pixel2m = float(header['PSize_2'])
    distance_m = float(header['SampleDistance'])
    wavelength_m = float(header['WaveLength'])

    # First configuration
    my_azin = pyFAI.AzimuthalIntegrator(
            dist=distance_m,
            poni1=centre1px*pixel1m,
            poni2=centre2px*pixel2m,
            wavelength=wavelength_m)
    # Second configuration
    my_azin.setFit2D(directDist=distance_m*1e3, centerX=centre1px,
                     centerY=centre2px,
                     tilt=0.0, tiltPlanRotation=0.0,
                     pixelX=pixel1m*1e6, pixelY=pixel2m*1e6)
    return my_azin

# Get edf header data
def simplify_header(par_edf_handle):
    """
    Return a dictionary containing only centre coordinates, pixel sizes,
    wavelength, sample-detector distance and individual name
    """
    # Get header
    header = par_edf_handle.header
    centre1px = float(header['Center_1'])
    centre2px = float(header['Center_2'])
    pixel1m = float(header['PSize_1'])
    pixel2m = float(header['PSize_2'])
    distance_m = float(header['SampleDistance'])
    wavelength_m = float(header['WaveLength'])
    ind_name = header['Title']
    return {"Center_1": centre1px,
            "Center_2": centre2px,
            "PSize_1": pixel1m,
            "PSize_2": pixel2m,
            "SampleDistance": distance_m,
            "WaveLength": wavelength_m,
            "IndName": ind_name}


def theta_radians(parq, parlambda):
    """
    Returns a measure of theta in radians given q and wavelength
    parq and parlambda are in nanometres
    """
    return np.arcsin(parlambda*parq/4*np.pi)


def metres_nanometres(par_m):
    """
    Converts metres to nanometres
    """
    return par_m * 1e9


def metres_micrometres(par_m):
    """
    Converts metres to micrometres
    """
    return par_m * 1e9


def pix_to_q(par_psize, par_npix, sd_dist, wavelength):
    """
    Arguments:
    ----------
    par_psize size of pixel in metres
    par_npix number of pixels
    sd_dist is the sample to detector distance in metres
    wavelength is the beam wavelength in metres
    PSize_1 is in the x axis
    PSize_2 is in the y axis
    And then use circle parameterization with x and y coords to get
    the circles drawn as plt.contour on top of 2D patterns
    Give also centre coordinates
    https://stackoverflow.com/questions/32092899/plot-equation-showing-a-circle

    Test values:
    ------------
    PSize_1: 2.36285e-05 m
    PSize_2: 2.39661e-05 m
    WavewLength: 9.95058e-11 m
    SampleDistance: 30.9883 m
    # Centre coordinates provided in case they are needed
    Center_1: 958.376 px
    Center_2: 939 px
    # TODO MAKE SURE THIS IS PRODUCING CORRECT RESULTS
    """
    pix_dist = par_psize * par_npix
    wl_nm = metres_nanometres(wavelength)
    sd_mu = metres_micrometres(sd_dist)
    px_mu = metres_micrometres(pix_dist)
    # q = 4*pi*sin(2*theta)/wl_nm
    # tan(2*theta) = px_mu/sd_mu
    theta = np.arctan(px_mu/sd_mu)/2
    q_nm = (4*np.pi*np.sin(theta))/wl_nm
    return q_nm


def q_to_pix(par_psize1, par_psize2, par_q, sd_dist, wavelength):
    """
    Converts q [nm-1] to number of pixels
    Can be improved with *args and **kwargs
    par_psize1 is the pixel size in x axis
    par_psize2 is the pixel size in y axis
    par_q is the scattering vector in inverse nanometres
    sd_dist is the sample to detector distance in metres
    wavelength is the beam wavelength in metres

    Return:
    -------

    A tuple of pixel coordinates (x, y)
    """

    wl_nm = metres_nanometres(wavelength)
    # Formula:
    # sd_dist (m) * tan(arcsin(wl*q/4*pi)) = pxdist (m)
    # q = 4*pi*sin(2*theta)/wl_nm
    # theta = 2*arcsin(q*wl_nm/4*pi)
    theta = 2*np.arcsin((wl_nm*par_q)/(4*np.pi))
    px_dist = sd_dist*(np.tan(theta))
    npix1 = int(np.ceil(px_dist/par_psize1))
    npix2 = int(np.ceil(px_dist/par_psize2))
    return (npix1, npix2)


def d_to_pix(par_psize1, par_psize2, par_d, sd_dist, wavelength):
    """
    Wrap around q_to_pix to get number of pixels in x and y for a
    particular distance par_d measured in nanometres
    """
    q_nm = 2*np.pi/par_d
    return q_to_pix(par_psize1, par_psize2, q_nm, sd_dist, wavelength)
