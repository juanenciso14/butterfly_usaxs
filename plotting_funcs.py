#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 11:25:34 2018

@author: juanenciso
"""
# %% Append current directory to path so custom modules can be imported
import sys
sys.path.append('.')
# %% Other imports

# built in imports
import re

# third party imports
import fabio
import matplotlib.pyplot as plt
import numpy as np

# local imports
import misc_funcs
# %% Fns


def plot_single_2d(par_edf_data):
    """
    # TODO: Docstring this function
    """
    my_fig = plt.figure()
    my_ax = my_fig.gca()
    # Transform data with arcsinh
    tr_edf = np.arcsinh(par_edf_data)
    low, up = np.min(tr_edf), np.max(tr_edf)
    lspace = np.linspace(low, up, 50)
    my_ax.imshow(tr_edf, origin="lower",
                 aspect="auto", cmap="jet",
                 vmin=lspace[3],
                 vmax=lspace[-12])


# TODO: Modify so that I can pass a dictionary of keyword arguments in **kwargs
# TODO: Modify so that I pass only a par_edf_list that should contain the full
# path to access edf files
def plot_multi_2d(par_edf_list, par_path, par_nrows, par_ncols,
                  par_dist, cmap="jet", show=False, title_override=True):
    """
    Parameters:
    -----------

    par_edf_list is a list with the names of the files that are going to
    be used to plot.

    par_path is the path to the edf files to process.

    par_nrows is the number of rows that the pdf file must have.

    par_ncols is the number of columns that the pdf file must have.

    par_dist is the sample to detector distance of the scattering experiment.

    cmap is the colormap to be used (default "jet").

    show whether to display the plots in a window (default False).

    title_override whether to override the automatic title given to the plot
    by matplotlib's pdf interface. When True, the name of the sample will be
    used in the title (default True).
    """
    edf_generator = ((fabio.open(par_path + edf),
                      edf) for edf in par_edf_list)
    # Compile regex first to save time
    rgx = re.compile(r'ls\d+_frelon_(.+)_raw.edf')
    ind_name = "NA"
    fig, ax = plt.subplots(par_nrows, par_ncols,
                           figsize=(par_ncols, par_nrows))
    fig.subplots_adjust(hspace=0, wspace=0)
    for i in range(par_nrows):
        for j in range(par_ncols):
            try:
                curr_edf, edf_name = next(edf_generator)
                # Get header details
                simple_head = misc_funcs.simplify_header(curr_edf)
                # simple_head has:
                # Center_1 Center_2 (px) PSize_1 PSize_2 (m)
                # SampleDistance (m) WaveLength (m) Title
                cent1 = int(np.rint(simple_head['Center_1']))
                cent2 = int(np.rint(simple_head['Center_2']))
                px1, px2 = misc_funcs.d_to_pix(simple_head['PSize_1'],
                                               simple_head['PSize_2'],
                                               par_dist,
                                               simple_head['SampleDistance'],
                                               simple_head['WaveLength'])
                # Get name of individual for title
                ind_name = simple_head['IndName']
                xmin, xmax = cent1 - px1, cent1 + px1
                ymin, ymax = cent2 - px2, cent2 + px2
                # Get frame name
                frame_id = rgx.sub(r'\1', edf_name)
                curr_edf = curr_edf.data[np.ix_(np.arange(ymin, ymax),
                                                np.arange(xmin, xmax))]
                # Commented lines below are useful to modify contrast of plots
                # tr_edf = np.arcsinh(curr_edf)
                # low, up = np.min(tr_edf), np.max(tr_edf)
                # lspace = np.linspace(low, up, 50)
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
                ax[i, j].yaxis.set_major_locator(plt.NullLocator())
                # ax[i, j].imshow(np.arcsinh(curr_edf), origin="lower",
                # aspect="auto", cmap=cmap, vmin=lspace[3], vmax=lspace[-12])
                ax[i, j].imshow(np.arcsinh(curr_edf), origin="lower",
                                aspect="auto", cmap=cmap)
                ax[i, j].text(0.25, 0.1, frame_id,
                              bbox={'facecolor': 'white',
                                    'edgecolor': 'none', 'pad': 2},
                              transform=ax[i, j].transAxes, fontsize=3,
                              ha='center', va='center')
            # Generator exhausted
            except StopIteration:
                ax[i, j].set_axis_off()
    if not title_override:
        fig.suptitle(ind_name, va='top')
    if show:
        plt.show()

# TODO: Add plotting 1D utilities

# %% Test
# TODO: DELETE LATER
# Oct 18 multiplot some erato and melpomene
# path_to_edf = ("/media/juanenciso/HD-LXU3/erato18009_FRAMES_A0/")
# Make a list of files in path
# edf_list = misc_funcs.get_filenames('.edf', path_to_edf)
# Decide which number of rows to use if we are to use only 5 columns
# nrows = (len(edf_list) // 5) + 1
# Make a generator to call when plotting each image
# edf_files = (fabio.open(path_to_edf + edf) for edf in edf_list)
# %% Plot several 2D edfs
# plot_multi_2d(edf_list, path_to_edf, 6, 300, cmap="jet")
