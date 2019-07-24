#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 19:13:55 2018

@author: juanenciso
"""

# %% Append current directory to path so modules can be imported
import sys
sys.path.append('.')
# %%

# builtin module imports
# third party imports
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
# local imports
import plotting_funcs
import misc_funcs
import matplotlib.pyplot as plt

# %% Functions


def multiple_pages(param_path=None, param_list=None, par_nrows=5, par_ncols=4,
                   par_dist=300, par_dpi=300, title_override=False):
    """
    Prints several edf frames in pdf files. The frames can be organised
    by individual in a single folder or they can come as a list that contains
    the full path plus edf names.

    TODO: Fully document function
    """
    # Assert that both param_path and param_list can't be both None and can't
    # be both different to None
    if param_path is not None and param_list is not None:
        raise ValueError(("Use either a path to edfs via param_path "
                          "or a list of edf files including their path "
                          "to call multiple_pages function"))
    elif param_path is None and param_list is None:
        raise ValueError(("Either param_path or param_list has to be "
                          "different from None when calling multiple_pages"))

    frames_per_page = par_nrows * par_ncols

    # Prefix should change dependig on whether a path containing all
    # the frames in one individual is passed or whether a list of single file
    # locations is passed
    # Name of file
    if param_path is not None:
        prefix = param_path.split("/")[-2]
        list_edfs = misc_funcs.get_filenames('.edf', param_path)
        path_to_edf = param_path
    else:
        prefix = param_list[0][1]
        list_edfs = [x[-1].split("/")[-1] for x in param_list]
        path_to_edf = "/".join(param_list[0][-1].split("/")[:-1]) + "/"
    # If all is recovered from a single path then list_edfs is obtained
    # via misc_funcs.get_filenames, otherwise it needs to be done from the
    # parameter list passed to the function
    name = prefix + ".pdf"
    num_pages = (len(list_edfs) // frames_per_page) + 1
    list_lists = np.array_split(list_edfs, num_pages)
    # Set context manager for PdfPages
    with PdfPages(name) as pdf:
        for list_files in list_lists:
            plotting_funcs.plot_multi_2d(list(list_files), path_to_edf,
                                         par_nrows, par_ncols, par_dist,
                                         title_override=title_override)
            if title_override:
                plt.suptitle(prefix, va='top')
            pdf.savefig(dpi=par_dpi)
            plt.close()
