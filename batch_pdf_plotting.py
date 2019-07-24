#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:53:45 2018

@author: juanenciso
"""

# %% Append current directory to path so modules can be imported
import sys

# builtin module imports
# local imports
from multi_page_plotting import multiple_pages

# %%
if __name__ == "__main__":
    if not sys.stdin.isatty():
        # There's data in stdin, should procceed using a list of files
        list_specs = [x.split() for x in sys.stdin.readlines()]
        multiple_pages(param_list=list_specs, par_dist=300)
    else:
        # There is not data in stdin, use first item in sys.argv as path
        # containing the edf files
        path_to_edf = sys.argv[1]
        multiple_pages(param_path=path_to_edf, par_dist=300)


# SOMEDAY: Extend this so more parameters of the plot can be modified when
# calling the script
