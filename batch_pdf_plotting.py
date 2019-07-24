#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:53:45 2018

@author: juanenciso
"""

import sys
from multi_page_plotting import multiple_pages

# %%
if __name__ == "__main__":
    if not sys.stdin.isatty():
        list_specs = [x.split() for x in sys.stdin.readlines()]
        multiple_pages(param_list=list_specs, par_dist=300)
    else:
        path_to_edf = sys.argv[1]
        multiple_pages(param_path=path_to_edf, par_dist=300)
