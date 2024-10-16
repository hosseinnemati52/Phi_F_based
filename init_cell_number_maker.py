#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:25:17 2024

@author: hossein
"""





import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import re
import subprocess
import os


N_runs = 10
WT_init_frac_interval = np.array([0.0, 0.0])

eps = 0.0001

if np.mean(WT_init_frac_interval) > (1.0 - eps):
    # pure WT
    
    init_num_cells = np.random.randint(10,70,N_runs)
    
    
elif np.mean(WT_init_frac_interval) < eps:
    # pure Cancer
    pass
else:
    # mixed
    pass