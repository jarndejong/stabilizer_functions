# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:45:56 2018

@author: Jarnd
"""
import numpy as np

import Stab_group_functions as stab

S = np.array([[1, 1, 1, 0], [1, 1, 0, 0],[0, 0, 1, 0]]);

print(stab.ech(S))