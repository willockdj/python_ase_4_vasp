# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 22:34:29 2019

@author: Dave
"""
from physical_constants import *

import numpy as np

def vib_energy(freq, units):
    
    if units[0] == "cm-1" :
        conv = c_light*100
        freq_SI = np.multiply(conv, freq)
    else:
        print("Unrecognised unit :", units[0], " in vib_energy_frm_list.py ")
    
    if units[1] == "eV" :
        planck_needed = h_planck / e_charge
    elif units[1] == "kJmol-1":
        planck_needed = h_planck * n_avo / 1000.0
    elif units[1] == "Ha":
        planck_needed = h_planck / (e_charge*Ha_to_eV)
    else:
        print("Unrecognised unit :", units[1], " in vib_energy_frm_list.py ")
        
 
    return np.multiply(planck_needed,freq_SI)
    