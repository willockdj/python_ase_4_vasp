# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 22:29:42 2019

@author: Dave Willock
"""
import csv
import sys
import os

from physical_constants import *
from vib_energy_frm_list import vib_energy

# Will use numpy for arrays
"""Some information """

import numpy as np

# And plot it
import matplotlib.pyplot as plt

from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.thermochemistry import HarmonicThermo
from ase.visualize import view
from ase.io import read, write, Trajectory
from ase.io.vasp import write_vasp

from ase import Atom
from ase import Atoms


physical_constants()
#
##def read_car(fileobj, index):
def read_outcar_freq(fileobj):
    from ase import Atoms
    lines = fileobj.readlines()
    nlines = len(lines)  
#    print('File contains:',nlines,' lines')
#
#     Initiate arrays for collecting information from gaussian output
#
    symbols = []
    positions = []
    masses = []
    type_masses=[]
    type_labels=[]
    num_ions_per_type=[]
    freq= []
    alatt=[]
    blatt=[]
    clatt=[]
    latt_vecs=[]
    nfreq=0
    have_masses=False
    have_types=False
    have_num_types=False
    have_labelled=False
    have_coords=False
    have_latt_vecs=False
    read_coords=False
    read_latt_vecs=False
##
## This loop will now run through all the lines of the file so need to
## extract information from the lines as appropriate, hence the if blocks
##
    have_coords = False
    coord_line = -10
    for iline in range(len(lines)):
        line=lines[iline]
        linesplit = line.split()
        nwords = len(linesplit)
#
# Read in the lattice vectors
#
        if read_latt_vecs and nwords < 3:
            read_latt_vecs=False

        if read_latt_vecs:
            latt_vecs.append(float(linesplit[0]))
            latt_vecs.append(float(linesplit[1]))
            latt_vecs.append(float(linesplit[2]))
            have_latt_vecs=True
            print(latt_vecs)

        if not have_latt_vecs and nwords > 5 and 'direct' in linesplit[0] and 'vectors' in linesplit[2]:
            print("Found Lattice vector lines")
            read_latt_vecs=True
        
#
# Get masses from POMAS line which lists them
#
        if nwords > 0 and 'POMAS' in linesplit[0] and 'ZVAL' not in linesplit[3]:
            print(line)
            num_species=nwords-2
            for ispec in range(0,num_species):
                type_masses.append(float(linesplit[ispec+2]))
            print(type_masses)
            have_masses=True
#
# Get species labels from TITEL lines
#
        if nwords > 0 and 'TITEL' in linesplit[0]:
            type_labels.append(linesplit[3])
            print(type_labels)
            have_types=True
#
# Get number of ions per type from :    ions per type =  
#
        if nwords > 0 and 'ions' in linesplit[0] and 'type' in linesplit[2]:
            print(line)
            num_species=nwords-4
            for ispec in range(0,num_species):
                num_ions_per_type.append(int(linesplit[ispec+4]))
            print(num_ions_per_type)
            have_num_types=True
            
        if have_masses and have_types and have_num_types and not have_labelled:
            for ispec in range(0,num_species):
                for iatom in range(0,num_ions_per_type[ispec]):
                    symbols.append(type_labels[ispec])
                    masses.append(type_masses[ispec])
            have_labelled=True
            print(symbols)
            print(masses)

#
# Coordinates start with "position of ions in cartesian coordinates  (Angst):"
        if nwords == 7:
            if 'position' in linesplit[0] and 'cartesian' in linesplit[4]:
                print('have coords line')
                coord_line=iline
                read_coords = True
                natoms = 0
                positions = []
                
        if read_coords and nwords == 3:
#            print('atom line:')
#            print(line)
            x = linesplit[0]
            y = linesplit[1]
            z = linesplit[2]
            have_coords=True
                            
            positions.append([float(x), float(y), float(z)])
            natoms= natoms+1
#
# Recognise end of atom co-ordinates as the blank line under the list and close off
# reading the coordinates.
#                   
        if read_coords and nwords == 0:
            read_coords=False
            print(positions)
#
# Get the calculated frequencies
# Lines like: 60 f  =    0.338796 THz     2.128720 2PiTHz   11.301027 cm-1
#
        if nwords > 7 and 'f' in linesplit[1] and 'THz' in linesplit[4] and 'cm-1' in linesplit[8]:
            print('Found freq line:',line)
            freq.append(float(linesplit[7]))
            print('Read freq: ', nfreq,'  ',freq[nfreq])
            nfreq=nfreq+1
# and imaginary modes reported like:
# 7 f/i=    0.455872 THz     2.864330 2PiTHz   15.206261 cm-1     
        if nwords > 7 and 'f/i' in linesplit[1] and 'THz' in linesplit[3] and 'cm-1' in linesplit[7]:
            print('Found freq line:',line)
            print("That one is imaginary")
            freq.append(-float(linesplit[6]))
            print('Read freq: ', nfreq,'  ',freq[nfreq])
            nfreq=nfreq+1
#
#    print(symbols)
#    print(positions)           
#    print(masses)

#    for iatom in range(0,natoms):
#            atoms.append(Atom(symbols[iatom],positions[iatom]))
#    atoms.set_masses(masses)
    bail_out=False
    if not have_latt_vecs:
        print("ERROR: OUTCAR file does not contain standard direct lattice vectors information")
        bail_out=True
    if not have_coords:
        print("ERROR: OUTCAR file does not contain atomic positions")
        bail_out=True
        
    if bail_out:
        print("DEBUG: Bailing out....because of above errors.")
        sys.exit()
        
    atoms = Atoms(symbols,positions)
#
# set ASE cell information
#
    alatt=(latt_vecs[0],latt_vecs[1],latt_vecs[2])
    blatt=(latt_vecs[3],latt_vecs[4],latt_vecs[5])
    clatt=(latt_vecs[6],latt_vecs[7],latt_vecs[8])
    atoms.set_cell((alatt,blatt,clatt))
#
# Take a look at the system
#
    view(atoms)
#
#
    return (atoms,natoms,freq, nfreq)
#
# Function to check frequency list for negative modes
#
def check_modes(freq_list, num):
    found = False
    new_list=[]
    num_new = 0
    num_imag = 0
    for ifreq in range(0,num):
        if freq_list[ifreq] < 0:
            print("Warning negative mode found in freq list, removing")
            found = True
            num_imag = num_imag + 1
        else:
            new_list.append(freq_list[ifreq])
            num_new=num_new+1
    
    return (new_list, num_new, num_imag, found)
#
# Read a gaussian output file
#
current_dir = os.getcwd()
#
# Set the files containing the frequency jobs, 
# num_react_files should be how many reactant files there are
# num_prod_files should be how many product files there are
#
num_react_files=2
num_prod_files=1
react_files=[]
prod_files=[]
react_stoichiometry=[]
prod_stoichiometry=[]
react_vaspE=[]
prod_vaspE=[]
#
# Directory for OUTCAR files and names of the OUTCAR files to use
# Note that the second is the molecule in gas phase
# stoichiometry defines the number of each species in the reaction you are
# trying to simulate
# vaspE is the total electronic energy from the VASP optimisation calculation
#
# Note: Always put the molecular adsorbate second. This should be the only
# Note: term that can have a stoichiometry not equal to 1.
#
# The clean surface
data_dir = current_dir + "\\data_files\\"
react_files.append("OUTCAR_111_freq")
react_stoichiometry.append(1)
react_vaspE.append(-1182.44057011)
#
# The molecular reference
react_files.append("OUTCAR_h2o_freq")
react_stoichiometry.append(16)
react_vaspE.append(-14.23575147)
#
# The combined system.
prod_files.append("OUTCAR_1ML_111_freq")
prod_stoichiometry.append(1)
prod_vaspE.append(-1422.23760559)
#
# define file for output
#
output_file = current_dir + "\\thermo_analysie_CeO2_111_1ML.out"
plot_file = current_dir + "\\thermo_plot_CeO2_111_1ML.csv"
#
# Define pressure for calculations, in Pa
# Temperature range in K
# pressure set to 0.6 mbar from Louise's estimate
pressure=101325.*0.0006
temp_range=(10,1000)
temp_step=10
#
out_file_pntr=open(output_file,'w+')
plot_file_pntr=open(plot_file,'w+')
out_file_pntr.write("Data files taken from directory: %s\n\n" % data_dir)
#
# input and output units for the vibrational energy calculation
#
units = [ "cm-1", "eV" ]
#
#
# Open the files and read the contents, reactants first
#
react=[]
num_react_atoms=[]
react_freq=[]
num_react_freq=[]
react_struct=[]
react_vib_energies=[]
react_ZPE=[]
react_H=[]
react_S=[]
react_G=[]
react_geom=[]
#
# Loop over reactant files
#
for ifile in range(0,num_react_files):
    file_pntr=open(data_dir+react_files[ifile],'rt')
    print("Reading from file: ", react_files[ifile])

    atoms, natoms, freq, nfreq = read_outcar_freq(file_pntr)
    react.append(atoms)
    num_react_atoms.append(natoms)
    react_freq.append(freq)
    print
    print('Back in main....ifile = ', ifile,' has ',num_react_atoms[ifile],' atoms.')
    print('Back in main....ifile = ', ifile,' has  ', nfreq, ' frequencies.')
#
# Write information to the summary file
#
    out_file_pntr.write("\n\nReactant number: %d from file %s\n" % ((ifile+1),react_files[ifile]))
    out_file_pntr.write("number of atoms: %d\n" % (num_react_atoms[ifile]))
#
    for iatom in range(0,num_react_atoms[ifile]):
        out_file_pntr.write("%s %f\n" % (react[ifile][iatom].symbol, react[ifile][iatom].mass))
#
# Note that linear is curremtly only detected based on atom count: diatomics are linear
# everything else is non-linear
#
    if num_react_atoms[ifile] == 2:
        react_geom.append('linear')
    else:
        react_geom.append('nonlinear')

#
# Second file for react is the water molecule
#
    if ifile == 1:
        react_struct.append('molecule')
        #
        out_file_pntr.write("This reactant molecule is %s\n\n" % react_geom[ifile])
    else:
        react_struct.append('slab')
        out_file_pntr.write("This reactant is a slab and will be treated in the harmonic approximation.\n\n")
# Write out the reactant slab
        write_vasp('slab_react_111_clean_check.vasp', atoms ,vasp5=True)
        

#
# For molecules that are non-linear remove last six frequencies that vasp will have calculated but are really translations
# and rotations
    if react_struct[ifile] == 'molecule':
        if react_geom[ifile] == 'nonlinear':
            out_file_pntr.write("For molecule species that is non-linear removing last six frequencies..nfreq was %d..\n" % nfreq)
            nfreq=nfreq-6
            out_file_pntr.write("Now.nfreq is %d..\n" % nfreq)
#
    num_react_freq.append(nfreq)
#
# Check for negative modes before going further
#
    react_freq[ifile], num_react_freq[ifile], num_imag, found = check_modes(react_freq[ifile], num_react_freq[ifile])
#
    if found:
        print("Note negative modes have been removed from frequencies list")
        out_file_pntr.write("\n\nNote %d negative modes have been removed from frequencies list\n\n" 
                                                                                            % num_imag)
#
# Work out vibrational energy states
#
    react_vib_energies.append(vib_energy(react_freq[ifile], units))
 
    for ifreq in range(0,num_react_freq[ifile]):
        out_file_pntr.write("mode: %d freq: %10.6f cm-1 energy: %10.6f %s\n" 
                            % ((ifreq+1),react_freq[ifile][ifreq], react_vib_energies[ifile][ifreq], units[1]))

    react_ZPE.append(0.5*np.sum(react_vib_energies[ifile]))
    out_file_pntr.write("\nTotal ZPE: %10.6f %s\n" % (react_ZPE[ifile],units[1]))
#
# carry out thermocalculations
#
    out_file_pntr.write("\nReactant absolute thermodynamic quantities at T = 298.15 K and P = 1 atm.\n")
    
    if react_struct[ifile] == 'molecule':
        thermo = IdealGasThermo(vib_energies=react_vib_energies[ifile],
                                atoms=react[ifile],
                                geometry=react_geom[ifile],
                                symmetrynumber=2, spin=0)
#
        react_H.append(thermo.get_enthalpy(temperature=298.15))
        react_S.append(thermo.get_entropy(temperature=298.15, pressure=pressure))
        react_G.append(thermo.get_gibbs_energy(temperature=298.15, pressure=pressure))
#
        out_file_pntr.write("\nVASP electronic energy, E : %10.6f %s\n" % (react_vaspE[ifile],units[1]))
        out_file_pntr.write("\nVibrational contribution to U : %10.6f %s\n" % (react_H[ifile],units[1]))
        react_H[ifile]=react_H[ifile]+react_vaspE[ifile]
        react_G[ifile]=react_G[ifile]+react_vaspE[ifile]

        out_file_pntr.write("\nTotal Enthalpy, H: %10.6f %s\n" % (react_H[ifile],units[1]))
        out_file_pntr.write("Entropy,  S: %10.6f %s K-1\n" % (react_S[ifile],units[1]))
        out_file_pntr.write("Gibbs FE, G: %10.6f %s\n" % (react_G[ifile],units[1]))
#

    else:
        thermo = HarmonicThermo(vib_energies=react_vib_energies[ifile],
                                potentialenergy=0)
#
        react_H.append(thermo.get_internal_energy(temperature=298.15))
        react_S.append(thermo.get_entropy(temperature=298.15))
        react_G.append(thermo.get_helmholtz_energy(temperature=298.15))
#
        out_file_pntr.write("\nVASP electronic energy, E : %10.6f %s\n" % (react_vaspE[ifile],units[1]))
        out_file_pntr.write("\nVibrational contribution to U : %10.6f %s\n" % (react_H[ifile],units[1]))
        react_H[ifile]=react_H[ifile]+react_vaspE[ifile]
        react_G[ifile]=react_G[ifile]+react_vaspE[ifile]
        out_file_pntr.write("\nTotal Enthalpy, H (strictly internal energy U for slab) : %10.6f %s\n" % (react_H[ifile],units[1]))
        out_file_pntr.write("Entropy,  S (defined using just the harmonic vibs): %10.6f %s K-1\n" % (react_S[ifile],units[1]))
        out_file_pntr.write("Gibbs FE, G (strictly Helmholtz for slab)         : %10.6f %s\n" % (react_G[ifile],units[1]))
#
# Now read products
#
prod=[]
num_prod_atoms=[]
prod_freq=[]
num_prod_freq=[]
prod_vib_energies=[]
prod_struct=[]
prod_geom=[]
prod_ZPE=[]
prod_H=[]
prod_S=[]
prod_G=[]
for ifile in range(0,num_prod_files):
    file_pntr=open(data_dir+prod_files[ifile],'rt')

    atoms, natoms, freq, nfreq = read_outcar_freq(file_pntr)
    prod.append(atoms)
    num_prod_atoms.append(natoms)
    prod_freq.append(freq)
    num_prod_freq.append(nfreq)
    print
    print('Back in main....ifile = ', ifile,' has ',num_prod_atoms[ifile],' atoms.')
    print('Back in main....ifile = ', ifile,' has  ', num_prod_freq[ifile], ' frequencies.')
#
# Write information to the summary file
#
    out_file_pntr.write("\n\nProduct number: %d from file %s\n" % ((ifile+1),prod_files[ifile]))
    out_file_pntr.write("number of atoms: %d\n" % (num_prod_atoms[ifile]))
#
    for iatom in range(0,num_prod_atoms[ifile]):
        out_file_pntr.write("%s %f\n" % (prod[ifile][iatom].symbol, prod[ifile][iatom].mass))
#
# Note that linear is curremtly only detected based on atom count: diatomics are linear
# everything else is non-linear
#
    if num_prod_atoms[ifile] == 2:
        prod_geom.append('linear')
    else:
        prod_geom.append('nonlinear')
#
    out_file_pntr.write("This molecule is %s\n\n" % prod_geom[ifile])
#
# All product structures are slabs in this case
#
    prod_struct.append('slab')
#
# Check for negative modes before going further
#
    prod_freq[ifile], num_prod_freq[ifile], num_imag, found = check_modes(prod_freq[ifile], num_prod_freq[ifile])
#
    if found:
        print("Note negative modes have been removed from frequencies list")
        out_file_pntr.write("\n\nNote %d negative modes have been removed from frequencies list\n\n" % num_imag)
#
    prod_vib_energies.append(vib_energy(prod_freq[ifile], units))
 
    for ifreq in range(0,num_prod_freq[ifile]):
        out_file_pntr.write("mode: %d freq: %10.6f cm-1 energy: %10.6f %s\n" 
                            % ((ifreq+1),prod_freq[ifile][ifreq], prod_vib_energies[ifile][ifreq], units[1]))

    prod_ZPE.append(0.5*np.sum(prod_vib_energies[ifile]))
    out_file_pntr.write("\nTotal ZPE: %10.6f %s\n" % (prod_ZPE[ifile],units[1]))
#
# carry out thermocalculations
#
    out_file_pntr.write("\nProduct absolute thermodynamic quantities at T = 298.15 K and P = 1 atm.\n")
#
    if prod_struct[ifile] == 'molecule':
        thermo = IdealGasThermo(vib_energies=prod_vib_energies[ifile],
                                atoms=prod[ifile],
                                geometry=prod_geom[ifile],
                                symmetrynumber=2, spin=0)
#
# Get standard state items
#
        prod_H.append(thermo.get_enthalpy(temperature=298.15))
        prod_S.append(thermo.get_entropy(temperature=298.15, pressure=pressure))
        prod_G.append(thermo.get_gibbs_energy(temperature=298.15, pressure=pressure))
#
        out_file_pntr.write("\nVASP electronic energy, E : %10.6f %s\n" % (prod_vaspE[ifile],units[1]))
        out_file_pntr.write("\nVibrational contribution to U : %10.6f %s\n" % (prod_H[ifile],units[1]))
        prod_H[ifile]=prod_H[ifile]+prod_vaspE[ifile]
        prod_G[ifile]=prod_G[ifile]+prod_vaspE[ifile]
        out_file_pntr.write("\nTotal Enthalpy, H: %10.6f %s\n" % (prod_H[ifile],units[1]))
        out_file_pntr.write("Entropy,  S: %10.6f %s K-1\n" % (prod_S[ifile],units[1]))
        out_file_pntr.write("Gibbs FE, G: %10.6f %s\n" % (prod_G[ifile],units[1]))
#
    else:
        thermo = HarmonicThermo(vib_energies=prod_vib_energies[ifile],
                                potentialenergy=0)
#
# Get standard state items
#
        prod_H.append(thermo.get_internal_energy(temperature=298.15))
        prod_S.append(thermo.get_entropy(temperature=298.15))
        prod_G.append(thermo.get_helmholtz_energy(temperature=298.15))
#
        out_file_pntr.write("\nVASP electronic energy, E : %10.6f %s\n" % (prod_vaspE[ifile],units[1]))
        out_file_pntr.write("\nVibrational contribution to U : %10.6f %s\n" % (prod_H[ifile],units[1]))
        prod_H[ifile]=prod_H[ifile]+prod_vaspE[ifile]
        prod_G[ifile]=prod_G[ifile]+prod_vaspE[ifile]
        out_file_pntr.write("\nTotal Enthalpy, H (strictly internal energy U for slab) : %10.6f %s\n" % (prod_H[ifile],units[1]))
        out_file_pntr.write("Entropy,  S (defined using just the harmonic vibs): %10.6f %s K-1\n" % (prod_S[ifile],units[1]))
        out_file_pntr.write("Gibbs FE, G (strictly Helmholtz for slab)         : %10.6f %s\n" % (prod_G[ifile],units[1]))
#
#
# Work out energy changes for the reaction
#
#
out_file_pntr.write("\nChanges for this reaction under standard conditions:\n\n")

tot_ZPE=0
tot_vaspE=0
tot_H=0
tot_S=0
tot_G=0
tot_dof=0
for ifile in range(0,num_react_files):
    out_file_pntr.write("Reactant %d read from %s\n" % (ifile,react_files[ifile]))
    if react_struct[ifile] == 'molecule':
        dof=num_react_freq[ifile]+6
        out_file_pntr.write("....is a molecule with\n%d vibrations, 3 rotations and 3 translations, i.e. %d dof.\n" %
                                            (num_react_freq[ifile], dof))
        if react_stoichiometry[ifile] > 1:
            out_file_pntr.write("There are %d of these in the reaction so %d dof in all\n" %
                                (react_stoichiometry[ifile],react_stoichiometry[ifile]*dof))                     
    else:
        dof=num_react_freq[ifile]
        out_file_pntr.write("....is a slab,\ntreated in the harmonic approx., only vibs, i.e. %d dof.\n" % 
                                            dof )
        if react_stoichiometry[ifile] > 1:
            out_file_pntr.write("There are %d of these in the reaction so %d dof in all\n" %
                                (react_stoichiometry[ifile],react_stoichiometry[ifile]*dof))                     

    tot_dof = tot_dof+react_stoichiometry[ifile]*dof        
    tot_ZPE=tot_ZPE-react_stoichiometry[ifile]*react_ZPE[ifile]
    tot_vaspE=tot_vaspE-react_stoichiometry[ifile]*react_vaspE[ifile]
    tot_H=tot_H-react_stoichiometry[ifile]*react_H[ifile]
    tot_S=tot_S-react_stoichiometry[ifile]*react_S[ifile]
    tot_G=tot_G-react_stoichiometry[ifile]*react_G[ifile]

for ifile in range(0,num_prod_files):
    out_file_pntr.write("Product %d read from %s\n" % (ifile,prod_files[ifile]))
    if prod_struct[ifile] == 'molecule':
        dof=num_prod_freq[ifile]+6
        out_file_pntr.write("....is a molecule with\n %d vibrations, 3 rotations and 3 translations, i.e. %d dof.\n" %
                                            (num_prod_freq[ifile], dof))
        if prod_stoichiometry[ifile] > 1:
            out_file_pntr.write("There are %d of these in the reaction so %d dof in all\n" %
                                (prod_stoichiometry[ifile],prod_stoichiometry[ifile]*dof))                     

    else:
# Write out the product slab
        write_vasp('slab_prod_111_H2O_check.vasp', atoms ,vasp5=True)
        dof=num_prod_freq[ifile]
        out_file_pntr.write("....is a slab,\ntreated in the harmonic approx., only vibs, i.e. %d dof.\n" % 
                                            dof )
        if prod_stoichiometry[ifile] > 1:
            out_file_pntr.write("There are %d of these in the reaction so %d dof in all\n" %
                                (prod_stoichiometry[ifile],prod_stoichiometry[ifile]*dof))                     

    tot_dof = tot_dof-prod_stoichiometry[ifile]*dof        

    tot_ZPE=tot_ZPE+prod_stoichiometry[ifile]*prod_ZPE[ifile]
    tot_vaspE=tot_vaspE+prod_stoichiometry[ifile]*prod_vaspE[ifile]
    tot_H=tot_H+prod_stoichiometry[ifile]*prod_H[ifile]
    tot_S=tot_S+prod_stoichiometry[ifile]*prod_S[ifile]
    tot_G=tot_G+prod_stoichiometry[ifile]*prod_G[ifile]
#
#
#
if tot_dof != 0:
    out_file_pntr.write("WARNING: degrees of freedom before and after reaction\n")
    out_file_pntr.write("WARNING: DO NOT BALANCE, tot_dof=%d\n" % tot_dof)
    out_file_pntr.write("WARNING: If this is what you want you should know why!!!!\n")
    out_file_pntr.write("WARNING: Otherwise check that there are the\n")
    out_file_pntr.write("WARNING: same number of free atoms in surface slabs\n")
    out_file_pntr.write("WARNING: and make sure no imaginary modes \n")
    out_file_pntr.write("WARNING: have had to be left out\n")

out_file_pntr.write("\n\nChanges in energy per set of molecules\n")
out_file_pntr.write("\ndelta ZPE: %10.6f %s\n" % (tot_ZPE,units[1]))
out_file_pntr.write("\ndelta E from VASP: %10.6f %s\n" % (tot_vaspE,units[1]))
out_file_pntr.write("delta Enthalpy, DH: %10.6f %s\n" % (tot_H,units[1]))
out_file_pntr.write("delta Entropy,  DS: %10.6f %s K-1\n" % (tot_S,units[1]))
out_file_pntr.write("delta Gibbs FE, DG: %10.6f %s\n" % (tot_G,units[1]))
#
out_file_pntr.write("\n\nChanges in energy per adsorbate molecule\n")
out_file_pntr.write("\ndelta ZPE: %10.6f %s\n" % (tot_ZPE/react_stoichiometry[1],units[1]))
out_file_pntr.write("\ndelta E from VASP: %10.6f %s\n" % (tot_vaspE/react_stoichiometry[1],units[1]))
out_file_pntr.write("delta Enthalpy, DH: %10.6f %s\n" % (tot_H/react_stoichiometry[1],units[1]))
out_file_pntr.write("delta Entropy,  DS: %10.6f %s K-1\n" % (tot_S/react_stoichiometry[1],units[1]))
out_file_pntr.write("delta Gibbs FE, DG: %10.6f %s\n" % (tot_G/react_stoichiometry[1],units[1]))
#
# Set up a csv file that generates the plotting data
#
plot_file_pntr.write("Temperature (K), delta H (eV), -T delta_S (eV), delta_G (eV)\n")
# 
print("temp min",temp_range[0]," temp_step: ", temp_step, " temp_max: ", temp_range[1])
for temp in range(temp_range[0],temp_range[1],temp_step):
    print("Temperature: ",temp)
    tot_ZPE=0
    tot_H=0
    tot_S=0
    tot_G=0
#
    react_H=[]
    react_S=[]
    react_G=[]
#
    prod_H=[]
    prod_S=[]
    prod_G=[]
#
# carry out product thermocalculations
#

    for ifile in range(0,num_prod_files):
#
        if prod_struct[ifile] == 'molecule':
            thermo = IdealGasThermo(vib_energies=prod_vib_energies[ifile],
                                    atoms=prod[ifile],
                                    geometry=prod_geom[ifile],
                                    symmetrynumber=2, spin=0)
#
# Get standard state items
#
            prod_H.append(thermo.get_enthalpy(temperature=temp))
            prod_S.append(thermo.get_entropy(temperature=temp, pressure=pressure))
            prod_G.append(thermo.get_gibbs_energy(temperature=temp, pressure=pressure))
#
        else:
            thermo = HarmonicThermo(vib_energies=prod_vib_energies[ifile],
                                    potentialenergy=0)
#
# Get standard state items
#
            prod_H.append(thermo.get_internal_energy(temperature=temp))
            prod_S.append(thermo.get_entropy(temperature=temp))
            prod_G.append(thermo.get_helmholtz_energy(temperature=temp))
#
# Add in electronic contribution
#
        prod_H[ifile]=prod_H[ifile]+prod_vaspE[ifile]
        prod_G[ifile]=prod_G[ifile]+prod_vaspE[ifile]
#
        tot_H=tot_H+prod_stoichiometry[ifile]*prod_H[ifile]
        tot_S=tot_S+prod_stoichiometry[ifile]*prod_S[ifile]
        tot_G=tot_G+prod_stoichiometry[ifile]*prod_G[ifile]
#
# carry out reactant thermocalculations
#
    for ifile in range(0,num_react_files):
#
        if react_struct[ifile] == 'molecule':
            thermo = IdealGasThermo(vib_energies=react_vib_energies[ifile],
                                    atoms=react[ifile],
                                    geometry=react_geom[ifile],
                                    symmetrynumber=2, spin=0)
#
# Get standard state items
#
            react_H.append(thermo.get_enthalpy(temperature=temp))
            react_S.append(thermo.get_entropy(temperature=temp, pressure=pressure))
            react_G.append(thermo.get_gibbs_energy(temperature=temp, pressure=pressure))
#
        else:
            thermo = HarmonicThermo(vib_energies=react_vib_energies[ifile],
                                    potentialenergy=0)
#
# Get standard state items
#
            react_H.append(thermo.get_internal_energy(temperature=temp))
            react_S.append(thermo.get_entropy(temperature=temp))
            react_G.append(thermo.get_helmholtz_energy(temperature=temp))
#
# Add in electronic contribution
#
        react_H[ifile]=react_H[ifile]+react_vaspE[ifile]
        react_G[ifile]=react_G[ifile]+react_vaspE[ifile]
#
        tot_H=tot_H-react_stoichiometry[ifile]*react_H[ifile]
        tot_S=tot_S-react_stoichiometry[ifile]*react_S[ifile]
        tot_G=tot_G-react_stoichiometry[ifile]*react_G[ifile]
       
    plot_file_pntr.write("%10.6f, %10.6f, %10.6f, %10.6f\n" 
                            % (temp, tot_H/react_stoichiometry[1], -temp*tot_S/react_stoichiometry[1], tot_G/react_stoichiometry[1]))
    
#
#
#
plot_file_pntr.close()
out_file_pntr.close()
