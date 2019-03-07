# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 14:27:54 2018

@author: dkorff
"""

import numpy as np
import cantera as ct

from li_ion_battery_p2d_inputs import Inputs

# Import Cantera objects:
anode_obj = ct.Solution(Inputs.ctifile,Inputs.anode_phase)
elyte_obj = ct.Solution(Inputs.ctifile,Inputs.elyte_phase)
cathode_obj = ct.Solution(Inputs.ctifile,Inputs.cathode_phase)
conductor_obj = ct.Solution(Inputs.ctifile,Inputs.metal_phase)

anode_surf_obj = ct.Interface(Inputs.ctifile,Inputs.anode_surf_phase,
    [anode_obj,elyte_obj,conductor_obj])
cathode_surf_obj = ct.Interface(Inputs.ctifile,Inputs.cathode_surf_phase,
    [cathode_obj,elyte_obj,conductor_obj])


"""========================================================================="""
"""========================================================================="""
"""========================================================================="""

class anode():

    # Set flag so solver knows whether to implement the anode:
    flag = Inputs.flag_anode

    # Number of nodes in the y-direction:
    npoints = Inputs.npoints_anode

    # Number of "shells" in anode particle:
    nshells = Inputs.nshells_anode

    # Initial conditions: UPDATE TO GENERALIZE FOR MULT. SPEC - DK 8/31/18
    LiC6_0 = Inputs.SOC_0; C6_0 = 1 - LiC6_0
    X_an_init = '{}:{}, {}:{}'.format(Inputs.Li_species_anode, LiC6_0,
        Inputs.Vac_species_anode, C6_0)

    # Anode variables for a given volume include X_Li for each shell, Phi_an,
    #       Phi_elyte, and rho_k_elyte for all elyte species.
    nVars = nshells + 2 + elyte_obj.n_species

    # Pointers
    ptr = {}
    ptr['iFar'] = 1
    ptr['X_ed'] = np.arange(0, Inputs.nshells_anode)
    ptr['rho_k_elyte'] = nshells + np.arange(0,elyte_obj.n_species)
    ptr['Phi_ed'] = nshells + elyte_obj.n_species
    ptr['Phi_dl'] = nshells + elyte_obj.n_species + 1

    # Anode/elyte interface area per unit volume
    #   [m^2 interface / m_3 total electrode volume]
    # For spherical particles, the total surface area per unit volume can be
    #   calculated from the geometry. Since some particles will overlap, we
    #   multiply by an 'overlap' input parameter.
    A_surf = (1-Inputs.overlap_an)*6*Inputs.eps_solid_an/Inputs.d_part_an

    # Create Cantera objects
    anode_obj.X = X_an_init
    anode_obj.TP = Inputs.T, ct.one_atm
    elyte_obj.TP = Inputs.T, ct.one_atm
    anode_surf_obj.TP = Inputs.T, ct.one_atm

    X_an_0 = anode_obj.X

    if hasattr(Inputs, 'X_elyte_init'):
        elyte_obj.X = X_elyte_init

    X_elyte_0 = elyte_obj.X

    # Set up solution vector
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])

    # Array of offsets to point to each node's variable:
    offsets = np.arange(0,int(nSV),int(nVars))
    for j in range(npoints):
        SV_0[offsets[j] + ptr['X_ed']] = np.ones([nshells])*X_an_0[0]
        SV_0[offsets[j] + ptr['rho_k_elyte']] = elyte_obj.Y*elyte_obj.density_mass
        SV_0[offsets[j] + ptr['Phi_ed']] = Inputs.Phi_anode_init
        SV_0[offsets[j] + ptr['Phi_dl']] = Inputs.Phi_elyte_init - Inputs.Phi_anode_init
        # Electrolyte electric potential = 0, so V_dl_init = V_an_init

    # Calculate the current density [A/m^2] corresponding to a C_rate of 1:
    oneC = geom['eps_ed']*anode_obj.density_mole*Inputs.H_an*ct.faraday/3600

    # Calculate the actual current density:
    # The minus sign is because we begin with the charging reaction, which
    #   delivers negative charge to the anode:
    i_ext = -Inputs.C_rate*oneC

    # Store parameters as class object attributes:
    T = Inputs.T
    C_dl = Inputs.C_dl_an
    X_Li_max = Inputs.SOC_max
    X_Li_min = Inputs.SOC_min
    D_Li_ed = Inputs.D_Li_an

    # Geometric parameters:
    eps_ed = Inputs.eps_solid_an
    eps_elyte = 1 - Inputs.eps_solid_an
    tau_ed = Inputs.tau_an
    r_pore = Inputs.r_p_an
    d_part = Inputs.d_part_an
    dyInv = npoints/Inputs.H_an
    dr = geom['d_part']*0.5/nshells

    # Calculate the percent volume of a single graphite particle that exists in
    #   each 'shell'. I.e. for shell j, what is the volume of that shell,
    #   relative to the total particle volume? The radius of the volume is
    #   currently discretized evenly (i.e. 'dr' the differential radius is
    #   constant). Certainly other discretizations (such as constant
    #   differential volume) are possible (and perhaps better).
    #
    #   Because the volume is 4/3 pi*r^3, the volume of the shell relative to
    #   the total volume is (r_shell/r_particle)^3, and the differential volume
    #   relative to the total, for shell 'j' is:
    #       (r_shell(j+1)^3 - r_shell(j)^3)/r_particle^3
    #   Because the radius is discretized evenly, the radius of shell j, r_j,
    #   relative to the total radius r_particle, is:
    #       r_j/r_particle = j/nshells

    params['V_shell'] = np.zeros([nshells])
    nsr3 = 1/nshells/nshells/nshells
    for j in np.arange(0, nshells, 1):
        params['V_shell'][j] = ((j+1)**3 - (j)**3)*nsr3

    # Electronic conductivity of the electrode phase:
    params['sigma_eff_ed'] = Inputs.sigma_an*geom['eps_ed']/geom['tau_ed']**3
    # Species mobilities of the electrolyte phase.  Converted from user input
    #   diffusion coefficients:
    params['u_Li_elyte'] = (Inputs.D_Li_elyte*geom['eps_elyte']/ct.gas_constant
          /Inputs.T/geom['tau_ed']**3)

    # Set up algebraic variable vector:
    algvar = np.zeros([nSV])
    for j in np.arange(0, npoints):
        # X_Li_an has a differential equation:
        algvar[offsets[j] + ptr['X_ed']] = 1
        # X_Li_elyte has a differential equation:
        algvar[offsets[j] + ptr['rho_k_elyte']] = 1
        # Double Layer elec potential has a differential equation involving
        #   phi_an and phi_elyte:
        algvar[offsets[j] + ptr['Phi_dl']] = 1

        # The anode electric potential is governed by an algebraic equation.

    t_flag = []


    """========================================================================="""
    """========================================================================="""
    """========================================================================="""

class separator():
    # Set a flag to let the solver know whether to implement this class:
    flag = Inputs.flag_sep

    # Number of nodes in the y-direction
    npoints = Inputs.npoints_elyte

    # Number of variables per node:
    nVars = 1 + elyte_obj.n_species

    H = Inputs.H_elyte  # Separator thickness [m]

    tau_sep = Inputs.tau_sep  # Tortuosity of separator

    # Geometric parameters:
    eps_elyte = Inputs.eps_elyte_sep
    dyInv = npoints/H
    tau_sep = tau_sep

    ptr = {}
    ptr['rho_k_elyte'] = np.arange(0,elyte_obj.n_species)
    ptr['Phi'] = elyte_obj.n_species

    # Set up the solution vector
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])

    # Array of offsets to point to each node's variable:
    offset = int(anode.nSV)
    offsets = np.arange(offset,offset+int(nSV),int(nVars))

    for j in np.arange(0, npoints):
        SV_0[j*nVars + ptr['rho_k_elyte']] = elyte_obj.Y*elyte_obj.density_mass
        SV_0[j*nVars + ptr['Phi']] = Inputs.Phi_elyte_init

    # Set up algebraic variable vector:
    algvar = np.zeros([nSV])

    for j in range(npoints):
        algvar[j*nVars + ptr['rho_k_elyte']] = 1
        algvar[j*nVars + ptr['Phi']] = 0


"""========================================================================="""
"""========================================================================="""
"""========================================================================="""

class cathode():
    # Set flag so solver knows it's looking at the anode class
    flag = Inputs.flag_cathode

    # Number of nodes in the y-direction:
    npoints = Inputs.npoints_cathode

    # Number of "shells" in the cathode particle
    nshells = Inputs.n_shells_cathode

    # Number of state variables per node:
    nVars = nshells + 2 + elyte_obj.n_species

    # Initial conditions: UPDATE TO GENERALIZE FOR MULT. SPECIES - DK 8/31/18
    LiCoO2_0 = 1- Inputs.SOC_0; CoO2_0 = 1 - LiCoO2_0
    X_cat_init = '{}:{}, {}:{}'.format(Inputs.Li_species_cathode, LiCoO2_0,
            Inputs.Vac_species_cathode, CoO2_0)

    # Pointers
    ptr = {}
    ptr['iFar'] = 0
    ptr['X_ed'] = np.arange(0, nshells)
    ptr['rho_k_elyte'] = nshells + np.arange(0,elyte_obj.n_species)
    ptr['Phi_ed'] = nshells + elyte_obj.n_species
    ptr['Phi_dl'] = nshells + elyte_obj.n_species + 1

    # Cathode/elyte interface area per unit volume
    #   [m^2 interface / m_3 total electrode volume]
    # For spherical particles, the total surface area per unit volume can be
    #   calculated from the geometry. Since some particles will overlap, we
    #   multiply by an 'overlap' input parameter.
    A_surf = (1-Inputs.overlap_ca)*6*Inputs.eps_solid_ca/Inputs.d_part_ca

    cathode_obj.X = X_cat_init
    cathode_obj.TP = Inputs.T, ct.one_atm
    cathode_surf_obj.TP = Inputs.T, ct.one_atm

    X_cat_0 = cathode_obj.X
    X_elyte_0 = elyte_obj.X

    # Set up the solution vector
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])

    # Array of offsets to point to each node's variable:
    offset = int(anode.nSV + separator.nSV)
    offsets = np.arange(offset,offset+int(nSV),int(nVars))

    for j in range(npoints):
        SV_0[j*nVars + ptr['X_ed']] = np.ones([nshells])*X_cat_0[0]
        SV_0[j*nVars + ptr['rho_k_elyte']] = elyte_obj.Y*elyte_obj.density_mass
        SV_0[j*nVars + ptr['Phi_ed']] = Inputs.Phi_anode_init \
            + Inputs.Delta_Phi_init
        SV_0[j*nVars + ptr['Phi_dl']] = Inputs.Phi_elyte_init \
                            - (Inputs.Phi_anode_init + Inputs.Delta_Phi_init)
        # Electrolyte electric potential = 0, so V_dl_init = V_an_init


    # Calculate the current density [A/m^2] corresponding to a C_rate of 1:
    oneC = geom['eps_ed']*cathode_obj.density_mole*Inputs.H_ca*ct.faraday/3600

    # Calculate the actual current density:
    # The minus sign is because we begin with the charging reaction, which
    #   delivers negative charge to the anode:
    i_ext = -Inputs.C_rate*oneC

    #  Store Parameters in the Class object:
    T = Inputs.T
    C_dl = Inputs.C_dl_ca
    X_Li_max = Inputs.SOC_max
    X_Li_min = Inputs.SOC_min
    D_Li_ed = Inputs.D_Li_ca

    # Geometric parameters:
    eps_ed = Inputs.eps_solid_ca
    eps_elyte = 1 - Inputs.eps_solid_ca
    tau_ed = Inputs.tau_ca
    r_p = Inputs.r_p_ca
    d_part = Inputs.d_part_ca
    dyInv = npoints/Inputs.H_ca
    dr = geom['d_part']*0.5/nshells

    # Calculate the percent volume of a single graphite particle that exists in
    #   each 'shell'. I.e. for shell j, what is the volume of that shell,
    #   relative to the total particle volume? The radius of the volume is
    #   currently discretized evenly (i.e. 'dr' the differential radius is
    #   constant). Certainly other discretizations (such as constant
    #   differential volume) are possible (and perhaps better).
    #
    #   Because the volume is 4/3 pi*r^3, the volume of the shell relative to
    #   the total volume is (r_shell/r_particle)^3, and the differential volume
    #   relative to the total, for shell 'j' is:
    #       (r_shell(j+1)^3 - r_shell(j)^3)/r_particle^3
    #   Because the radius is discretized evenly, the radius of shell j, r_j,
    #   relative to the total radius r_particle, is:
    #       r_j/r_particle = j/nshells

    params['V_shell'] = np.zeros([nshells])
    nsr3 = 1/nshells/nshells/nshells
    for j in np.arange(0, nshells, 1):
        params['V_shell'][j] = ((j+1)**3 - (j)**3)*nsr3

    params['sigma_eff_ed'] = Inputs.sigma_ca*geom['eps_ed']/geom['tau_ed']**3
    params['u_Li_elyte'] = (Inputs.D_Li_elyte*geom['eps_elyte']/ct.gas_constant
          /Inputs.T/geom['tau_ed']**3)

    # Set up algebraic variable vector:
    algvar = np.zeros([nSV])
    for j in np.arange(0, npoints):
        # X_Li_an has a differential equation:
        algvar[j*nVars + ptr['X_ed']] = 1
        # X_Li_elyte has a differential equation:
        algvar[j*nVars + ptr['rho_k_elyte']] = 1
        # Double Layer elec potential has a differential equation involving
        #   phi_an and phi_elyte:
        algvar[j*nVars + ptr['Phi_dl']] = 1

        # The anode electric potential is governed by an algebraic equation.

    t_flag = []

"""========================================================================="""
"""========================================================================="""
"""========================================================================="""
