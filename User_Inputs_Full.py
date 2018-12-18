# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 14:27:54 2018

@author: dkorff
"""

import numpy as np
import cantera as ct

class Inputs():
    
    # The C-rate is the rate of charge/discharge - how many charges/discharges
    #   can be carried out in 1 hour? This sets the current density:
    C_rate = 1
    
    # Simulation temperature (or initial temperature)
    T = 300  # [K]
    
    # Faraday's constant 
    F = 96485e3  # [C/kmol]
    
    # Set initial SOC to generalize both electrod initial lithiation
    SOC_0 = 50/100.
    
"""========================================================================="""
"""========================================================================="""
"""========================================================================="""

class anode():
    # Set flag so solver knows it's looking at the anode class
    ed_flag = 1
    
    # Number of nodes in the y-direction:
    npoints = 5
    
    # Number of "shells" in anode particle:
    nshells = 5
    
    # Initial conditions: UPDATE TO GENERALIZE FOR MULT. SPEC - DK 8/31/18
    LiC6_0 = Inputs.SOC_0; C6_0 = 1 - LiC6_0
    X_an_init = 'LiC6:0.05, C6:0.95'
    X_elyte_init = 'Li(e):0.98, Solvent:0.02'  # PF6-(e):7.8730103237e-2, EC(e):2.8328131770e-1'
    V_an_init = 0
    
    # Cti file info:
    ctifile = 'LiBatteryFull.cti'
    anode_b = 'graphite'
    anode_s = 'anode_surf'
    elytephase = 'electrolyte'
    
    # Cutoff Values for lithiation and delithiation of anode:
    X_Li_an_max = 1 - 1e-2
    X_Li_an_min = 1 - X_Li_an_max
    
    # Microstructure
    phi_an = 0.6  # Graphite volume fraction [-]
    tau_an = 1.6  # Tortuosity - assume equal values for carbon and elyte [-]
    r_p = 5e-6    # Average pore radius [m]
    d_part = 5e-6 # Average particle diameter for graphite [m]
    H = 50e-6     # Anode thickness [m]
    
    # Other Parameters
    C_dl = 1.5e-2       # Double-layer capacitance [F/m^2]
    sigma_an = 75.0     # Bulk anode electrical conductivity [S/m]
    D_Li_elyte = 1e-10  # Bulk diffusion coefficient for Li+ in elyte [m^2/s]
    D_Li_an = 7.5e-16   # Bulk diffusion coefficient for Li in graphite [m^2/s]
    
    nVars = nshells + 3
    
    # Pointers
    ptr = {}
    ptr['iFar'] = 1
    ptr['X_ed'] = np.arange(0, nshells)
    ptr['X_elyte'] = nshells
    ptr['V_ed'] = nshells + 1
    ptr['V_dl'] = nshells + 2
    
    # For spherical particles, the total surface area per unit volume can be
    #   calculated from the geometry. Since some particles will overlap, we take
    #   a percentage (60%, below) of this theoretical value:
    A_surf = 0.6*6*phi_an/d_part  # Anode/elyte interface area [m^2/m_3]
    
    # Create Cantera objects
    anode = ct.Solution(ctifile, anode_b)
    anode.X = X_an_init
    anode.TP = Inputs.T, ct.one_atm
    
    X_an_0 = anode.X
    
    elyte = ct.Solution(ctifile, elytephase)
    elyte.X = X_elyte_init
    elyte.TP = Inputs.T, ct.one_atm
    
    X_elyte_0 = elyte.X
    
    an_surf = ct.Interface(ctifile, anode_s, [anode, elyte])
    an_surf.TP = Inputs.T, ct.one_atm
    
    # Set up solution vector
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])
    
    for j in range(npoints):
        offset = (j)*nVars
        SV_0[offset + ptr['X_ed']] = np.ones([nshells])*X_an_0[0]
        SV_0[offset + ptr['X_elyte']] = X_elyte_0[0]
        SV_0[offset + ptr['V_ed']] = V_an_init
        SV_0[offset + ptr['V_dl']] = V_an_init
        # Electrolyte electric potential = 0, so V_dl_init = V_an_init
        
    # Create dictionary to hold geometrical parameters
    geom = {}
    geom['phi_ed'] = phi_an
    geom['phi_elyte'] = 1 - phi_an
    geom['tau_ed'] = tau_an
    geom['r_p'] = r_p
    geom['d_part'] = d_part
    geom['A_surf'] = A_surf
    geom['dyInv'] = npoints/H
    geom['dr'] = geom['d_part']*0.5/nshells
    dr = geom['d_part']*0.5/nshells
    
    # Calculate the current density [A/m^2] corresponding to a C_rate of 1:
    oneC = geom['phi_ed']*anode.density_mole*H*Inputs.F/3600
    
    # Calculate the actual current density:
    # The minus sign is because we begin with the charging reaction, which 
    #   delivers negative charge to the anode:        
    i_ext = -Inputs.C_rate*oneC
    
    # Create dict to hold 'other' parameters
    params = {}
    params['npoints'] = npoints
    params['nshells'] = nshells
    params['nVars'] = nVars
    params['T'] = Inputs.T
    params['C_dl'] = C_dl
    params['X_Li_max'] = X_Li_an_max
    params['X_Li_min'] = X_Li_an_min
    params['D_Li_ed'] = D_Li_an
    params['i_ext'] = i_ext
    params['ed_flag'] = ed_flag
    
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
    for j in np.arange(0, nshells, 1):
        params['V_shell'][j] = ((j+1)**3 - (j)**3)/nshells/nshells/nshells
        
    params['sigma_eff_ed'] = sigma_an*geom['phi_ed']/geom['tau_ed']**3
    params['u_Li_elyte'] = (D_Li_elyte*geom['phi_elyte']/ct.gas_constant
          /Inputs.T/geom['tau_ed']**3)
    
    # Set up algebraic variable vector:
    algvar = np.zeros([nSV])
    for j in np.arange(0, npoints):
        offset = (j)*nVars
        # X_Li_an has a differential equation:
        for i in range(nshells):
            algvar[offset + ptr['X_ed'][i]] = 1
        
        # X_Li_elyte has a differential equation:
        algvar[offset + ptr['X_elyte']] = 1
        
        # Anode elec potential has a differential equation involving phi_an and
        #   phi_elyte:
#        algvar[offset + ptr['V_dl']] = 1
        
    # For points 1 to npoints-1, the electrolyte electric potential is governed
    #   by an algebraic equation. For the final point, the electrolyte electric
    #   potential is our reference value, and is held constant at zero. Hence
    #   it is governed by a differential equation, with dSVdt = 0
#    algvar[ptr['V_ed']] = 1
    
    # Store Cantera objects
    obj = {}
    obj['electrode'] = anode
    obj['elyte'] = elyte
    obj['surf'] = an_surf
    
    t_flag = []
    
"""========================================================================="""
"""========================================================================="""
"""========================================================================="""

class cathode():
    # Set flag so solver knows it's looking at the anode class
    ed_flag = 2
    
    # Number of nodes in the y-direction:
    npoints = 0

    # Number of "shells" in the cathode particle
    nshells = 5
    
    # Initial conditions: UPDATE TO GENERALIZE FOR MULT. SPECIES - DK 8/31/18
    LiCoO2_0 = Inputs.SOC_0; CoO2_0 = 1 - LiCoO2_0
    X_cat_init = 'LiCoO2:0.999, CoO2:0.0001'
    X_elyte_init = 'Li(e):0.98, Solvent:0.02'  #PF6-(e):7.8730103237e-2, EC(e):2.8328131770e-1'
    V_cat_init = 4
    
    ctifile = 'LiBatteryFull.cti'
    cathode_b = 'cathode'
    cathode_s = 'cathode_surf'
    elytephase = 'electrolyte'
    
    # Cutoff values for lithiation and delithiation of cathode
    X_Li_cat_max = 1 - 1e-6
    X_Li_cat_min = 1e-6
    
    # Microstructure:
    phi_cat = 0.5  # LiCoO2 volume fraction [-]
    tau_cat = 1.6  # Tortuosity - assume equal values for LiCoO2 and elyte [-]
    r_p = 5e-6     # Average pore radius [m]
    d_part = 30e-6  # Average particle diameter for LiCoO2 [m]
    H = 30e-6      # Cathode thickness [m]
    
    # Other parameters:
    C_dl = 1.5       # Double-layer capacitance [F/m^2]
    sigma_cat = 7.50    # Bulk cathode electrical conductivity [S/m]
    D_Li_elyte = 1e-10  # Bulk diffusion coefficient for Li+ in elyte [m^2/s]
    D_Li_cat = 7.5e-16  # Bulk diffusion coefficient for Li in LiCoO2 [m^2/s]
    
    nVars = nshells + 3
    
    # Pointers
    ptr = {}
    ptr['iFar'] = 0
    ptr['X_ed'] = np.arange(0, nshells)
    ptr['X_elyte'] = nshells
    ptr['V_ed'] = nshells + 1
    ptr['V_dl'] = nshells + 2
    
    # For spherical particles, the total surface area per unit volume can be
    #   calculated from the geometry. Since some particles will overlap, we take
    #   a percentage (60%, below) of this theoretical value:
    A_surf = 0.6*6*phi_cat/d_part  # Cathode/elyte interface area [m^2/m_3]

    # Create Cantera objects
    cathode = ct.Solution(ctifile, cathode_b)
    cathode.X = X_cat_init
    cathode.TP = Inputs.T, ct.one_atm
    
    X_cat_0 = cathode.X
    
    # Electrolyte
    elyte = ct.Solution(ctifile, elytephase)
    elyte.X = X_elyte_init
    elyte.TP = Inputs.T, ct.one_atm
    
    X_elyte_0 = elyte.X
    
    cat_surf = ct.Interface(ctifile, cathode_s, [cathode, elyte])
    cat_surf.TP = Inputs.T, ct.one_atm
    
    # Set up the solution vector
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])
    
    for j in range(npoints):
        offset = (j)*nVars
        SV_0[offset + ptr['X_ed']] = np.ones([nshells])*X_cat_0[1]
        SV_0[offset + ptr['X_elyte']] = X_elyte_0[0]
        SV_0[offset + ptr['V_ed']] = V_cat_init
        SV_0[offset + ptr['V_dl']] = V_cat_init
        # Electrolyte electric potential = 0, so V_dl_init = V_an_init
    
    # Create dict to hold geometrical parameters: 
    geom = {}
    geom['phi_ed'] = phi_cat
    geom['phi_elyte'] = 1 - phi_cat
    geom['tau_ed'] = tau_cat
    geom['r_p'] = r_p
    geom['d_part'] = d_part
    geom['A_surf'] = A_surf
    geom['dyInv'] = npoints/H
    geom['dr'] = geom['d_part']*0.5/nshells
    dr = geom['d_part']*0.5/nshells
        
    # Calculate the current density [A/m^2] corresponding to a C_rate of 1:
    oneC = geom['phi_ed']*cathode.density_mole*H*Inputs.F/3600
    
    # Calculate the actual current density:
    # The minus sign is because we begin with the charging reaction, which 
    #   delivers negative charge to the anode:        
    i_ext = -Inputs.C_rate*oneC
    
    # Create dict to hold 'other' parameters
    params = {}
    params['npoints'] = npoints
    params['nshells'] = nshells
    params['nVars'] = nVars
    params['T'] = Inputs.T
    params['C_dl'] = C_dl
    params['X_Li_max'] = X_Li_cat_max
    params['X_Li_min'] = X_Li_cat_min
    params['D_Li_ed'] = D_Li_cat
    params['i_ext'] = i_ext
    params['ed_flag'] = ed_flag
    
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
    for j in np.arange(0, nshells, 1):
        params['V_shell'][j] = ((j+1)**3 - (j)**3)/nshells/nshells/nshells
        
    params['sigma_eff_ed'] = sigma_cat*geom['phi_ed']/geom['tau_ed']**3
    params['u_Li_elyte'] = (D_Li_elyte*geom['phi_elyte']/ct.gas_constant
          /Inputs.T/geom['tau_ed']**3)
    
    # Set up algebraic variable vector:
    algvar = np.zeros([nSV])
    for j in range(npoints):
        offset = (j)*nVars
        # X_Li_an has a differential equation:
        for i in range(nshells):
            algvar[offset + ptr['X_ed'][i]] = 1
        
        # X_Li_elyte has a differential equation:
        algvar[offset + ptr['X_elyte']] = 1
        
        # Anode elec potential has a differential equation involving phi_an and
        #   phi_elyte:
        algvar[offset + ptr['V_dl']] = 1
        
    # For points 1 to npoints-1, the electrolyte electric potential is governed
    #   by an algebraic equation. For the final point, the electrolyte electric
    #   potential is our reference value, and is held constant at zero. Hence
    #   it is governed by a differential equation, with dSVdt = 0
#    algvar[offset + ptr['V_ed']] = 1
    
    # Store Cantera objects
    obj = {}
    obj['electrode'] = cathode
    obj['elyte'] = elyte
    obj['surf'] = cat_surf

"""========================================================================="""
"""========================================================================="""
"""========================================================================="""

class separator():
    # Number of nodes in the y-direction
    npoints = 0
    
    nVars = 2
    
    phi_elyte_init = 0
    X_Li_elyte_init = 0.98
    
    H = 100e-6  # Separator thickness [m]
    
    tau_sep = 1.6  # Tortuosity of separator
    sigma_sep = 50.0  # Bulk ionic conductivity of separator
    
    geom = {}
    geom['phi_sep'] = 0.6
    geom['phi_elyte'] = 1 - geom['phi_sep']
    geom['dyInv'] = npoints/H
    geom['tau_sep'] = tau_sep
    
    ptr = {}
    ptr['X_elyte'] = 0
    ptr['V_elyte'] = 1
    
    # Set up the solution vector
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])
    for i in np.arange(0, npoints):
        offset = i*nVars
        SV_0[offset + ptr['X_elyte']] = X_Li_elyte_init
        SV_0[offset + ptr['V_elyte']] = phi_elyte_init
    
    # Set up algebraic variable vector:
    algvar = np.zeros([nSV])
    
    for j in range(npoints):
        offset = j*nVars
        algvar[offset + ptr['X_elyte']] = 1
        algvar[offset + ptr['V_elyte']] = 0
        
    an_max = anode.npoints
    sep_max = anode.npoints + npoints
    cat_max = sep_max + cathode.npoints
        
    offset_vec = np.zeros([anode.npoints + npoints + cathode.npoints])
#    offset_vec[0:cat_max] = cathode.nVars*np.arange(0, cathode.npoints)
    offset_vec[0:an_max] = anode.nVars*np.arange(0, anode.npoints)
#    offset_vec[an_max:sep_max] = offset_vec[an_max-1]+anode.nVars+nVars*np.arange(0, npoints)
#    offset_vec[sep_max:cat_max] = offset_vec[sep_max-1]+nVars+cathode.nVars*np.arange(0, cathode.npoints)
    
    
    