# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:37:06 2018

@author: dkorff
"""

import numpy as np
import cantera as ct
from assimulo.problem import Implicit_Problem

from li_ion_battery_p2d_init import anode as an
from li_ion_battery_p2d_init import cathode as cat
from li_ion_battery_p2d_init import separator as sep
from li_ion_battery_p2d_init import Inputs
from li_ion_battery_p2d_init import anode_obj as anode
from li_ion_battery_p2d_init import anode_surf_obj as anode_s
from li_ion_battery_p2d_init import elyte_obj as elyte
from li_ion_battery_p2d_init import cathode_surf_obj as cathode_s
from li_ion_battery_p2d_init import cathode_obj as cathode
from li_ion_battery_p2d_init import conductor_obj as conductor




class Extended_Problem(Implicit_Problem):
    sw0 = True
    def Battery_Func(t, SV, SV_dot, sw):

        """================================================================="""
        """==========================INITIALIZE============================="""
        """================================================================="""

        print(t)
        nSV = len(SV)
        res = np.zeros([nSV])

        offset_vec = sep.offsets

        """        anode = an.obj['electrode']
        anode_s = an.obj['surf']
        elyte = an.obj['elyte']
        cathode = cat.obj['electrode']
        cathode_s = cat.obj['surf']"""

        nsp_an = anode.n_species; nsp_cat = cathode.n_species

        F = ct.faraday; R = ct.gas_constant; T = Inputs.T
        #sigma_eff_an = an.params['sigma_eff_ed']; dyInv = an.geom['dyInv']
        #u_Li_elyte = an.params['u_Li_elyte']; D_Li_an = an.params['D_Li_ed']
        #dr = an.dr
# %%
        """================================================================="""
        """============================ANODE================================"""
        """================================================================="""
        #  --------------------------------
        #  ANODE CURRENT COLLECTOR BOUNDARY
        #  --------------------------------

        # Looking at node 1, j=0, set THIS node conditions
        offset = an.offsets
        ptr = an.ptr
        j = 0

        N_io_m = 0
        i_io_m = 0
        i_el_m = an.i_ext

        X_an_1 = SV[offset[j] + ptr['X_ed'][-1]]
        rho_k_elyte_1 = SV[offset[j] + ptr['rho_k_elyte']]

        phi_elec_an_1 = SV[offset[j] + ptr['Phi_ed']]
        phi_elec_elyte_1 = phi_elec_an_1 - SV[offset[j] + ptr['Phi_dl']]

        anode.X = [X_an_1, 1-X_an_1]
        anode.electric_potential = phi_elec_an_1
        conductor.electric_potential = phi_elec_an_1

        #elyte.TDY = Inputs.T, np.sum(rho_k_elyte_1), rho_k_elyte_1
        elyte.Y = rho_k_elyte_1/np.sum(rho_k_elyte_1)
        elyte.electric_potential = phi_elec_elyte_1

        sdot_1 = anode_s.net_production_rates

        # Shift forward to node 2, j=1, set NEXT node conditions
        j = 1; offset = int(offset_vec[j])

        X_an_2 = SV[offset + an.ptr['X_ed'][-1]]
        rho_k_elyte_2 = SV[offset + an.ptr['rho_k_elyte']]

        phi_elec_an_2 = SV[offset + an.ptr['Phi_ed']]
        phi_elec_elyte_2 = phi_elec_an_2 - SV[offset + an.ptr['Phi_dl']]

        anode.X = [X_an_2, 1-X_an_2]
        conductor.electric_potential = phi_elec_an_2
        anode.electric_potential = phi_elec_an_2

        #elyte.TDY = Inputs.T, np.sum(rho_k_elyte_2), rho_k_elyte_2
        elyte.Y = rho_k_elyte_2/np.sum(rho_k_elyte_2)
        elyte.electric_potential = phi_elec_elyte_2

        sdot_2 = anode_s.net_production_rates

        # Shift back to node 1, j=0, set THIS node outlet conditions
        j = 0; offset = int(offset_vec[j])

        i_el_p = an.sigma_eff_ed*(phi_elec_an_1-phi_elec_an_2)*an.dyInv

        N_io_p = (-an.u_Li_elyte*elyte.density_mole*(R*T*(rho_k_elyte_2 - rho_k_elyte_1)
        + F*(phi_elec_elyte_2 - phi_elec_elyte_1))*an.dyInv)

        i_io_p = np.dot(N_io_p,Inputs.z_k_elyte)*F

        i_Far_1 = sdot_1[an.ptr['iFar']]*F*an.A_surf/an.dyInv

        X_Li = 1 - SV[offset + an.ptr['X_ed']]
        DiffFlux = np.zeros([an.nshells+1])
        DiffFlux[1:-1] = an.D_Li_ed*(X_Li[0:-1] - X_Li[1:])/an.dr
        DiffFlux[-1] = sdot_1[0]/anode.density_mole

        k_m = np.arange(0, an.nshells)/an.nshells
        k_p = np.arange(1, an.nshells+1)/an.nshells

#        print(anode_s.forward_rate_constants, phi_elec_an_1, sdot_1[an.ptr['iFar']])

        """Calculate the change in X_C6 in the particle interior.
            Note that the DiffFlux is the diffusion of lithium
            toward the particle surface, and that diffusion of Li
            into the shell decreases the amount of C6. The fluxes
            must be scaled by the shell interface surface area
            relative to the total particle surface area"""
        res[offset + an.ptr['X_ed']] = (SV_dot[offset + an.ptr['X_ed']]
        - ((DiffFlux[1:]*k_p**2 - DiffFlux[0:-1]*k_m**2)
        * an.A_surf/an.eps_ed/an.V_shell))

        """Change in electrolyte_composition"""
        res[offset + an.ptr['rho_k_elyte']] = (SV_dot[offset + an.ptr['rho_k_elyte']]
        - (((N_io_m - N_io_p)*an.dyInv + sdot_1[nsp_an]*an.A_surf)
        /elyte.density_mole/an.eps_elyte))

        """Double-layer voltage"""
        res[offset + an.ptr['Phi_dl']] = (SV_dot[offset + an.ptr['Phi_dl']]
        - (i_Far_1 + i_el_m - i_el_p)*an.dyInv/an.C_dl/an.A_surf)

        """Algebraic equation for ANODE electric potential boundary condition"""
        res[offset + an.ptr['Phi_ed']] = SV[offset + an.ptr['Phi_ed']]
#        (i_el_m - i_el_p + i_io_m - i_io_p)
#        SV_dot[offset + an.ptr['V_ed']]
# %%
        """============================ANODE================================"""
        """INTERIOR NODES"""
        for j in np.arange(2, an.npoints):

            # Save previous node outlet conditions as new inlet conditions
            N_io_m = N_io_p
            i_io_m = i_io_p
            i_el_m = i_el_p
            X_an_1 = X_an_2
            rho_k_elyte_1 = rho_k_elyte_2
            phi_elec_an_1 = phi_elec_an_2
            phi_elec_elyte_1 = phi_elec_elyte_2
            sdot_1 = sdot_2

            # Shift forward to NEXT node
            offset = int(an.offsets[j])

            X_an_2 = SV[offset + an.ptr['X_ed'][-1]]
            rho_k_elyte_2 = SV[offset + an.ptr['rho_k_elyte']]

            phi_elec_an_2 = SV[offset + an.ptr['Phi_ed']]
            phi_elec_elyte_2 = phi_elec_an_2 - SV[offset + an.ptr['Phi_dl']]

            anode.X = [X_an_2, 1-X_an_2]
            anode.electric_potential = phi_elec_an_2
            conductor.electric_potential = phi_elec_an_2

            elyte.Y = rho_k_elyte_2/np.sum(rho_k_elyte_2)
            elyte.electric_potential = phi_elec_elyte_2

            sdot_2 = anode_s.net_production_rates

            # Shift back to THIS node, set THIS node outlet conditions
            offset = int(an.offsets[j - 1])

            i_el_p = an.sigma_eff_ed*(phi_elec_an_1-phi_elec_an_2)*an.dyInv

            N_io_p = (-an.u_Li_elyte*elyte.density_mole*(R*T*(rho_k_elyte_2 - rho_k_elyte_1)
            + F*(phi_elec_elyte_2 - phi_elec_elyte_1))*an.dyInv)

            i_io_p = np.dot(N_io_p,Inputs.z_k_elyte)*F

            i_Far_1 = sdot_1[an.ptr['iFar']]*F*an.A_surf/an.dyInv

            X_Li = 1 - SV[offset + an.ptr['X_ed']]
            DiffFlux = np.zeros([an.nshells+1])
            DiffFlux[1:-1] = an.D_Li_ed*(X_Li[0:-1] - X_Li[1:])/an.dr
            DiffFlux[-1] = sdot_1[0]/anode.density_mole

            """Calculate the change in X_C6 in the particle interior."""
            res[offset + an.ptr['X_ed']] = (SV_dot[offset + an.ptr['X_ed']]
            - ((DiffFlux[1:]*k_p**2 - DiffFlux[0:-1]*k_m**2)
            * an.A_surf/an.eps_ed/an.V_shell))

            """Change in electrolyte_composition"""
            res[offset + an.ptr['rho_k_elyte']] = (SV_dot[offset + an.ptr['rho_k_elyte']]
            - (((N_io_m - N_io_p)*an.dyInv + sdot_1[nsp_an]*an.A_surf)
            /elyte.density_mole/an.eps_elyte))

            """Double-layer voltage"""
            res[offset + an.ptr['Phi_dl']] = (SV_dot[offset + an.ptr['Phi_dl']]
            - (i_Far_1 + i_el_m - i_el_p)*an.dyInv/an.C_dl/an.A_surf)

            """Algebraic equation for ANODE electric potential boundary condition"""
            res[offset + an.ptr['Phi_ed']] = (i_el_m - i_el_p + i_io_m - i_io_p)
# %%
        """============================ANODE================================"""
        """Separator boundary"""
        # Save previous node outlet conditions as new inlet conditions
        N_io_m = N_io_p
        i_io_m = i_io_p
        i_el_m = i_el_p
        X_an_1 = X_an_2
        rho_k_elyte_1 = rho_k_elyte_2
        phi_elec_an_1 = phi_elec_an_2
        phi_elec_elyte_1 = phi_elec_elyte_2
        sdot_1 = sdot_2

        # Shift forward to NEXT node (first separator node)
#        j = an.npoints; offset = int(offset_vec[j])
#
#        X_elyte_2 = SV[offset + sep.ptr['X_elyte']]
#
#        phi_elec_elyte_2 = SV[offset + sep.ptr['V_elyte']]

        # Shift back to THIS node, set THIS node outlet conditions
        i_el_p = 0

#        N_io_p = (-u_Li_elyte*elyte.density_mole*(R*T*(X_elyte_2 - X_elyte_1)
#        + F*(phi_elec_elyte_2 - phi_elec_elyte_1))*dyInv)
#
#        i_io_p = N_io_p*F

        # Set j to final ANODE node
        j = an.npoints-1; offset = int(an.offsets[j])

        i_Far_1 = sdot_1[an.ptr['iFar']]*F*an.A_surf/an.dyInv

        i_io_p = an.i_ext
        #THIS IS TEMPORARY, NON-GENERALIZED CODE:
        N_io_p = np.zeros_like(N_io_p)
        N_io_p[2] = i_io_p/F

        X_Li = 1 - SV[offset + an.ptr['X_ed']]
        DiffFlux = np.zeros([an.nshells+1])
        DiffFlux[1:-1] = an.D_Li_ed*(X_Li[0:-1] - X_Li[1:])/an.dr
        DiffFlux[-1] = sdot_1[0]/anode.density_mole

        """Calculate the change in X_C6 in the particle interior."""
        res[offset + an.ptr['X_ed']] = (SV_dot[offset + an.ptr['X_ed']]
        - ((DiffFlux[1:]*k_p**2 - DiffFlux[0:-1]*k_m**2)
        * an.A_surf/an.eps_ed/an.V_shell))

        """Change in electrolyte_composition"""
        res[offset + an.ptr['rho_k_elyte']] = (SV_dot[offset + an.ptr['rho_k_elyte']]
        - (((N_io_m - N_io_p)*an.dyInv + sdot_1[nsp_an]*an.A_surf)
        /elyte.density_mole/an.eps_elyte))

        """Double-layer voltage"""
        res[offset + an.ptr['Phi_dl']] = (SV_dot[offset + an.ptr['Phi_dl']]
        - (i_Far_1 + i_el_m - i_el_p)*an.dyInv/an.C_dl/an.A_surf)

        """Algebraic equation for ANODE electric potential boundary condition"""
        res[offset + an.ptr['Phi_ed']] = (i_el_m - i_el_p + i_io_m - i_io_p)
# %%
        """================================================================="""
        """==========================SEPARATOR=============================="""
        """================================================================="""
#        for j in np.arange(an.npoints+1, sep.sep_max):
#
#            X_elyte_1 = X_elyte_2
#            phi_elec_elyte_1 = phi_elec_elyte_2
#            N_io_m = N_io_p
#            i_io_m = i_io_p

            # Shift forward to NEXT node
#            offset = int(offset_vec[j])
#
#            X_elyte_2 = SV[offset + sep.ptr['X_elyte']]
#            phi_elec_elyte_2 = SV[offset + sep.ptr['V_elyte']]

            # Step back to THIS node to calculate outlet flux
#            offset = int(offset_vec[j-1])

#            N_io_p = (-u_Li_elyte*elyte.density_mole*(R*T*(X_elyte_2 - X_elyte_1)
#            + F*(phi_elec_elyte_2 - phi_elec_elyte_1))*sep.geom['dyInv'])
#
#            i_io_p = N_io_p*F

#            i_io_p = an.params['i_ext']
#            N_io_p = i_io_p/F
#
#            """Change in electrolyte_composition"""
#            res[offset + sep.ptr['X_elyte']] = (SV_dot[offset + sep.ptr['X_elyte']]
#            - (((N_io_m - N_io_p)*dyInv)/elyte.density_mole/sep.geom['phi_elyte']))
#
#            """Charge neutrality enforced"""
#            res[offset + sep.ptr['V_elyte']] = (i_io_m - i_io_p)
# %%
        # Looking at LAST node in separator
#        X_elyte_1 = X_elyte_2
#        phi_elec_elyte_1 = phi_elec_elyte_2
#        N_io_m = N_io_p
#        i_io_m = i_io_p

        # Shift forward to NEXT node, first cathode node
#        j = sep.sep_max; offset = int(offset_vec[j])
#
#        X_cat_2 = SV[offset + cat.ptr['X_ed'][-1]]
#        X_elyte_2 = SV[offset + cat.ptr['X_elyte']]
#
#        phi_elec_cat_2 = SV[offset + cat.ptr['V_ed']]
#        phi_elec_elyte_2 = phi_elec_cat_2 - SV[offset + cat.ptr['V_dl']]
#
#        cathode.X = [1-X_cat_2, X_cat_2]
#        cathode.electric_potential = phi_elec_cat_2
#
#        elyte.X = [X_elyte_2, 7.8730103237e-2, 2.8328131770e-1]
#        elyte.electric_potential = phi_elec_elyte_2
#
#        sdot_2 = cathode_s.net_production_rates

        # Shift back to THIS node (last separator node)
#        j = sep.sep_max-1; offset = int(offset_vec[j])
#
#        i_el_p = 0

#        N_io_p = (-u_Li_elyte*elyte.density_mole*(R*T*(X_elyte_2 - X_elyte_1)
#            + F*(phi_elec_elyte_2 - phi_elec_elyte_1))*sep.geom['dyInv'])
#
#        i_io_p = N_io_p*F

#        i_io_p = an.params['i_ext']
#        N_io_p = i_io_p/F
#
#        """Change in electrolyte_composition"""
#        res[offset + sep.ptr['X_elyte']] = (SV_dot[offset + sep.ptr['X_elyte']]
#        - (((N_io_m - N_io_p)*dyInv)/elyte.density_mole/sep.geom['phi_elyte']))
#
#        """Charge neutrality enforced"""
#        res[offset + sep.ptr['V_elyte']] = (i_io_m - i_io_p)
#        print(SV, res)
#        SV[offset + sep.ptr['V_elyte']]
#       (i_io_m - i_io_p)

# %%
        """================================================================="""
        """===========================CATHODE==============================="""
        """================================================================="""

                # Alrighty, CATHODE time

#        sigma_eff_cat = cat.params['sigma_eff_ed']; dyInv = cat.geom['dyInv']
#        D_Li_cat = cat.params['D_Li_ed']
#
#        i_io_m = i_io_p
#        N_io_m = N_io_p
#        i_el_m = i_el_p
#        X_cat_1 = X_cat_2
#        X_elyte_1 = X_elyte_2
#        phi_elec_cat_1 = phi_elec_cat_2
#        phi_elec_elyte_1 = phi_elec_elyte_2
#        sdot_1 = sdot_2
#        j = sep.cat_max-1; offset = int(offset_vec[j])
#        i_el_p = -an.params['i_ext']
#        i_io_p = 0
#        N_io_p = i_io_p/F
#
#        i_Far_1 = sdot_1[cat.ptr['iFar']]*F*cat.geom['A_surf']/dyInv
#        print(cathode_s.forward_rate_constants, phi_elec_cat_1, sdot_1[cat.ptr['iFar']])
#        X_Li = 1 - SV[offset + cat.ptr['X_ed']]
#        DiffFlux = np.zeros([cat.nshells+1])
#        DiffFlux[1:-1] = D_Li_cat*(X_Li[0:-1] - X_Li[1:])/dr
#        DiffFlux[-1] = sdot_1[0]/cathode.density_mole
#
#        """Calculate the change in CoO2 in the particle interior."""
#        res[offset + cat.ptr['X_ed']] = (SV_dot[offset + cat.ptr['X_ed']]
#        - ((DiffFlux[1:]*k_p**2 - DiffFlux[0:-1]*k_m**2)
#        * cat.geom['A_surf']/cat.geom['phi_ed']/cat.params['V_shell']))
#
#        """Change in electrolyte_composition"""
#        res[offset + cat.ptr['X_elyte']] = (SV_dot[offset + cat.ptr['X_elyte']]
#        - (((N_io_m - N_io_p)*dyInv + sdot_1[nsp_cat]*cat.geom['A_surf'])
#        /elyte.density_mole/cat.geom['phi_elyte']))
#
#        """Double-layer voltage"""
#        res[offset + cat.ptr['V_dl']] = (SV_dot[offset + cat.ptr['V_dl']]
#        - (i_Far_1 + i_el_m - i_el_p)*dyInv/cat.params['C_dl']/cat.geom['A_surf'])
#
#        """Algebraic equation for CATHODE electric potential boundary condition"""
#        res[offset + cat.ptr['V_ed']] = (i_el_m - i_el_p + i_io_m - i_io_p)
#        print(SV, res)


#        for j in np.arange(an.npoints + sep.npoints, sep.cat_max-1):
#            N_io_m = N_io_p
#            i_io_m = i_io_p
#            i_el_m = i_el_p
#            X_cat_1 = X_cat_2
#            X_elyte_1 = X_elyte_2
#            phi_elec_cat_1 = phi_elec_cat_2
#            phi_elec_elyte_1 = phi_elec_elyte_2
#            sdot_1 = sdot_2
#
#            # Look at NEXT node
#            offset = int(offset_vec[j])
#
#            X_cat_2 = SV[offset + cat.ptr['X_ed'][-1]]
#            X_elyte_2 = SV[offset + cat.ptr['X_elyte']]
#
#            phi_elec_cat_2 = SV[offset + cat.ptr['V_ed']]
#            phi_elec_elyte_2 = phi_elec_cat_2 - SV[offset + cat.ptr['V_dl']]
#
#            cathode.X = [1-X_cat_2, X_cat_2]
#            cathode.electric_potential = phi_elec_cat_2
#
#            elyte.X = [X_elyte_2, 1-X_elyte_2]
#            elyte.electric_potential = phi_elec_elyte_2
#
#            sdot_2 = cathode_s.net_production_rates
#
#            # Shift back to THIS node, set THIS node outlet conditions
#            offset = int(offset_vec[j-1])
#
#            i_el_p = sigma_eff_cat*(phi_elec_cat_1 - phi_elec_cat_2)*dyInv
#
#            N_io_p = (-u_Li_elyte*elyte.density_mole*(R*T*(X_elyte_2 - X_elyte_1)
#            + F*(phi_elec_elyte_2 - phi_elec_elyte_1))*dyInv)
#
#            i_io_p = N_io_p*F
#
#            i_Far_1 = sdot_1[cat.ptr['iFar']]*F*cat.geom['A_surf']/dyInv
#
#            X_Li = 1 - SV[offset + cat.ptr['X_ed']]
#            DiffFlux = np.zeros([cat.nshells+1])
#            DiffFlux[1:-1] = D_Li_cat*(X_Li[0:-1] - X_Li[1:])/dr
#            DiffFlux[-1] = sdot_1[1]/cathode.density_mole
#
#            """Calculate the change in CoO2 in the particle interior."""
#            res[offset + cat.ptr['X_ed']] = (SV_dot[offset + cat.ptr['X_ed']])
#            """- ((DiffFlux[1:]*k_p**2 - DiffFlux[0:-1]*k_m**2)
#            * cat.geom['A_surf']/cat.geom['phi_ed']/cat.params['V_shell']))"""
#
#            """Change in electrolyte_composition"""
#            res[offset + cat.ptr['X_elyte']] = (SV_dot[offset + cat.ptr['X_elyte']])
#            """- (((N_io_m - N_io_p)*dyInv + sdot_1[nsp_cat]*cat.geom['A_surf'])
#            /elyte.density_mole/cat.geom['phi_elyte']))"""
#
#            """Double-layer voltage"""
#            res[offset + cat.ptr['V_dl']] = (SV_dot[offset + cat.ptr['V_dl']]
#            - (i_Far_1 + i_el_m - i_el_p)*dyInv/cat.params['C_dl']/cat.geom['A_surf'])
#
#            """Algebraic equation for CATHODE electric potential boundary condition"""
#            res[offset + cat.ptr['V_ed']] = (i_el_m - i_el_p + i_io_m - i_io_p)

# %%
        """=========================CATHODE============================="""
        """current collector boundary"""

#        N_io_m = N_io_p
#        i_io_m = i_io_p
#        i_el_m = i_el_p
#        X_cat_1 = X_cat_2
#        X_elyte_1 = X_elyte_2
#        phi_elec_cat_1 = phi_elec_cat_2
#        phi_elec_elyte_1 = phi_elec_elyte_2
#        sdot_1 = sdot_2
#
#        # Set THIS node outlet conditions (last node BCs)
#        j = sep.cat_max - 1; offset = int(offset_vec[j])
#        i_io_p = 0
#        N_io_p = 0
#        i_el_p = cat.params['i_ext']
#
#        i_Far_1 = sdot_1[cat.ptr['iFar']]*F*cat.geom['A_surf']/dyInv
#
#        X_Li = 1 - SV[offset + cat.ptr['X_ed']]
#        DiffFlux = np.zeros([cat.nshells+1])
#        DiffFlux[1:-1] = D_Li_cat*(X_Li[0:-1] - X_Li[1:])/dr
#        DiffFlux[-1] = sdot_1[1]/cathode.density_mole
#
#        """Calculate the change in CoO2 in the particle interior."""
#        res[offset + cat.ptr['X_ed']] = (SV_dot[offset + cat.ptr['X_ed']])
#        """- ((DiffFlux[1:]*k_p**2 - DiffFlux[0:-1]*k_m**2)
#        * cat.geom['A_surf']/cat.geom['phi_ed']/cat.params['V_shell']))"""
#
#        """Change in electrolyte_composition"""
#        res[offset + cat.ptr['X_elyte']] = (SV_dot[offset + cat.ptr['X_elyte']])
#        """- (((N_io_m - N_io_p)*dyInv + sdot_1[nsp_cat]*cat.geom['A_surf'])
#        /elyte.density_mole/cat.geom['phi_elyte']))"""
#
#        """Double-layer voltage"""
#        res[offset + cat.ptr['V_dl']] = (SV_dot[offset + cat.ptr['V_dl']]
#        - (i_Far_1 + i_el_m - i_el_p)*dyInv/cat.params['C_dl']/cat.geom['A_surf'])
#
#        """Algebraic equation for CATHODE electric potential boundary condition"""
#        res[offset + cat.ptr['V_ed']] = (i_el_m - i_el_p + i_io_m - i_io_p)

        return res

    """====================================================================="""
    """====================================================================="""
    """====================================================================="""
# %%
    def state_events(self, t, y, yd, sw):
        event1 = np.zeros([an.params['npoints']])
        event2 = np.zeros([an.params['npoints']])
        event3 = np.zeros([an.params['nshells']])
        event4 = np.zeros([an.params['nshells']])

        for j in np.arange(0, an.params['npoints']):
            offset = j*an.params['nVars']

            event1[j] = (y[offset + an.ptr['V_dl']])
            event2[j] = (1 - y[offset + an.ptr['V_dl']])

            for i in np.arange(0, an.params['nshells']):
                event3[i] = y[offset + an.ptr['X_ed'][i]] - (1 - an.params['X_Li_max'])
                event4[i] = (((1 - an.params['X_Li_min']) - y[offset + an.ptr['X_ed'][i]]))

        event5 = np.zeros([cat.params['npoints']])
        event6 = np.zeros([cat.params['npoints']])
        event7 = np.zeros([cat.params['nshells']])
        event8 = np.zeros([cat.params['nshells']])

        for j in np.arange(0, cat.params['npoints']):
            offset = j*cat.params['nVars'] + an.npoints*an.nVars + sep.npoints*sep.nVars

            event5[j] = (y[offset + cat.ptr['V_dl']])
            event6[j] = (y[offset + cat.ptr['V_dl']] - 5)

            for i in np.arange(0, cat.params['nshells']):
                event7[i] = y[offset + cat.ptr['X_ed'][i]] - (1 - cat.params['X_Li_max'])
                event8[i] = (1 - cat.params['X_Li_min']) - y[offset + cat.ptr['X_ed'][i]]

        event9 = np.zeros([sep.npoints])
        event10 = np.zeros([sep.npoints])
        for j in np.arange(0, sep.npoints):
            offset = an.npoints*an.nVars
            event9[j] = 1 - y[offset + sep.ptr['X_elyte']]
            event10[j] = y[offset + sep.ptr['X_elyte']]

        events = np.concatenate((event1, event2, event3, event4, event5, event6,
                                 event7, event8, event9, event10))

        return events

    """====================================================================="""

    def handle_event(self, solver, event_info):
        while True:
            self.event_switch(solver, event_info)
            self.init_mode(solver)

            if not True in event_info:
                break

    def event_switch(self, solver, event_info):
        if not all(event_info):
            solver.sw = [not solver.sw]

    def init_mode(self, solver):
        an.t_flag = solver.t
        if an.params['i_ext'] != 0:
            an.params['i_ext'] = 0
            cat.params['i_ext'] = 0
