# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 12:47:07 2018

@author: Owner

"""

###############################################################################
#
#   1-D Battery model
#
#       -S. DeCaluwe, Colorado School of Mines
#       -November, 2014
#
#       -D. Korff, Colorado School of Mines
#       -2018, Translated to Python and updated
#
###############################################################################

import numpy as np
import time
from matplotlib import pyplot as plt

from Battery_Func_module_full import Equilibrate
from Battery_Func_module_full import Charge
from Battery_Func_module_full import Re_equilibrate
from Battery_Func_module_full import Discharge
from Battery_Func_module_full import Extended_Problem

from Post_process_full import Label_Columns
from Post_process_full import tag_strings
from Post_process_full import plot_sims

from assimulo.solvers import IDA

from User_Inputs_Full import Inputs as inp
from User_Inputs_Full import anode
from User_Inputs_Full import cathode
from User_Inputs_Full import separator

def main():
    plt.close()
    t_count = time.time()
    
    atol1 = 1e-9
    rtol1 = 1e-5
    
    atol2 = 1e-3
    rtol2 = 1e-9
    
    atol3 = 1e-9
    rtol3 = 1e-5
    
    atol4 = atol2
    rtol4 = rtol2
    
#    inp.C_rate = C_rate
    
    # Calculate the time span, which should be enough to charge or discharge fully
    #   (i.e. 3600 seconds divided by the C-rate):
    t_f = 3600/inp.C_rate
    t_0 = 0
#    SV_0 = anode.SV_0
#    SV_0 = np.concatenate((anode.SV_0, separator.SV_0, cathode.SV_0))
    SV_0 = np.concatenate((anode.SV_0, separator.SV_0)) 
    SV_dot_0 = SV_0*0
            
    """----------Equilibration----------"""
    
    # Equilibrate by integrating at zero current:
    print('\nEquilibrating...')

    # Create problem instance
    Battery_equil = Extended_Problem(Equilibrate, SV_0, SV_dot_0, t_0)
    Battery_equil.external_event_detection = True
#    algvar = anode.algvar
#    algvar = np.concatenate((anode.algvar, separator.algvar, cathode.algvar))
    algvar = np.concatenate((anode.algvar, separator.algvar))
#    Battery_equil.algvar = algvar
    
    # Simulation parameters
    equil_sim = IDA(Battery_equil)           # Create simulation instance
    equil_sim.atol = atol1                 # Solver absolute tolerance
    equil_sim.rtol = rtol1                  # Solver relative tolerance
    equil_sim.verbosity = 50
    equil_sim.make_consistent('IDA_YA_YDP_INIT')
        
    t_eq, SV_eq, SV_dot_eq = equil_sim.simulate(t_f)
    
    # Put solution into pandas dataframe with labeled columns
    SV_eq_df = Label_Columns(t_eq, SV_eq)
    
    # Obtain tag strings for dataframe columns
    tags = tag_strings(SV_eq_df)

    print('Done equilibrating\n')
    """---------------------------------"""
    
    """------------Charging-------------"""
    print('\nCharging...')

    # New initial conditions are the final equilibrium conditions
    t_0 = t_eq[-1]
    t_f1 = t_0 + t_f
    SV_0 = SV_eq[-1, :]
    SV_dot_0 = SV_dot_eq[-1 :]
    
    # Charge the battery
    anode.params['i_ext'] = anode.i_ext
    cathode.params['i_ext'] = -cathode.i_ext
    
    # Create problem instance 
    Battery_charge = Extended_Problem(Charge, SV_0, SV_dot_0, t_0)
    Battery_charge.external_event_detection = True
    Battery_charge.algvar = algvar    
    
    # Simulation parameters
    charge_sim = IDA(Battery_charge)
    charge_sim.atol = atol2
    charge_sim.rtol = rtol2
    charge_sim.verbosity = 50
    charge_sim.make_consistent('IDA_YA_YDP_INIT')
#    charge_sim.maxh = 0.5
        
    t_charge, SV_charge, SV_dot_charge = charge_sim.simulate(t_f1)
    
    t_flag1 = anode.t_flag
    
    SV_charge_df = Label_Columns(t_charge, SV_charge) 
    
    print('Done charging\n')
    """---------------------------------"""
    
    """------------Re_equilibrating-------------"""
    
    # New initial conditions are the final charge conditions
    t_0 = t_charge[-1]  #anode.t_flag
    t_f2 = t_0 + t_f
    SV_0 = SV_charge[-1, :]
    SV_dot_0 = SV_dot_charge[-1, :]
    
    # Equilibrate again. Note - this is a specific choice to reflect 
    #   equilibration after the charging steps. We may want, at times, to
    #   simulate a situation where the battery is not equilibrated between
    #   charge and discharge, or is equilibrated for a shorter amount of time.
    
    print('\nRe-equilibrating...')
    
    Battery_re_equil = Extended_Problem(Re_equilibrate, SV_0, SV_dot_0, t_0)
    Battery_re_equil.external_event_detection = True
    Battery_re_equil.algvar = algvar
    
    # Simulation parameters
    re_equil_sim = IDA(Battery_re_equil)
    re_equil_sim.atol = atol3
    re_equil_sim.rtol = rtol3
    re_equil_sim.verbosity = 50
    re_equil_sim.make_consistent('IDA_YA_YDP_INIT')
        
    t_req, SV_req, SV_dot_req = re_equil_sim.simulate(t_f2)
    
    SV_req_df = Label_Columns(t_req, SV_req) 

    print('Done re-equilibrating\n')
    
    """---------------------------------"""
    
    """------------Discharging-------------"""
    
    print('\nDischarging...')
    
    t_0 = t_req[-1]
    t_f3 = t_f2 + t_f
    SV_0 = SV_req[-1, :]
    SV_dot_0 = SV_dot_req[-1, :]
    
    anode.params['i_ext'] = -anode.i_ext
    cathode.params['i_ext'] = cathode.i_ext
    
    Battery_discharge = Extended_Problem(Discharge, SV_0, SV_dot_0, t_0)
    Battery_discharge.external_event_detection = True
    Battery_discharge.algvar = algvar
    
    # Simulation parameters
    Battery_discharge = IDA(Battery_discharge)
    Battery_discharge.atol = atol4
    Battery_discharge.rtol = rtol4
    Battery_discharge.verbosity = 50
    Battery_discharge.make_consistent('IDA_YA_YDP_INIT')

    t_discharge, SV_discharge, SV_dot_discharge = Battery_discharge.simulate(t_f3)
    
    t_flag2 = anode.t_flag
    
    SV_discharge_df = Label_Columns(t_discharge, SV_discharge)
    
    print('Done discharging\n')
    
    """---------------------------------"""
    
    SV_dict = {}
    SV_dict['SV_eq'] = SV_eq_df
    SV_dict['SV_charge'] = SV_charge_df
    SV_dict['SV_req'] = SV_req_df
    SV_dict['SV_discharge'] = SV_discharge_df
    
#    SV_charge_df = SV_charge_df.loc[SV_charge_df['Time'] <= t_flag1]
#    SV_discharge_df = SV_discharge_df.loc[SV_discharge_df['Time'] <= t_flag2]
    
#    dt_cap = t_flag1 - t_eq[-1]
#    t_80 = 0.8*dt_cap + t_eq[-1]
    
#    elyte80_charge = SV_charge_df.loc[SV_charge_df['Time'] <= t_80]
#    elyte75_charge = SV_charge_df.loc[SV_charge_df['X_an15'] <= 0.20]
#    el80 = elyte80_charge[X_elyte].iloc[0]
#    el80 = list(el80)
#    el80.append(inp.C_rate)
    
#    import openpyxl
#    book = openpyxl.load_workbook('Elyte_depth_profiles.xlsx')
#    sheet = book.active
#    sheet.append(el80)
#    book.save('Elyte_depth_profiles.xlsx')
       
#    print('Elyte Li+ at 80% SOC = ', elyte75_charge, '\n')
    
    plot_sims(tags['V_an'], tags['X_an'], tags['X_el_an'], SV_dict, t_f1, t_f3, t_flag1, t_flag2)
#    plot_sims(tags['V_cat'], tags['X_cat'], tags['X_el_cat'], SV_dict, t_f1, t_f3, t_flag1, t_flag2)
    
#    fig, axes = plt.subplots(sharey = "row", figsize = (18, 8), nrows=1, ncols=4)
#    plt.subplots_adjust(wspace=0, hspace=0.5)
#    
#    V_cell_plot = V_cell_eq.plot(ax = axes[0])
#    V_cell_plot.set_title('Equilibration') 
#    V_cell_plot.set_ylabel('Cell Voltage')
#    V_cell_plot.legend().set_visible(False)
#    
#    V_cell_plot = V_cell_charge.plot(ax = axes[1])
#    V_cell_plot.set_title('Charging') 
#    V_cell_plot.set_ylabel('Cell Voltage')
#    V_cell_plot.legend().set_visible(False)
#    
#    V_cell_plot = V_cell_req.plot(ax = axes[2])
#    V_cell_plot.set_title('Re-Equilibration') 
#    V_cell_plot.set_ylabel('Cell Voltage')
#    V_cell_plot.legend().set_visible(False)
#    
#    V_cell_plot = V_cell_discharge.plot(ax = axes[3])
#    V_cell_plot.set_title('Discharge') 
#    V_cell_plot.set_ylabel('Cell Voltage')
#    V_cell_plot.legend().set_visible(False)
    
    
    """---------------------------------"""    
    # Post processing
    V_charge = np.array(SV_charge_df['V_dl1'])
    V_discharge = np.array(SV_discharge_df['V_dl1'])
    t_charge = np.array(SV_charge_df['Time'])
    t_discharge = np.array(SV_discharge_df['Time'])
    dt_charge = t_charge - t_charge[0]
    dt_discharge = t_discharge - t_discharge[0]
    
    # Plot charge-discharge curve
    Capacity_charge = -dt_charge*anode.i_ext/3600  # A-h/m^2
    Capacity_discharge = -dt_discharge*anode.i_ext/3600  # A-h/m^2
    
    fig1, ax1 = plt.subplots()
    ax1.plot(Capacity_charge, V_charge)
    ax1.plot(Capacity_discharge, V_discharge)
#    ax1.set_ylim((0, 1.2))
#    ax1.set_xlim((-0.1, 18.5))
    ax1.set_xlabel("Capacity [Ah/m^2]")
    ax1.set_ylabel("Voltage [V]")
    ax1.legend(("Charge capacity", "Discharge capacity"))
    
    # Calculate battery energy storage/recovery and calculate round-trip
    #   efficiency. Anode voltage is referenced to its initial equilibrium
    #   value (i.e. in the discharged state).
    
    # NOTE: This is in W-h/m^2, per unit area of battery. For the specific
    #   capacity, you want W-h/g of anode material.
#    E_stored = 0
#    E_recovered = 0
    
#    for k in np.arange(1, len(t_charge)):
#        E_stored = (E_stored - (anode.V_cathode - 0.5*(V_charge[k] + V_charge[k-1]))
#                    *ep.i_ext*(dt_charge[k] - dt_charge[k-1]))
#        
#    for k in np.arange(1, len(t_discharge)):
#        E_recovered = (E_recovered - (ep.V_cathode - 
#                        0.5*(V_discharge[k] + V_discharge[k-1]))
#                        *ep.i_ext*(dt_discharge[k] - dt_discharge[k-1]))
#    
#    Cap_recovered = Capacity_discharge[-1]
#    Eta_RT = E_recovered/E_stored
    
#    print('Cap_recovered = ', Cap_recovered, '\n')
#    print(E_stored, '\n')
#    print(E_recovered, '\n')
#    print('Eta_RT = ', Eta_RT, '\n')
    
    elapsed = time.time() - t_count
    print('t_cpu=', elapsed, '\n')
    
    return t_req, SV_req, SV_dot_req

if __name__ == "__main__":
    t, SV, SV_dot = main()
