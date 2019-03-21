# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:59:27 2018

@author: dkorff
"""

from li_ion_battery_p2d_init import anode
from li_ion_battery_p2d_init import cathode
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from li_ion_battery_p2d_init import separator as sep

def plot_sims(V, X, X_el, SV_dict, t_f1, t_f3, t_flag1, t_flag2):

    SV_eq = SV_dict['SV_eq']; t_eq = SV_eq['Time']
    SV_charge = SV_dict['SV_charge']; t_charge = SV_charge['Time']
    SV_req = SV_dict['SV_req']; t_req = SV_req['Time']
    SV_discharge = SV_dict['SV_discharge']; t_discharge = SV_discharge['Time']

    """====================================================================="""
    """Plot equilibration"""
    fig, axes = plt.subplots(sharey = "row", figsize = (18, 8), nrows=3, ncols=4)
    plt.subplots_adjust(wspace=0, hspace=0.5)

    # Plot anode and double-layer voltage
    SV_plot = SV_eq.plot(x = 'Time', y = V[1],
                            ax = axes[0,0], xlim = [0, t_eq.iloc[-1]])
    SV_plot.set_title('Equilibration')
    SV_plot.set_ylabel('Voltages')
    SV_plot.legend().set_visible(False)
    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    #Plot anode composition
    SV_plot = SV_eq.plot(x = 'Time', y = X,
                            ax = axes[1,0], xlim = [0, t_eq.iloc[-1]], ylim = [0, 1.2])
    SV_plot.set_ylabel('Anode composition')
    SV_plot.legend().set_visible(False)
    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Plot elyte composition
    SV_plot = SV_eq.plot(x = 'Time', y = X_el,
                            ax = axes[2, 0], xlim = [0, t_eq.iloc[-1]], ylim = [0, 1.2])
    SV_plot.set_ylabel('Electrolyte composition')
    SV_plot.legend().set_visible(False)
    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    """====================================================================="""
    """Plot charging"""
    # Plot anode and double-layer voltage
    SV_charge_plot = SV_charge.plot(x = 'Time', y = V,
                                        ax = axes[0,1],
                                       xlim = [t_eq.iloc[-1], t_charge.iloc[-1]])
    SV_charge_plot.legend().set_visible(False)
    SV_charge_plot.set_ylabel('Potential [V]')
    SV_charge_plot.set_title('Charging')
    SV_charge_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


    # Plot anode composition
    SV_charge_plot = SV_charge.plot(x = 'Time', y = X,
                            ax = axes[1, 1], xlim = [t_eq.iloc[-1], t_charge.iloc[-1]], ylim = [0, 1.2])
    SV_charge_plot.legend().set_visible(False)
    SV_charge_plot.set_ylabel('[X_C6 in anode]')
    SV_charge_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Plot elyte composition
    SV_charge_plot = SV_charge.plot(x = 'Time', y = X_el,
                            ax = axes[2, 1], xlim = [t_eq.iloc[-1], t_charge.iloc[-1]], ylim = [0, 1.2])
    SV_charge_plot.legend().set_visible(False)
    SV_charge_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    """====================================================================="""
    """Plot re-equilibration"""
    # Plot anode and double-layer voltage
    SV_plot = SV_req.plot(x = 'Time', y = V,
                             ax = axes[0, 2],
                             xlim = [t_req.iloc[0], t_req.iloc[-1]])
    SV_plot.legend().set_visible(False)
    SV_plot.set_title('Re-equilibration')
    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Plot anode composition
    SV_plot = SV_req.plot(x = 'Time', y = X,
                            ax = axes[1,2], xlim = [t_req.iloc[0], t_req.iloc[-1]], ylim = [0, 1.2])
    SV_plot.legend().set_visible(False)
    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Plot elyte composition
    SV_plot = SV_req.plot(x = 'Time', y = X_el,
                            ax = axes[2, 2], xlim = [t_req.iloc[0], t_req.iloc[-1]], ylim = [0, 1.2])
    SV_plot.legend().set_visible(False)
    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    """====================================================================="""
    """Plot discharging"""
    # Plot anode and double-layer voltage
    SV_discharge_plot = SV_discharge.plot(x = 'Time',
                                             y = V,
                                             ax = axes[0,3],
                                             xlim = [t_req.iloc[-1], t_discharge.iloc[-1]])
    SV_discharge_plot.legend(loc = 2, bbox_to_anchor = (1.05, 1),
                             ncol = 1, borderaxespad = 0, frameon = False)
    SV_discharge_plot.set_title('Discharging')
    SV_discharge_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Plot anode composition
    SV_discharge_plot = SV_discharge.plot(x = 'Time', y = X,
                                          ax = axes[1,3],
                                          xlim = [t_req.iloc[-1], t_discharge.iloc[-1]], ylim = [0, 1.2])
    SV_discharge_plot.legend(loc = 'best', bbox_to_anchor = (0.99, 1.05),
                             ncol = 1, borderaxespad = 1.0, frameon = False)
    SV_discharge_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    # Plot elyte composition
    SV_discharge_plot = SV_discharge.plot(x = 'Time', y = X_el,
                                          ax = axes[2, 3],
                                          xlim = [t_req.iloc[-1], t_discharge.iloc[-1]], ylim = [0, 1.2])
    SV_discharge_plot.legend(loc = 2, bbox_to_anchor = (1.05, 1),
                             ncol = 1, borderaxespad = 0, frameon = False)
    SV_discharge_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


"""========================================================================="""

def tag_strings(SV):

    SV_eq_labels = SV.columns.values.tolist()

    X_an = []
    X_cat = []
    X_el_an = []
    X_el_cat = []
    V_an = []
    V_cat = []

    for j in np.arange(0, anode.npoints):
        offset = int(sep.offset_vec[j])
        V_an[0+offset:1+offset] = SV_eq_labels[anode.ptr['V_ed']+offset:anode.ptr['V_dl']+offset+1]
        X_an[0+offset:anode.params['nshells']+offset] = SV_eq_labels[0+offset:anode.params['nshells']+offset]
        X_el_an.append(SV_eq_labels[offset + anode.params['nshells']])

    for j in np.arange(sep.sep_max, sep.cat_max):
        offset = int(sep.offset_vec[j])
        V_cat[0+offset:1+offset] = SV_eq_labels[cathode.ptr['V_ed']+offset:cathode.ptr['V_dl']+offset+1]
        X_cat[0+offset:cathode.params['nshells']+offset] = SV_eq_labels[0+offset:cathode.params['nshells']+offset]
        X_el_cat.append(SV_eq_labels[offset + cathode.params['nshells']])

    tags = {}
    tags['V_an'] = V_an; tags['V_cat'] = V_cat; tags['X_an'] = X_an; tags['X_cat'] = X_cat
    tags['X_el_an'] = X_el_an; tags['X_el_cat'] = X_el_cat

    return tags

"""========================================================================="""

def Label_Columns(t, SV):
    # For readability, store local copies of some variables:
    params = anode.params

    t_df = pd.DataFrame(t)
    newcols = {0: anode.nVars*anode.npoints + sep.nVars*sep.npoints + cathode.nVars*cathode.npoints}
    t_df.rename(columns = newcols, inplace = True)

    SV_df = pd.DataFrame(SV)
    SV_df = pd.concat((SV_df, t_df), axis = 1)

    for j in np.arange(0, anode.npoints):
        offset = sep.offset_vec[j]
        newcols0 = {
                offset: 'X_an'+str(j+1)+str(1)
                }
        for k in np.arange(1, params['nshells'], 1):
            newcols = {
                    k+offset: 'X_an'+str(j+1)+str(k+1)
                    }
            newcols0.update(newcols)
        newcols2 = {
                0+offset+params['nshells']: 'X_elyte'+str(j+1),
                1+offset+params['nshells']: 'V_an'+str(j+1),
                2+offset+params['nshells']: 'V_dl'+str(j+1),
                }
        newcols0.update(newcols2)
        SV_df.rename(columns=newcols0, inplace=True)

    for j in np.arange(anode.npoints, sep.sep_max):
        offset = sep.offset_vec[j]
        newcols0 = {
                0+offset: 'X_elyte'+str(j+1),
                1+offset: 'V_elyte'+str(j+1)
                }
        newcols0.update(newcols0)
        SV_df.rename(columns=newcols0, inplace=True)

    for j in np.arange(sep.sep_max, sep.cat_max):
        offset = sep.offset_vec[j]
        newcols0 = {
                offset: 'X_cat'+str(j+1)+str(1)
                }
        for k in np.arange(1, params['nshells'], 1):
            newcols = {
                    k+offset: 'X_cat'+str(j+1)+str(k+1)
                    }
            newcols0.update(newcols)
        newcols2 = {
                0+offset+params['nshells']: 'X_elyte'+str(j+1),
                1+offset+params['nshells']: 'V_cat'+str(j+1),
                2+offset+params['nshells']: 'V_dl'+str(j+1),
                }
        newcols0.update(newcols2)
        SV_df.rename(columns=newcols0, inplace=True)

    newcols3 = {anode.npoints*anode.nVars + sep.npoints*sep.nVars + cathode.npoints*cathode.nVars: 'Time'}
    SV_df.rename(columns=newcols3, inplace=True)

    return SV_df
