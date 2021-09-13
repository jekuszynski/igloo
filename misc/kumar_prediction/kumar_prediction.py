import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def plot_spectra(df,title,xdata,ydata):
    plt.clf()
    fig,ax=plt.subplots(figsize=(5,3))

    ##plot all lines
    ax.scatter(df[xdata],df[ydata],c='black')
    m, b, r = linregress(df[xdata],df[ydata])[:3]

    if xdata == 'Mobility':
        ax.set_xlabel(r'Mobility, $\mu$ ($cm^2 V^{-1} s^{-1}$)')
        new_x = 0.16
        ax.set_xlim(0,3.0)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5)) # auto set minor ticks
        ax.set_ylim(1.5,2.8)
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2)) # auto set minor ticks
        ax.text(2.3, 2.6, 'y=%0.3fx+%0.3f\n' % (m, b) + r'$R^2$=%0.5f' % (r),color='black',size=8)

    elif xdata == 'Resistivity':
        ax.set_xlabel(r'Resistivity, $\rho$ ($m\Omega cm$)')
        new_x = 21.4
        ax.set_xlim(5.0,23.0)
        ax.xaxis.set_major_locator(MultipleLocator(5.0))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5)) # auto set minor ticks
        ax.set_ylim(1.25,6.25)
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5)) # auto set minor ticks
        ax.text(6, 5.5, 'y=%0.3fx+%0.3f\n' % (m, b) + r'$R^2$=%0.5f' % (r),color='black',size=8)

    elif xdata == 'R&N':
        ax.set_xlabel(r'$\rho*n$ ($m\Omega$' + ' ' + r'$cm^{-2}$)')
        new_x = 21.4 * 1.9e21
        ax.set_xlim(0.1e21,4.3e22)
        ax.xaxis.set_major_locator(MultipleLocator(0.5e22))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5)) # auto set minor ticks
        ax.set_ylim(1.4,5.0)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2)) # auto set minor ticks
        # ax.text(5e21, 4.0, 'y=%0.3fx+%0.3f\n' % (m, b) + r'$R^2$=%0.5f' % (r),color='black',size=8)

    new_y = new_x*m+b
    print('\n' + xdata + r' gives $m^{*}$ of %0.3f' % (new_y))

    new_xdata = df[xdata].values.tolist()
    new_ydata = df[ydata].values.tolist()
    new_xdata.insert(0,new_x)
    new_ydata.insert(0,new_y)

    poly1d_fn = np.poly1d([m,b])
    ax.scatter(new_x,new_y,c='blue',zorder=4)
    ax.plot(new_xdata, poly1d_fn(new_xdata), '--', lw=1.2, c='black')
    
    ax.set_ylabel(r'Effective Mass, $m^{*}$/$m_e$')

    #Set tick parameters
    ax.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', top='on') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

    # plt.tight_layout()

    plt.style.use('seaborn-paper')
    plt.savefig(title,dpi=300,transparent=False,bbox_inches='tight')
    plt.show()

if __name__ == '__main__':

    working_path = '/home/jkusz/github/igloo/misc'
    os.chdir(working_path)

    kumar_data = [[2.56,0.36,9.3,1.35e21],[2.42,1.15,7.2,6.24e20],[1.62,2.80,6.0,3.65e20]]
    df = pd.DataFrame(kumar_data,columns=['Effective Mass','Mobility','Resistivity','Density'])
    df['R&N'] = df['Resistivity'] * df['Density']
    # print(df['R&N'])

    ## For some reason these wont work in tandem... run one at a time
    # plot_spectra(df,'u_vs_m','Mobility','Effective Mass')
    # plot_spectra(df,'r_vs_m','Resistivity','Effective Mass')
    plot_spectra(df, 'rn_vs_m','R&N','Effective Mass')