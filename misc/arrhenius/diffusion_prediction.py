import sys
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from math import exp
from scipy import constants as const
from matplotlib import gridspec, pyplot as plt

def arrhenius(A, E_a, T):
    k = []
    R = const.physical_constants['Boltzmann constant in eV/K'][0]
    for t in T:
        k.append(A*exp(-E_a/R/t))
    return k

if __name__ == '__main__':

    working_path = '/home/jkusz/github/igloo/misc/arrhenius'
    os.chdir(working_path)

    A = 5e-11 #cm^2/s
    Ea_substitutional = 0.50 #eV
    Ea_interstial = 1.63 #eV

    temperature_data = list(range(30,506,5))
    inverse_temperature_data = []
    for t in temperature_data:
        inverse_temperature_data.append(1/t)
    substitutional_arrhenius_data = arrhenius(A, Ea_substitutional, temperature_data)
    interstitial_arrhenius_data = arrhenius(A,Ea_interstial,temperature_data)

    plt.clf()
    fig = plt.figure(figsize=(6,4))
    gs = gridspec.GridSpec(2,2)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3], sharex=ax3)
    gs.update(hspace=0.5)

    ##plot all lines
    ax1.plot(temperature_data,substitutional_arrhenius_data,'k-', label='substitutional')
    ax2.plot(temperature_data,interstitial_arrhenius_data,'r-', label='interstitial')
    ax3.plot(inverse_temperature_data,substitutional_arrhenius_data,'k--', label='substitutional')
    ax4.plot(inverse_temperature_data,interstitial_arrhenius_data,'r--', label='interstitial')
    plt.setp([ax1,ax2], xlabel=r'Temperature, T ($K$)')
    plt.setp([ax3,ax4], yscale='log', xlabel=r'Inverse Temperature, 1/T ($K^{-1}$)')
    plt.setp([ax1,ax2,ax3,ax4], ylabel=r'Diffusion, D ($cm^2/s$)')

    ax1.set_xlim(0,500)
    ax1.xaxis.set_major_locator(MultipleLocator(100))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(5)) # auto set minor ticks
    ax3.set_ylim(0,0.035)
    ax3.xaxis.set_major_locator(MultipleLocator(0.01))
    ax3.xaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks   

    #Set tick parameters
    # ax1.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', top='on') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    # ax1.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', top='on')
    # ax1.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
    # ax1.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

    gs.tight_layout(fig)
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    plt.style.use('seaborn-paper')
    plt.savefig('arrhenius_test_1.png',dpi=300,transparent=False,bbox_inches='tight')
    plt.show()