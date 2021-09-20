from math import exp
import sys
import os
from numpy.ma import compress_cols
import pandas as pd
import numpy as np
from matplotlib import gridspec, pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from scipy import optimize
from scipy.stats import linregress
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from lmfit.models import ExpressionModel

def getWavelength(eV):
    return 1240/eV

def getEnergy(nm):
    return 1240/nm

def fwhm(fwhm):
    b = 4 * np.log(2) / fwhm**2
    return b

def twin_axis(ax,xlabel=r'Wavelength (nm)'):
    new_ax = ax.twiny()
    new_ax.set_xlabel(xlabel)
    new_ax.set_xscale('function',functions=(getWavelength,getEnergy)) # set twin scale (convert degree eV to nm)
    xmin, xmax = ax.get_xlim() # get left axis limits
    new_ax.set_xlim((getWavelength(xmax),getWavelength(xmin))) # apply function and set transformed values to right axis limits
    new_ax.xaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks
    new_ax.plot([],[]) # set an invisible artist to twin axes to prevent falling back to initial values on rescale events
    return new_ax    

def parse_csv(path):
    # use_mdeg = np.r_[0:2,#:#] #not sure if I'll ever want to use this
    use_deltaA = np.r_[0:2,13:24]
    df = pd.read_csv(path,usecols=use_deltaA)
    return df

def plot_spectra(df,title):
    plt.clf()
    fig,ax=plt.subplots(figsize=(8,4))

    ##setup color bar
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(-10,10,11),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    fig.colorbar(sm,ticks=range(-10,11,2),label='H (T)') #make color bar based on H (T) for plot

    ##convert dataframe into lengthwise plotting for using sns
    df_melt = df.melt(id_vars=['wavelength','energy'],var_name='field',value_name='deltaA')
    df_melt['field'] = df_melt['field'].str.split('_').str[0].astype(int)
    # print(df_melt)
    ##plot all lines
    sns.lineplot(data=df_melt,x='energy',y='deltaA', linewidth=0.6,
                hue='field',hue_norm=(-10,10),
                palette=sns.color_palette('coolwarm_r',as_cmap=True),
                legend=None)

    ax.set_xlabel('Energy (eV)')
    ax.plot([-10,10],[0,0],color='black',linestyle='-',linewidth='1') #add 0T baseline

    ax.set_ylabel(r'$\Delta$A/A$_{\mathrm{max}}$ (x $10^{-3}$)')
    # ax.set_ylabel('MCD (mdeg)')
    ax.set_xlim(1.6,0.8)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks
    # ax.set_ylim(-1.5,1.5)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks
    
    ax2 = ax.twiny() # creates a new axis with invisible y and independent x on opposite side of first x-axis
    ax2.set_xlabel(r'Wavelength (nm)')
    ax2.set_xscale('function',functions=(getWavelength,getEnergy)) # set twin scale (convert degree eV to nm)
    xmin, xmax = ax.get_xlim() # get left axis limits
    ax2.set_xlim((getWavelength(xmax),getWavelength(xmin))) # apply function and set transformed values to right axis limits
    ax2.xaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks
    
    ax2.plot([],[]) # set an invisible artist to twin axes to prevent falling back to initial values on rescale events

    #Set tick parameters    
    ax.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in')
    ax.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

    ax2.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in')
    ax2.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in')

    # plt.style.use('classic') #Doesn't seem to affect plot.
    plt.tight_layout()

    plt.style.use('seaborn-paper')
    plt.savefig('mcd_fit' + title,dpi=300,transparent=False,bbox_inches='tight')
    plt.show()

    new_df = df_melt.pivot(index=['wavelength','energy'],columns='field')
    return new_df

def _mcdtermsAB(x, cen, ampA, ampB, wid):
    return ampA*(x-cen)*exp(-wid*(x-cen)**2) +\
                ampB*exp(-wid*(x-cen)**2)

def aterm(dif, fwhm, abs_peak, field):
    return (exp(1/2)*dif*fwhm)/(2*2.35*0.4671*abs_peak*field)
    
def bterm(mcd_peak, abs_peak, field):
    return mcd_peak/(0.4671*abs_peak*field)

def plot_mcdterms(df,fname,num=1,xlabel='wavelength',ylabel='deltaA'):
    fig = plt.figure(figsize=(4,4))
    gs = gridspec.GridSpec(2,1, height_ratios=[1,0.25])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    gs.update(hspace=0)
    p0 = [0, 0, 0] #amplitude, center, width; adjust as needed
    bounds = []
    field_list=[]
    width1_list=[]
    width1_err_list=[]
    for field in df:
        xdata = df[xlabel]
        ydata = df[ylabel]
        ax1.plot(xdata, ydata, "ko", ms=0.5)
        # use one peak
        if num==1:
            popt_mcdterms, pcov_mcdterms = optimize.curve_fit(_mcdtermsAB, xdata, ydata, p0=p0)
            perr_mcdterms = np.sqrt(np.diag(pcov_mcdterms))
            pars_1 = popt_mcdterms[0:3]
            mcd_term_peak_1 = _mcdtermsAB(xdata, *pars_1)
            ax1.plot(xdata, _mcdtermsAB(xdata, *popt_mcdterms), 'r--', lw=0.4)
            residual_mcd_terms = ydata - (_mcdtermsAB(xdata, *popt_mcdterms))
            ax2.plot(xdata, residual_mcd_terms, "bo", ms=0.4)
            ax2.fill_between(xdata,0,residual_mcd_terms,facecolor='blue',alpha=0.1)
            field_list.append(field)
            width1_list.append(pars_1[2])
            width1_err_list.append(perr_mcdterms[2])
            peak_params_df = pd.DataFrame(list(zip(field_list, width1_list, width1_err_list)),columns=['Pressure','FWHM','err'])

    # ax1.set_xlim(-10000,10000)
    # ax2.set_xlim(-10000,10000)
    # ax2.set_ylim(-25,25)

    ax2.set_xlabel(r'Frequency, $\nu-\nu_0$ (MHz)')
    ax1.set_ylabel(r'Signal, $\Delta$I$_{pr}$ (mV)',labelpad=10)
    ax2.set_ylabel("Residual",labelpad=4.5)

    ax1.xaxis.set_major_formatter(plt.NullFormatter())

    ax1.xaxis.set_major_locator(MultipleLocator(5000))
    ax2.xaxis.set_major_locator(MultipleLocator(5000))
    # ax2.yaxis.set_major_locator(MultipleLocator(10))

    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())

    ax1.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom='off') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax1.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom='off')
    ax1.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on', bottom='off')
    ax1.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on', bottom='off')

    ax2.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', top='off') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax2.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', top='off')
    ax2.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on', top='off')
    ax2.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on', top='off')

    fig.tight_layout()
    fig.savefig(fname+'.png', format="png",dpi=1000)

    return peak_params_df

def plot_peak_params(df,fname,num=1):
    df.sort_values('Pressure',inplace=True)
    fig = plt.figure(figsize=(3,4))
    gs = gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0])

    if num==1:
        xdata = df['Pressure']
        ydata = df['FWHM']  
        yerr = df['err']
        ax1.errorbar(xdata, ydata, yerr=yerr, ms=8, fmt='.', color='maroon')
        m, b, r, p, err = linregress(xdata,ydata)[:5]
        poly1d_fn = np.poly1d([m,b])
        ax1.plot(xdata, poly1d_fn(xdata), 'k--', lw=1.2, c='maroon')
        ax1.text(5, 4000, 'y=%0.3fx+%0.3f\n' % (m, b) + r'$R^2$=%0.5f' % (r) + '\n' + r'$m\pm$%0.4f' % (err),color='maroon')

    ax1.set_xlim(-5,105)
    ax1.set_ylim(0,5000)
    
    ax1.set_xlabel(r'Pressure, (Torr)')
    ax1.set_ylabel(r'Frequency, $\Delta\nu$ (MHz)')

    ax1.xaxis.set_major_locator(MultipleLocator(20))
    ax1.yaxis.set_major_locator(MultipleLocator(1000))
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())

    ax1.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax1.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in')
    ax1.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
    ax1.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

    fig.tight_layout()
    fig.savefig(fname+'.png', format="png",dpi=1000)

if __name__ == '__main__':
    
    working_path = '/home/jkusz/github/igloo/mcd/fitting/temp/'
    os.chdir(working_path)

    data_path = '/home/jkusz/github/strouse-data/for_publication/3-1CFS/adj_nir_pub/3-1CFS_worked_up_diff_mcd.csv'

    '''Parse Data'''
    data = parse_csv(data_path)

    '''Plot Data'''
    new_data = plot_spectra(data,'test_plot') 
    new_data.reset_index(level=(0,1),inplace=True) #change out of MultiIndex. Probably more powerful, but easier to work without for now.

    x = new_data['wavelength']
    y = new_data['deltaA',10]

    # a_term_model = ExpressionModel('ampA * (x-cen) * exp(-wid*(x-cen)**2)')
    # b_term_model = ExpressionModel('ampB * exp(-wid*(x-cen)**2)')

    ab_term_model1 = ExpressionModel('ampA1 * (x-cen1) * exp(-wid1*(x-cen1)**2) + ampB1 * exp(-wid1*(x-cen1)**2)')
    ab_term_model2 = ExpressionModel('ampA2 * (x-cen2) * exp(-wid2*(x-cen2)**2) + ampB2 * exp(-wid2*(x-cen2)**2)')
    ab_term_model3 = ExpressionModel('ampA3 * (x-cen3) * exp(-wid3*(x-cen3)**2) + ampB3 * exp(-wid3*(x-cen3)**2)')
    ab_term_model4 = ExpressionModel('ampA4 * (x-cen4) * exp(-wid4*(x-cen4)**2) + ampB2 * exp(-wid4*(x-cen4)**2)')

    pars = ab_term_model1.make_params()
    pars['ampA1'].set(value=0.002, min=0, max=2)
    pars['ampB1'].set(value=0, vary=False)
    pars['cen1'].set(value=1150,vary=False)
    pars['wid1'].set(fwhm(500))

    pars.update(ab_term_model2.make_params())
    pars['ampA2'].set(value=0.002, min=0, max=2)
    pars['ampB2'].set(value=0, vary=False)
    pars['cen2'].set(value=1100,vary=False)
    pars['wid2'].set(fwhm(500))

    ab_term_model_total = ab_term_model1+ab_term_model2

    init_fit = ab_term_model_total.eval(pars, x=x)
    fitting_result = ab_term_model_total.fit(y, pars, x=x)

    print(fitting_result.fit_report(min_correl=0.5))

    comps = fitting_result.eval_components(x=x)
    print(comps)

    plt.clf()
    fig = plt.figure(figsize=(6,4))
    gs = gridspec.GridSpec(2,1, height_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    ax1.plot([0,2000],[0,0],'k-',lw=0.8)
    ax1.plot(x, y, 'k-', label='raw data')
    # ax1.plot(x, fitting_result.init_fit, 'k--', label='initial fit')
    ax1.plot(x, fitting_result.best_fit, 'r-', label='best fit')
    ax1.plot(x, comps[0],'g--',label='ab_1')
    ax1.plot(x, comps[1],'m--',label='ab_2')
    ax1.legend(loc='best')
    ax2.set_xlabel(r'Wavelength, $\lambda$ (nm)')
    ax1.set_ylabel(r'MCD, $\Delta$A/A$_{\mathrm{max}}$ (x $10^{-3}$)',labelpad=6, size=9)
    ax2.set_ylabel(r'Absorbance, A (a.u.)',labelpad=14,size=9)

    ax1.set_xlim(400, 1600)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(0.1,1.2)
    ax1.set_ylim(-1.4,1.4)

    ax3 = twin_axis(ax1)
    ax4 = twin_axis(ax2)

    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax4.xaxis.set_major_formatter(plt.NullFormatter())
    ax3.set_xlabel(r'Energy, $h\nu$ (eV)')
    ax4.set_xlabel('')
    
    ax1.xaxis.set_major_locator(MultipleLocator(200))
    ax2.xaxis.set_major_locator(MultipleLocator(200))
    # ax2.yaxis.set_major_locator(MultipleLocator(10))

    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax4.xaxis.set_minor_locator(AutoMinorLocator())

    ax1.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom='off') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax1.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom='off')
    ax1.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on', bottom='off')
    ax1.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on', bottom='off')

    ax2.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom=True) #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax2.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom=True)
    ax2.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right=True, top=False)
    ax2.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right=True, top=False)

    ax3.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom=False, top=True) 
    ax3.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom=False, top=True)
    ax4.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom=False, top=True) 
    ax4.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom=False, top=True)


    gs.tight_layout(fig)
    gs.update(hspace=0)
    plt.show()
    plt.savefig('fitting.png', format='png',dpi=300)

    sys.exit() 

    # '''Save worked up data in CSV'''
    # dict_to_csv(data,'worked_up_data')