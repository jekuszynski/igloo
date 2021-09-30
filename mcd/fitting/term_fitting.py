from math import exp
import sys
import os
from lmfit.model import save_modelresult
import pandas as pd
import numpy as np
from matplotlib import gridspec, pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from scipy import optimize
from scipy.stats import linregress
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from lmfit.models import GaussianModel, ExponentialModel
from lmfit import Model
import random

def getWavelength(eV):
    return 1240/eV

def getEnergy(nm):
    return 1240/nm

def fwhm(fwhm):
    b = 4 * np.log(2) / (fwhm**2)
    return b

def getFWHM(b):
    return np.sqrt(4 * np.log(2) / b)

def getSigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))

def getWavenumber(eV):
    return eV * 8065.73

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

def plot_spectra(df,title,melt=True):
    plt.clf()
    fig,ax=plt.subplots(figsize=(8,4))

    ##setup color bar
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(-10,10,11),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    fig.colorbar(sm,ticks=range(-10,11,2),label='H (T)') #make color bar based on H (T) for plot

    ##convert dataframe into lengthwise plotting for using sns
    if melt == True:
        df_melt = df.melt(id_vars=['wavelength','energy'],var_name='field',value_name='deltaA')
        df_melt['field'] = df_melt['field'].str.split('_').str[0].astype(int)
    elif melt == False:
        df_melt = df
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
    ax.set_xlim(3.2,0.8)
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

    # print(df_melt)
    if melt == True:
        new_df = df_melt.pivot(index=['wavelength','energy'],columns='field')
    elif melt == False:
        new_df = df_melt.pivot(index=['wavelength','energy'],columns='field').reset_index()

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
    
    material = '5-1CFS'
    peak_count = '3'

    working_path = '/home/jkusz/github/igloo/mcd/fitting/temp/'
    os.chdir(working_path)

    # data_path = '/home/jkusz/github/strouse-data/for_publication/3-1CFS/adj_nir_pub/3-1CFS_worked_up_diff_mcd.csv'
    mcd_data_path = '/home/jkusz/github/igloo/mcd/fitting/temp/3-1CFS_merged_mcd_spectra.csv'
    abs_data_path = '/home/jkusz/github/strouse-data/abs/mcd/3-1_CFS_merged_abs_FOR_PUB.csv'

    '''Parse Data'''
    # data = parse_csv(data_path)
    mcd_data = pd.read_csv(mcd_data_path)
    abs_data = pd.read_csv(abs_data_path,na_values='--').dropna(how='all')

    '''Plot Data'''
    new_data = plot_spectra(mcd_data,'test_plot',melt=False) 
    # new_data.reset_index(level=(0,1),inplace=True) #change out of MultiIndex. Probably more powerful, but easier to work without for now. Use if melting prior.
    
    new_data = new_data.loc[new_data['wavelength'] <= 1600] #remove extraneous data for fitting
    # print(new_data)
    x = new_data['energy']
    y = new_data['deltaA',10]

    def ab_term_model1(x, ampA1, ampB1, cen1, wid1): #np.exp needed, because math.exp only expected single value insert, while np.exp can use array as is needed for lmfit.
        return (ampA1 * (x-cen1) * np.exp(-wid1*(x-cen1)**2)) + (ampB1 * np.exp(-wid1*(x-cen1)**2))

    def ab_term_model2(x, ampA2, ampB2, cen2, wid2):
        return (ampA2 * (x-cen2) * np.exp(-wid2*(x-cen2)**2)) + (ampB2 * np.exp(-wid2*(x-cen2)**2))
    
    def ab_term_model3(x, ampA3, ampB3, cen3, wid3):
        return (ampA3 * (x-cen3) * np.exp(-wid3*(x-cen3)**2)) + (ampB3 * np.exp(-wid3*(x-cen3)**2))

    def ab_term_model4(x, ampA4, ampB4, cen4, wid4):
        return (ampA4 * (x-cen4) * np.exp(-wid4*(x-cen4)**2)) + (ampB4 * np.exp(-wid4*(x-cen4)**2))

    # ab_term_model1 = ExpressionModel('ampA1 * (x-cen1) * exp(-wid1*(x-cen1)**2) + ampB1 * exp(-wid1*(x-cen1)**2)')
    # ab_term_model2 = ExpressionModel('ampA2 * (x-cen2) * exp(-wid2*(x-cen2)**2) + ampB2 * exp(-wid2*(x-cen2)**2)')
    # ab_term_model3 = ExpressionModel('ampA3 * (x-cen3) * exp(-wid3*(x-cen3)**2) + ampB3 * exp(-wid3*(x-cen3)**2)')

    # found from merge_mcd.py derivative.csv. 0 is nm, 1 is ev
    peak1 = [1060, 1240/1060]
    peak2 = [] # allow to float to find a/b term. Leave open to quantify what behavior is here. 
    peak3 = [482, 1240/482]
    peak4 = [350, 1240/350]
    
    def monte_carlo(n=10):
        ab_term_model_total = Model(ab_term_model1) + Model(ab_term_model2) + Model(ab_term_model3)
        pars = ab_term_model_total.make_params()
        param_list = []
        for i in range (0,n):
            valampA1 = random.uniform(-10, -0.002)
            # valampB1 = random.uniform(-20,20)
            # valcen1 = random.uniform(-10,10)
            valwid1 = random.uniform(fwhm(0.9),fwhm(0.2))
            # valampA2 = random.uniform(-np.inf, -0.001)
            valampB2 = random.uniform(-10,10)
            valcen2 = random.uniform(1.2,1.9)
            valwid2 = random.uniform(fwhm(1),fwhm(0.1))
            valampA3 = random.uniform(-10,10)
            valampB3 = random.uniform(-10,10)
            # valcen3 = random.uniform(-10,10)
            valwid3 = random.uniform(fwhm(2.4),fwhm(0.2))
            # valampA4 = random.uniform(-np.inf, -0.001)
            # valampB4 = random.uniform(-10,10)
            # valcen4 = random.uniform(-10,10)
            # valwid4 = random.uniform(-10,10)

            pars['ampA1'].set(value=valampA1, max=-0.001)
            pars['ampB1'].set(value=0, vary=False)
            pars['cen1'].set(value=peak1[1], vary=False)
            pars['wid1'].set(value=valwid1, min=fwhm(1), max=fwhm(0.1))
            pars['ampA2'].set(value=0, vary=False)
            pars['ampB2'].set(value=valampB2, vary=True)
            pars['cen2'].set(value=valcen2, vary=True, min=1.3, max=2.0)
            pars['wid2'].set(value=valwid2, min=fwhm(1), max=fwhm(0.1))
            pars['ampA3'].set(value=valampA3)
            pars['ampB3'].set(value=valampB3)
            pars['cen3'].set(value=peak3[1], vary=False)
            pars['wid3'].set(value=valwid3, min=fwhm(2.5), max=fwhm(0.1))
            # pars['ampA4'].set(value=2)
            # pars['ampB4'].set(value=-4, vary=True, min=-8, max=4)
            # pars['cen4'].set(value=peak4[1], vary=True, min=2.8, max=5)
            # pars['wid4'].set(fwhm(1.5), min=fwhm(2), max=fwhm(0.5))
            
            # print('Parameter    Value       Stderr')
            fit = ab_term_model_total.fit(y, pars, x=x)

            if fit.redchi < 0.05:
                print('redchi < 0.05 found!: ' + str(fit.redchi))
                peak_list = []
                ampA_list = []
                ampB_list = []
                wid_list = []
                cen_list = []
                for name, param in fit.params.items():
                    # print('{:7s} {:11.5f} {:11.5f}'.format(name, param.value, param.stderr))
                    if 'ampA' in name:
                        ampA_list.append([param.value, param.stderr])
                    elif 'ampB' in name:
                        ampB_list.append([param.value, param.stderr])
                    elif 'cen' in name:
                        cen_list.append([param.value, param.stderr])
                    elif 'wid' in name:
                        wid_list.append([param.value, param.stderr])
                for ampA, ampB, wid, cen in zip(ampA_list, ampB_list, wid_list, cen_list):
                    param_list.append([i,fit.redchi,ampA[0],ampA[1],ampB[0],ampB[1],cen[0],cen[1],wid[0],wid[1]])
            else:
                print('reduced chi not < 1: ' + str(fit.redchi))
            
        parameter_testing_df = pd.DataFrame(param_list, columns = ['iteration','redchi','ampA','ampA_stderr','ampB','ampB_stderr','cen','cen_stderr','wid','wid_stderr'])
        print(parameter_testing_df)
        parameter_testing_df.to_csv(working_path + 'parameter_testing.csv')

    monte_carlo(2)
    sys.exit()


    #increasing peak number is an increase in ev, decrease in nm.
    ab_term_model_total = Model(ab_term_model1) + Model(ab_term_model2) + Model(ab_term_model3)
    pars = ab_term_model_total.make_params()
    pars['ampA1'].set(value=-2, max=-0.001)
    pars['ampB1'].set(value=0, vary=False)
    pars['cen1'].set(value=peak1[1], vary=False)
    pars['wid1'].set(fwhm(0.6), min=fwhm(0.8), max=fwhm(0.2))
    pars['ampA2'].set(value=0, vary=False)
    pars['ampB2'].set(value=-0.4, vary=True)
    pars['cen2'].set(value=1.8, vary=True, min=1.3, max=2.0)
    pars['wid2'].set(fwhm(0.7), min=fwhm(0.8), max=fwhm(0.2))
    pars['ampA3'].set(value=0)
    pars['ampB3'].set(value=-2, vary=True)
    pars['cen3'].set(value=peak3[1], vary=False)
    pars['wid3'].set(fwhm(0.5), min=fwhm(1.5), max=fwhm(0.1))
    # pars['ampA4'].set(value=2)
    # pars['ampB4'].set(value=-4, vary=True, min=-8, max=4)
    # pars['cen4'].set(value=peak4[1], vary=True, min=2.8, max=5)
    # pars['wid4'].set(fwhm(1.5), min=fwhm(2), max=fwhm(0.5))

    print(pars)

    # init_fit = ab_term_model_total.eval(pars, x=x)
    fitting_result = ab_term_model_total.fit(y, pars, x=x)
    save_modelresult(fitting_result, material + '_' + peak_count + 'peak_params' + '.sav')
    print(fitting_result.fit_report(min_correl=0.5))

    wid_list = []
    cen_list = []
    param_list = []

    # print('-------------------------------')
    # print('Parameter    Value       Stderr')
    for name, param in fitting_result.params.items():
        # print('{:7s} {:11.5f} {:11.5f}'.format(name, param.value, param.stderr))
        if 'cen' in name:
            cen_list.append(param.value)
        elif 'wid' in name:
            wid_list.append(getSigma(getFWHM(param.value)))

    for wid, cen in zip(wid_list,cen_list):
        param_list.append([cen,wid])

    mcd_comps = fitting_result.eval_components(x=x)

    abs_x = abs_data['energy']
    abs_y = abs_data['normalized_absorbance']

    gauss1 = GaussianModel(prefix='g1_')
    gauss2 = GaussianModel(prefix='g2_')
    gauss3 = GaussianModel(prefix='g3_')
    gauss4 = GaussianModel(prefix='g4_')
    # exp = ExponentialModel(prefix='exp_')

    abs_model = gauss1 + gauss2 + gauss3 + gauss4
    pars = abs_model.make_params()

    pars['g1_center'].set(value=param_list[0][0], vary=False)
    pars['g1_sigma'].set(value=param_list[0][1], vary=False)
    pars['g1_amplitude'].set(value=.409, min=0, vary=True)
    pars['g2_center'].set(value=param_list[1][0], vary=False)
    pars['g2_sigma'].set(value=param_list[1][1], vary=False)
    pars['g2_amplitude'].set(value=.3668, min=0, vary=True)
    pars['g3_center'].set(value=param_list[2][0], vary=False)
    pars['g3_sigma'].set(value=param_list[2][1], vary=False)
    pars['g3_amplitude'].set(value=.2, min=0)
    pars['g4_center'].set(value=4, min=2.5, max=6, vary=True)
    pars['g4_sigma'].set(value=fwhm(1), min=fwhm(3), max=fwhm(.2), vary=True)
    pars['g4_amplitude'].set(value=.5, min=0)
    # pars['exp_amplitude'].set(value=0.01)
    # pars['exp_decay'].set(value=-5, max=0)
    
    # init = mod.eval(pars, x=x)
    abs_fit_output = abs_model.fit(abs_y, pars, x=abs_x)

    print(abs_fit_output.fit_report())
    # abs_fit_output.plot()

    abs_comps = abs_fit_output.eval_components(x=abs_x)

    '''Fit all data'''
    plt.clf()
    fig = plt.figure(figsize=(6,4))
    gs = gridspec.GridSpec(2,1, height_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    ax1.plot([0,2000],[0,0],'k-',lw=0.8)
    ax1.plot(x, y, 'k-', label='raw data')
    # ax1.plot(x, fitting_result.init_fit, 'k--', label='initial fit')

    colors = ["#ff3c38","#fff275","#e76f51","#6699cc"]

    ax1.plot(x, fitting_result.best_fit, 'r-', label='best fit')
    ax1.fill_between(x, list(mcd_comps.values())[0], color=colors[-1], lw=1.5, ls='--', label='lspr', alpha=0.2)
    ax1.fill_between(x, list(mcd_comps.values())[1], color=colors[-2], lw=1.5, ls='--', label='B-term-mixing', alpha=0.2)
    ax1.fill_between(x, list(mcd_comps.values())[2], color=colors[-3], lw=1.5, ls='--', label='IB?', alpha=0.2)
    # ax1.fill_between(x, list(mcd_comps.values())[3], color=colors[-4], lw=1.5, ls='--', label='UV', alpha=0.2)

    ax2.plot(abs_x, abs_y, 'k-', label='raw data')
    ax2.plot(abs_x, abs_fit_output.best_fit, 'r-', label='best fit')
    ax2.fill_between(abs_x, list(abs_comps.values())[0], color=colors[-1], lw=1.5, ls='--', label='lspr', alpha=0.2)
    ax2.fill_between(abs_x, list(abs_comps.values())[1], color=colors[-2], lw=1.5, ls='--', label='B-term-mixing', alpha=0.2)
    ax2.fill_between(abs_x, list(abs_comps.values())[2], color=colors[-3], lw=1.5, ls='--', label='IB?', alpha=0.2)
    ax2.fill_between(abs_x, list(abs_comps.values())[3], color=colors[-4], lw=1.5, ls='--', label='UV', alpha=0.2)
    # ax2.plot(abs_x, list(abs_comps.values())[3], color='purple', lw=1.5, ls='--', label='Continuum')

    ax1.legend(loc=0)
    ax2.set_xlabel(r'Wavelength, $\lambda$ (nm)')
    ax1.set_ylabel(r'MCD, $\Delta$A/A$_{\mathrm{max}}$ (x $10^{-3}$)',labelpad=6, size=9)
    ax2.set_ylabel(r'Absorbance, A (a.u.)',labelpad=14,size=9)

    ax1.set_xlim(3.09, 0.75)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(0,1.1)
    ax1.set_ylim(-1.4,0.5)

    ax3 = twin_axis(ax1)
    ax4 = twin_axis(ax2)

    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax4.xaxis.set_major_formatter(plt.NullFormatter())
    ax3.set_xlabel(r'Energy, $h\nu$ (eV)')
    ax4.set_xlabel('')
    
    ax1.xaxis.set_major_locator(MultipleLocator(0.5))
    ax2.xaxis.set_major_locator(MultipleLocator(0.5))
    ax3.xaxis.set_major_locator(MultipleLocator(500))
    ax4.xaxis.set_major_locator(MultipleLocator(500))
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
    plt.savefig(material + '_CFS_' + peak_count + 'peak_fit_zoomed.png', format='png',dpi=300)

    sys.exit() 

    # '''Save worked up data in CSV'''
    # dict_to_csv(data,'worked_up_data')