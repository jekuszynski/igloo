import sys
import os
import pandas as pd
import numpy as np
from matplotlib import gridspec, pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import savgol_filter
from lmfit.models import GaussianModel, VoigtModel


def getWavelength(eV):
    return 1240/eV

def getEnergy(nm):
    return 1240/nm

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

def plot_spectra(df,abs_df,title):
    fig = plt.figure(figsize=(6,4))
    gs = gridspec.GridSpec(2,1, height_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    ### Setup Data for Abs top fig
    x = abs_df['wavelength']
    y = abs_df['normalized_absorbance']
    ax1.plot(x, y, 'k-', label='abs data')

    ### Setup Data for MCD bottom fig
    ##setup color bar
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(-10,10,11),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    cbaxes = inset_axes(ax2, width="30%", height="4%", loc='upper right',borderpad=1.2) 
    # cbaxes.xaxis.set_ticks_position('top') # allegedly this should work, but doesn't. Instead using borderpad to keep away from edges.
    cb = plt.colorbar(sm, ticks=range(-10,11,10), cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    sns.lineplot(data=df,x='wavelength',y='deltaA', linewidth=0.6,
                hue='field', hue_norm=(-10,10),
                palette=sns.color_palette('coolwarm_r',as_cmap=True),
                legend=None, ax=ax2)
    ax2.plot([0,2000],[0,0],'k-',lw=0.8)

    # ax1.legend(loc='best')
    ax2.set_xlabel(r'Wavelength, $\lambda$ (nm)')
    ax2.set_ylabel(r'MCD, $\Delta$A/A$_{\mathrm{max}}$ (x $10^{-3}$)',labelpad=5,size=9)
    ax1.set_ylabel(r'Absorbance, A (a.u.)',labelpad=14,size=9)

    ax1.set_xlim(400,1600)
    ax1.set_ylim(0.1,1.2)
    ax2.set_xlim(400,1600)
    ax2.set_ylim(-1.4,1.4)

    ax3 = twin_axis(ax1,xlabel=r'Energy, $h\nu$ (eV)')
    ax4 = twin_axis(ax2,xlabel=r'')

    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax4.xaxis.set_major_formatter(plt.NullFormatter())

    ax1.xaxis.set_major_locator(MultipleLocator(200))
    ax2.xaxis.set_major_locator(MultipleLocator(200))
    ax3.xaxis.set_major_locator(MultipleLocator(0.5))
    ax4.xaxis.set_major_locator(MultipleLocator(0.5))

    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())

    ax1.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom='on', top=False) #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
    ax1.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom='on', top=False)
    ax1.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
    ax1.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

    ax2.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom='on', top=False) 
    ax2.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom='on', top=False)
    ax2.yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
    ax2.yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

    ax3.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom=False) 
    ax3.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom=False)

    ax4.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', bottom=False)
    ax4.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', bottom=False)

    gs.tight_layout(fig)
    gs.update(hspace=0)
    plt.style.use('seaborn-paper')
    plt.savefig(title,format='png',dpi=300)
    plt.show()

def convert_columns_to_field(df):
    df_melt = df.melt(id_vars=['wavelength','energy'],var_name='field',value_name='deltaA')
    df_melt['field'] = df_melt['field'].str.split('_').str[0].astype(int)
    new_df = df_melt.pivot(index=['wavelength','energy'],columns='field')
    return new_df

def merge_spectra(vis_data,nir_data,wavelength_match):
    # rename field titles to actual field values
    vis_field_data = convert_columns_to_field(vis_data)
    nir_field_data = convert_columns_to_field(nir_data)

    # change out of MultiIndex. Probably more powerful, but brain too smol to understand
    vis_field_data.reset_index(level=(0,1),inplace=True) 
    nir_field_data.reset_index(level=(0,1),inplace=True)

    # print(vis_field_data)
    # sys.exit()

    # use ### nm as conversion of NIR -> VIS data.
    for field in vis_field_data:
        if 'deltaA' in field[0]:
            vis_ref_point = float(vis_field_data.loc[vis_field_data[('wavelength','')] == wavelength_match, field])        
            nir_ref_point = float(nir_field_data.loc[nir_field_data[('wavelength','')] == wavelength_match, field])
            try:
                nir_to_vis_conversion = vis_ref_point/nir_ref_point
            except ZeroDivisionError: 
                nir_to_vis_conversion = 0
            vis_field_data[field] = vis_field_data.loc[vis_field_data[('wavelength','')] <= wavelength_match, field] # feel free to comment out this line or next to add "error" gradient to merged regions.
            nir_field_data[field] = nir_field_data.loc[nir_field_data[('wavelength','')] > wavelength_match, field]
            nir_field_data[field] = nir_field_data[field] * nir_to_vis_conversion
    merged_spectra = vis_field_data.merge(nir_field_data,how='outer')

    merged_melt = merged_spectra.melt(id_vars=['wavelength','energy'],value_name='deltaA')
    merged_melt = merged_melt.loc[:,merged_melt.columns.notnull()]

    return merged_melt

### Use to test SSR for fits
# win_list = list(range(5,252,4))
# poly_list = list(range(2,10))
# testing_array = np.empty((0,3))

# for window_length in win_list:
#     print("Testing window length of: " + window_length)
#     for polyorder in poly_list:
#         if window_length > polyorder:
#             line_fit = savgol_filter(y, window_length=window_length, polyorder=polyorder, deriv=0)
#             lsr = 0
#             for new_val, orig_val in zip(line_fit, y):
#                 lsr += np.sqrt((orig_val - new_val)**2) + lsr
#                 testing_array = np.append(testing_array,np.array([[window_length,polyorder,lsr]]),axis=0)
#         else:
#             pass
# # print(testing_array)
# lsr_df = pd.DataFrame(testing_array, columns=['window_length','polyorder','lsr'])
# best_fit_lsr = lsr_df['lsr'].idxmin()
# print(lsr_df.iloc[best_fit_lsr])

if __name__ == '__main__':
    working_path = '/home/jkusz/github/igloo/mcd/fitting/temp/'
    os.chdir(working_path)

    '''Setup Data Paths'''
    vis_data_path = '/home/jkusz/github/strouse-data/for_publication/3-1CFS/vis_pub/3-1CFS_worked_up_diff_mcd.csv'
    nir_data_path = '/home/jkusz/github/strouse-data/for_publication/3-1CFS/adj_nir_pub/3-1CFS_worked_up_diff_mcd.csv'
    abs_data_path = '/home/jkusz/github/strouse-data/abs/mcd/3-1_CFS_merged_abs_FOR_PUB.csv'

    '''Parse Data'''
    vis_mcd_data = parse_csv(vis_data_path)
    nir_mcd_data = parse_csv(nir_data_path)
    absorption_spectra_data = pd.read_csv(abs_data_path,na_values='--').dropna(how='all')
    print(absorption_spectra_data)

    '''Merge and Scale Data'''
    merged_mcd_spectra = merge_spectra(vis_mcd_data,nir_mcd_data,900)
    merged_mcd_spectra.to_csv('3-1CFS_merged_mcd_spectra.csv')
    
    '''Plot Abs with Merged MCD'''
    plot_spectra(merged_mcd_spectra,absorption_spectra_data,'merged_spectra.png')

    '''Find 2nd Der of Abs'''
    x = absorption_spectra_data['wavelength']
    y = absorption_spectra_data['normalized_absorbance']
    x_ev = absorption_spectra_data['energy']
    x_set = [x, x_ev]

    fit = savgol_filter(y, window_length=71, polyorder=3, deriv=0)
    dy  = savgol_filter(y, window_length=71, polyorder=3, deriv=1)
    ddy = savgol_filter(y, window_length=71, polyorder=3, deriv=2)
    ddyy= savgol_filter(y, window_length=15, polyorder=2, deriv=2)

    for num, x in enumerate(x_set):
        plt.clf()
        fig = plt.figure(figsize=(6,4))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0])
        ax2 = ax1.twinx()
        line1 = ax1.plot(x, y, 'ko', ms=1, label='raw data')
        # x_range = np.linspace(x[0],x[-1],1000)
        line2 = ax1.plot(x, fit, 'k-', lw=1, label='data fit')
        # plt.plot(x,first_der_fit, 'bo',label='1st der')
        # plt.plot(x,second_der_fit, '--', label='2nd der')
        # plt.plot(x,dy,label='1st der')
        line3 = ax2.plot(x, ddy, 'r-', lw=1, label='2nd der')
        line4 = ax2.plot(x, ddyy, 'r--', lw=0.5, label='raw 2nd der')
        all_lines = line1+line2+line3+line4
        line_labels = [l.get_label() for l in all_lines]
        ax1.legend(all_lines, line_labels, loc='best')
        ax1.set_ylabel(r'Absorbance, A (a.u.)')
        ax2.tick_params(axis='y',colors='red')
        ax2.set_ylabel(r'2$^{nd}$ Derivative, $\frac{\partial ^{2}A}{\partial x^{2}}$', color='red')
        if num == 0:
            ax1.set_xlim(400,1700)
        if num == 1:
            ax1.set_xlim(0.73,3.08)
        gs.tight_layout(fig)
        plt.show()
        plt.savefig('2nd_deriv_abs_' + str(num) + '.png',dpi=300)

    '''Find Peaks onto Abs'''

    # a_term_model = ExpressionModel('ampA * (x-cen) * exp(-wid*(x-cen)**2)')
    # b_term_model = ExpressionModel('ampB * exp(-wid*(x-cen)**2)')

    gauss1 = GaussianModel(prefix='g1_')
    pars = gauss1.make_params()

    pars['g1_center'].set(value=1.0, min=0.8, max=2)
    pars['g1_sigma'].set(value=0.2, min=0.05, max=2)
    pars['g1_amplitude'].set(value=.5, min=.2, max=2)

    gauss2 = GaussianModel(prefix='g2_')
    pars.update(gauss2.make_params())

    pars['g2_center'].set(value=1.5, min=0.8, max=1.6)
    pars['g2_sigma'].set(value=0.2, min=0.05, max=2)
    pars['g2_amplitude'].set(value=.5, min=.2, max=2)

    mod = gauss1 + gauss2
    
    init = mod.eval(pars, x=x)
    output = mod.fit(y, pars, x=x)

    print(output.fit_report())

    x = x_ev

    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
    axes[0].plot(x, y, 'b')
    axes[0].plot(x, init, 'k--', label='initial fit')
    axes[0].plot(x, output.best_fit, 'r-', label='best fit')
    axes[0].legend(loc='best')

    comps = output.eval_components(x=x)
    axes[1].plot(x, y, 'b')
    axes[1].plot(x, comps['g1_'], 'g--', label='Gaussian component 1')
    axes[1].plot(x, comps['g2_'], 'm--', label='Gaussian component 2')
    axes[1].legend(loc='best')

    plt.show()
    plt.savefig('gauss_fit_on_abs.png',format='png',dpi=300)