import sys
import os
import pandas as pd
import numpy as np
from matplotlib import gridspec, pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from scipy import optimize
from scipy.stats import linregress
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy.polynomial.polynomial as poly


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
    plt.savefig('fitting.png',format='png',dpi=300)
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

    raw_coefficients=poly.polyfit(x,y,9)
    first_derivative_coefficients = np.polyder(raw_coefficients)
    second_derivative_coefficients = np.polyder(raw_coefficients,2)
    fit=poly.polyval(x,raw_coefficients)
    first_der_fit = poly.polyval(x,first_derivative_coefficients)
    second_der_fit = poly.polyval(x,second_derivative_coefficients)

    plt.clf()
    plt.plot(x, y, 'ro',label='raw_data')
    # x_range = np.linspace(x[0],x[-1],1000)
    plt.plot(x,fit,'-',label='fit')
    plt.plot(x,first_der_fit, 'bo',label='1st der')
    plt.plot(x,second_der_fit, '--', label='2nd der')
    plt.legend(loc='best')
    plt.show()
    plt.savefig('2nd_deriv_abs.png',dpi=300)