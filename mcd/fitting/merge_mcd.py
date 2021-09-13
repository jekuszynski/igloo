import sys
import os
import pandas as pd
import numpy as np
from matplotlib import gridspec, pyplot as plt
import matplotlib.colors as colors
from pandas.core.arrays.categorical import contains
import seaborn as sns
from scipy import optimize
from scipy.stats import linregress
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def getWavelength(eV):
    return 1240/eV

def getEnergy(nm):
    return 1240/nm

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

    # use ### nm as conversion of NIR -> VIS data.
    for field in vis_field_data:
        if 'deltaA' in field[0]:
            vis_ref_point = float(vis_field_data.loc[vis_field_data[('wavelength','')] == wavelength_match, field])        
            nir_ref_point = float(nir_field_data.loc[nir_field_data[('wavelength','')] == wavelength_match, field])
            nir_to_vis_conversion = vis_ref_point/nir_ref_point
            print(nir_to_vis_conversion)
            nir_field_data[field] = nir_field_data[field] * nir_to_vis_conversion
        # merged_spectra = vis_field_data.merge(nir_field_data,????)

    sys.exit()

    print(nir_converted_data)
    sys.exit()
    x = new_data['wavelength']
    y = new_data['deltaA',10]
    return merged_spectra

if __name__ == '__main__':

    working_path = '/home/jkusz/github/igloo/mcd/fitting/temp/'
    os.chdir(working_path)

    vis_data_path = '/home/jkusz/github/strouse-data/for_publication/3-1CFS/vis_pub/3-1CFS_worked_up_diff_mcd.csv'
    nir_data_path = '/home/jkusz/github/strouse-data/for_publication/3-1CFS/adj_nir_pub/3-1CFS_worked_up_diff_mcd.csv'

    '''Parse Data'''
    vis_mcd_data = parse_csv(vis_data_path)
    nir_mcd_data = parse_csv(nir_data_path)

    '''Merge & Scale Data'''
    merge_spectra(vis_mcd_data,nir_mcd_data,900)

    sys.exit()

    '''Plot Data to check'''

    