import os
import re
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

plt.rcParams.update({'figure.max_open_warning': 0}) #Remove figure creation RuntimeWarning.

def getWavelength(eV):
    return 1240/eV

def getEnergy(nm):
    return 1240/nm

def parse_mcd(path):
    d={}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "test" not in str.lower(name): #remove any test files performed during data acqusition
                field_name = re.search('(.*)(?=T_)..',name).group(0) #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                f=field_name + str(num%3) #differentiate repeated scans at same field
                # print("Adding", f + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['field']=int(re.search('(.*)(?=T)',name).group(0)) #add column for field
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                d[f] = df #send dataframe to dictionary
    return d

def plot_mcd(dic,op='avg',x_axis='Energy (eV)',title='[PH]',xdata='energy',ydata='mdeg'):
    plt.clf()
    fig,ax=plt.subplots(figsize=(6,4))
    norm=plt.Normalize(-10,10) #use to remove discrete H bar divisions on scale.
    # norm=colors.BoundaryNorm(np.linspace(-10,10,11),ncolors=256) #use to set discrete H bar divisions on scale.
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    fig.colorbar(sm,ticks=range(-10,11,2),label='H (T)') #make color bar based on H (T) for plot
    for df in dic.values():
        #Dr. Seaborn or: How I Learned to Stop Worrying and Love sns.lineplot. Such efficiency. Much wow.
        sns.lineplot(data=df,x=xdata,y=ydata, linewidth=0.6,
                    hue='field',hue_norm=(-10,10),
                    palette=sns.color_palette('coolwarm_r',as_cmap=True),
                    legend=None)
    
    ax.plot([-10,10],[0,0],color='black',linestyle='-',linewidth='1') #add 0T baseline

    ax.set_xlabel(r'Energy (eV)')
    ax.set_ylabel(r'$\Delta$A/A$_{\mathrm{max}}$ (x $10^{-3}$)')
    ax.set_xlim(2.7,1.2)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks
    ax.set_ylim(-1.5,1.5)
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
    plt.savefig(op + '_mcd_' + title + '.pdf',transparent=False,bbox_inches='tight')
    plt.savefig(op + '_mcd_' + title + '.png',transparent=False,bbox_inches='tight',dpi=300)
    plt.show()

def calc_raw_avg_mcd(dic): #need to define this before finding the mcd difference
    df_avgs={}
    for name, df in dic.items():
        field = re.search('(.*)(?=T)',name).group(0) #set variable 'field' equal to int(field) from original dictionary
        if field not in df_avgs: 
            df_avgs[field] = pd.DataFrame() #if field is not in new dictionary, create an empty dataframe
        if field in df_avgs:
            df_concat=pd.concat([df_avgs[field], df]).groupby(['wavelength'], as_index=False) #concatenate field entry with new df entry, so long as field is matching
            df_avgs[field] = df_concat #update dictionary entry with newly concatenated one
        df_avgs[field]=df_avgs[field].mean() #take the average of the concatenated df
    return df_avgs



'''-------------------------------FUNCTIONAL CODE BELOW-------------------------------'''


if __name__ == '__main__':
    working_path = '/home/jkusz/github/igloo/mcd/fitting/temp/'
    os.chdir(working_path)

    plt.clf() #Clear all previous plots

    '''parse all data files'''
    raw_mcd_dic = parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/7-1 CFS/VIS/MCD 03-30-21 VIS/") #raw mcd data in dictionary

    '''fit raw and avg mcd straight from datafile - no workup'''
    df_avgs = calc_raw_avg_mcd(raw_mcd_dic)

    '''mcd difference (no blank)'''
    for name, df in df_avgs.items():
        df_avgs[name]['avg-0T'] = (df_avgs[name]['mdeg'] - df_avgs['0']['mdeg']) / 32982 * 1000

    # Uncomment below if need to merge separate data sets
    raw_mcd_dic2 = parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/7-1 CFS/VIS/MCD 04-06-21 VIS Neg/")
    df_avgs2 = calc_raw_avg_mcd(raw_mcd_dic2)
    for name, df in df_avgs2.items():
        df_avgs2[name]['avg-0T'] = (df_avgs2[name]['mdeg'] - df_avgs2['0']['mdeg']) / 32982 * 1000
    df_avgs.update(df_avgs2)

    plot_mcd(df_avgs,'avg',title='7-1 CFS VIS',ydata='avg-0T')

    xldf = pd.DataFrame()
    for field, d in df_avgs.items():
        df = pd.DataFrame(d)
        df.dropna(axis=1,how='all',inplace=True)
        df.dropna(how='any',inplace=True)
        df.drop(['chpx','chpy','field','pemx','pemy','energy'],axis=1,inplace=True)
        df['avg-0T-deltaA']=df['avg-0T'] # / 32982 * 1000 #give back deltaA x10^-3
        rename_list = {'deltaA':'{}_deltaA'.format(field),'mdeg':'{}_mdeg'.format(field),'avg-0T':'{}_avg-0T'.format(field),'avg-0T-deltaA':'{}_avg-0T-deltaA'.format(field)}
        df.rename(columns=rename_list,inplace=True)
        # print(df)
        try:
            xldf = xldf.merge(df, how='inner', on='wavelength')
        except KeyError:
            xldf = df

    #Make this a function
    xldf = xldf.reindex(sorted(list(xldf), key=lambda x: x.split('_')[-1]), axis=1)
    xldf.insert(0,'energy', [1240/x for x in xldf['wavelength']])
    xldf.set_index('wavelength',inplace=True)
    xldf.to_csv('CFS'+'_worked_up_diff_mcd.csv')

    print("...\nDone!")


