import os
import re
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

def plot_all_mcd(dic_list):
    #prepare figure with subplot spacings
    fig,ax=plt.subplots(2,2,figsize=(7,4.5))

    #add color bar
    fig.subplots_adjust(right=0.8, wspace=0.05, hspace=0.3)
    # use as needed (top=0.9, bottom=0.1, left=0.1, )
    norm=plt.Normalize(-10,10)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    cbar_ax=fig.add_axes([0.87, 0.095, 0.015, 0.8]) #change for position [0(x),1(y)] and size [2(w),3(h)]
    cb=fig.colorbar(sm,ticks=range(-10,11,2),label='H (T)',cax=cbar_ax) #make color bar based on H (T) for plot
    # Below: right aligns ticks successfully, however since graph won't work - commenting out for now.
    # for tick in cb.ax.yaxis.get_majorticklabels():
    #             tick.set_horizontalalignment("right")
    #             tick.set_x(2.75) 

    #add absorption peak matching
    # ax[1,0].plot([1.19,1.19],[-3,3],color='black',linestyle='--',linewidth='1') 
    # ax[1,1].plot([1.12,1.12],[-3,3],color='black',linestyle='--',linewidth='1') 
    # ax[1,2].plot([1.19,1.19],[-3,3],color='black',linestyle='--',linewidth='1') 

    #make all subplots
    for num, dic in enumerate(dic_list):
        row = math.floor(num/2)
        col = num%2
        for df in dic.values():
            sns.lineplot(ax=ax[row,col],data=df,x='energy',y='avg-0T',
                    linewidth=0.6,hue='field',hue_norm=(-10,10),
                    palette=sns.color_palette('coolwarm_r',as_cmap=True),
                    legend=None)
        ax[row,col].plot([-10,10],[0,0],color='black',linestyle='-',linewidth='1') #add 0T baseline
        if col == 0:
            ax[row,col].set_ylabel(r'$\Delta$A/A$_{\mathrm{max}}$ (x $10^{-3}$)')
        if col != 0: 
            ax[row,col].set_ylabel('')
        if col == 1:
            ax[row,col].yaxis.tick_right()
            ax[row,col].yaxis.set_ticks_position('both')

            ### One day I'll figure out how to align right align right axis ticks.
            # for tick in ax[row,col].yaxis.get_major_ticks():
            #     tick.label2.set_horizontalalignment('right')

        if row == 0:
            ax[row,col].set_xlim(2.7,1.2)
            ax[row,col].set_ylim(-1.5,1.5)
            ax[row,col].set_xlabel('') #set xlabel for top row as empty.
            ax[row,col].xaxis.set_major_locator(MultipleLocator(0.4))
            ax[row,col].xaxis.set_minor_locator(MultipleLocator(0.1)) # auto set minor ticks
        if row == 1:
            ax[row,col].set_xlim(1.65,0.8) 
            ax[row,col].set_ylim(-1.0,1.0)
            ax[row,col].set_xlabel(r'Energy (eV)') #set xlabel for bottom row only.
            ax[row,col].xaxis.set_major_locator(MultipleLocator(0.2))
            ax[row,col].xaxis.set_minor_locator(MultipleLocator(0.1)) # auto set minor ticks
        ax[row,col].yaxis.set_major_locator(MultipleLocator(0.5))
        ax[row,col].yaxis.set_minor_locator(AutoMinorLocator()) # auto set minor ticks
        
        ax2 = ax[row,col].twiny() # creates a new axis with invisible y and independent x on opposite side of first x-axis
        if row == 0:
            ax2.set_xlabel(r'Wavelength (nm)') #set wavelength label for top row only.
        ax2.set_xscale('function',functions=(getWavelength,getEnergy)) # set twin scale (convert degree eV to nm)
        xmin, xmax = ax[row,col].get_xlim() # get left axis limits
        ax2.set_xlim((getWavelength(xmax),getWavelength(xmin))) # apply function and set transformed values to right axis limits
        ax2.xaxis.set_major_locator(MultipleLocator(300)) #set wavelength x-axis spacing
        ax2.xaxis.set_minor_locator(MultipleLocator(50)) #set minor ticks
        
        ax2.plot([],[]) # set an invisible artist to twin axes to prevent falling back to initial values on rescale events

        #Set tick parameters    
        ax[row,col].xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in') #axes.linewidth by default for matplotlib is 0.8, so that value is used here for aesthetic matching.
        ax[row,col].xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in')
        ax[row,col].yaxis.set_tick_params(which='major', size=5, width=0.8, direction='in', right='on')
        ax[row,col].yaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in', right='on')

        ax2.xaxis.set_tick_params(which='major', size=5, width=0.8, direction='in')
        ax2.xaxis.set_tick_params(which='minor', size=2, width=0.8, direction='in')

    ax[0,0].text(0.07,0.9,'A',horizontalalignment='center',verticalalignment='center',transform=ax[0,0].transAxes,fontweight='bold',fontsize=14)
    ax[1,0].text(0.07,0.9,'B',horizontalalignment='center',verticalalignment='center',transform=ax[1,0].transAxes,fontweight='bold',fontsize=14)
    ax[0,0].text(0.5,0.9,r'Cu$_7$FeS$_4$',horizontalalignment='center',verticalalignment='center',transform=ax[0,0].transAxes,fontstyle='normal',fontsize=12)
    ax[0,1].text(0.5,0.9,r'Cu$_3$FeS$_4$',horizontalalignment='center',verticalalignment='center',transform=ax[0,1].transAxes,fontstyle='normal',fontsize=12)

    # plt.tight_layout()
    plt.savefig('Figure2_no5-1.pdf',transparent=False,bbox_inches='tight')
    plt.savefig('Figure2_no5-1.png',transparent=False,bbox_inches='tight',dpi=300)
    plt.show()

'''-------------------------------FUNCTIONAL CODE BELOW-------------------------------'''


plt.clf() #Clear all previous plots

df_7_CFS_VIS = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/7-1 CFS/VIS/MCD 04-06-21 VIS Neg/")) 
# df_5_CFS_VIS = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/5-1 CFS/MCD 05-17-21 VIS 5-1/"))
df_3_CFS_VIS = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/3-1 CFS/VIS/MCD 04-08-21 VIS Both/"))
df_7_CFS_NIR = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/7-1 CFS/NIR/MCD 04-02-21/"))
# df_5_CFS_NIR = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/5-1 CFS/MCD 05-18-21 NIR 5-1/"))
df_3_CFS_NIR = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/3-1 CFS/NIR/MCD 04-13-21 NIR 3-1/"))
 
df_7_CFS_VIS_2 = calc_raw_avg_mcd(parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/7-1 CFS/VIS/MCD 03-30-21 VIS/")) #use for second data set for merging

mcd_list=[df_7_CFS_VIS,
            # df_5_CFS_VIS,
            df_3_CFS_VIS,
            df_7_CFS_NIR,
            # df_5_CFS_NIR,
            df_3_CFS_NIR]

for df_avgs in mcd_list:
    for name, df in df_avgs.items():
        df_avgs[name]['avg-0T'] = (df_avgs[name]['mdeg'] - df_avgs['0']['mdeg']) / 32982 * 1000

for name, df in df_7_CFS_VIS_2.items():
    df_7_CFS_VIS_2[name]['avg-0T'] = (df_7_CFS_VIS_2[name]['mdeg'] - df_7_CFS_VIS_2['0']['mdeg']) / 32982 * 1000
df_7_CFS_VIS.update(df_7_CFS_VIS_2) #change variable to whichever df needs to be updated above.

plot_all_mcd(mcd_list)

print("...\nDone!")
