import sys
import os
import re
import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly 
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as lines
import seaborn as sns
from scipy import signal,optimize
import math

d_abs={}
df_avgs={}
df_diff={}
df_blank_subtracted={}

def read_mcd_file(path):
    d={}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
            try: #try; except for if there is no field in file name
                field_name = re.search('(.*)(?=T_)..',name).group(0) #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                df['field']=int(re.search('(.*)(?=T)',name).group(0)) #add column for field
                f=field_name + str(num%3) #differentiate repeated scans at same field
                # print("Adding", f + "...") #uncomment to check files are being added/named correctly
            except:
                f=name
                # print("Adding", f + "...") #uncomment to check files are being added/named correctly
            df['energy']=1240/df['wavelength'] #calculate energy from wavelength
            df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
            d[f] = df #send dataframe to dictionary
    return d

def tick_function(x): #will be used to find ticks for top x-axis
    nm = 1240/x
    return ['%.1f' % z for z in nm]

def plot_mcd(units='energy',data_filter_list=None):
    if data_filter_list is None:
        data_filter_list=[]

    fig=plt.figure(figsize=(10,8))
    ax1=fig.add_subplot(111)
    # ax2=ax1.twiny()

    if units is 'ev':
        for name, df in data.items():
            if name not in data_filter_list:
                sns.lineplot(data=df,x='energy',y='mdeg',legend='auto',label=name)
        plt.legend(ncol=3,frameon=False)
        ax1.set_xlabel('Energy (eV)')
        ax1.set_xlim(3,1.1)

    elif units is'nm':
        for name, df in data.items():
            if name not in data_filter_list:
                sns.lineplot(data=df,x='wavelength',y='mdeg',legend='auto',label=name)
        plt.legend(loc=0,ncol=3,frameon=False)
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_xlim(400,1100)

    ax1.set_ylabel('MCD (mdeg)')
    ax1.set_ylim()
    # ax2.set_xlim(ax1.get_xlim())
    # ax2_tick_locations = np.array(1240/data['Test_3']['wavelength'])
    # ax2.set_xticks(ax2_tick_locations)
    # ax2.set_xticklabels(tick_function(ax2_tick_locations))
    # ax2.set_xlabel("Wavelength (nm)")
    sns.set_context('talk')


path = "/mnt/c/Users/roflc/Downloads/MCD 21-02-03/"
data = read_mcd_file(path)

data_list = data.keys()
filter_list1 = ['0T_Test_1', 'Test_2']
filter_list2 = ['Test_3', 'Test_4', 'Test_5', 'Test_6', 'Test_7', 'Test_8', 'Test_9', 'Test_10', 'Test_11']

plt1 = plot_mcd('nm',filter_list1)
plt.savefig('/mnt/c/Users/roflc/Downloads/210203_MCD_1.png',dpi=300,facecolor='w',transparent=True,bbox_inches='tight')
plt2 = plot_mcd('nm',filter_list2)
plt.savefig('/mnt/c/Users/roflc/Downloads/210203_MCD_2.png',dpi=300,facecolor='w',transparent=True,bbox_inches='tight')
plt.show()
