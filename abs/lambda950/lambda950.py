import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def parse_folder(path):
    d={}
    d2={}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if ".sp" in str.lower(name):
                pass
            elif ".sample" in str.lower(name):
                df=pd.read_csv(path+name, sep=',') #create dataframe for each file
                df['eV']=1240/df['nm'] #calculate energy from wavelength
                df.set_index('nm')
                d[name] = df #send dataframe to dictionary
                pass
            elif ".correction" or "table" in str.lower(name):
                df=pd.read_csv(path+name, sep=',', index_col=0) #create dataframe for each file
                d2[name] = df
                pass
            else:
                print("Unexpected filetype found. Skipping.")
                pass
    return d, d2

def plot_abs(dic,xdata='nm',title='[PH]',save=''):
    fig,ax=plt.subplots(figsize=(8,4))
    for name, df in dic.items():
        sns.lineplot(data=df,x=xdata,y=' A', linewidth=1, label=name)
        wavelength_peak = df.loc[df[' A'] == df[' A'].loc[150:1000].max(),'nm'].iloc[0]
        print(str(wavelength_peak) + ' at ' + name.replace('.', ',').split(',')[1])
        if xdata == 'nm':
            plt.plot([wavelength_peak,wavelength_peak],[0,10],'k--',lw=0.8, label=str(wavelength_peak) + ' nm')
        elif xdata == 'eV':
            wavelength_peak = round(1240/wavelength_peak,4)
            plt.plot([wavelength_peak,wavelength_peak],[0,10],'k--',lw=0.8, label=str(wavelength_peak) + ' eV')
    if xdata is 'nm':
        plt.xlabel('Wavelength (nm)')
        plt.xlim(400,3000)
        ax.xaxis.set_major_locator(MultipleLocator(200))
    if xdata is 'eV':
        plt.xlabel('Energy (eV)')
        plt.xlim(4,0.5)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(0,1.5)
    plt.title(title +' Absorbance')
    plt.ylabel('Absorbance (a.u.)')

    plt.legend(loc='best')
    plt.style.use('seaborn-paper')
    plt.savefig(save + title + '.png',dpi=200,bbox_inches='tight')
    plt.show()

if __name__ == '__main__':

    working_path = '/home/jkusz/github/igloo/abs/lambda950/testing/'
    os.chdir(working_path)

    data_path = '/mnt/c/users/roflc/Dropbox/Research/FSU/Strouse/Projects/WO3-x/UVVISNIR/WO2,9 scan/'
    save_directory='/mnt/c/users/roflc/Dropbox/Research/FSU/Strouse/Projects/WO3-x/UVVISNIR/figures/'
    # save_directory = '/home/jkusz/github/igloo/abs/lambda950/testing/'

    # sample_list_path = directory + 'Sample Table.csv'

    sample_files, ref_files = parse_folder(data_path)
    # print(sample_files)

    plot_abs(sample_files,title='WO2,9',save=save_directory)
    plot_abs(sample_files,xdata='eV',title='WO2,9_eV',save=save_directory)

    # plt.plot(x,y,color='black')
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # plt.xlim(320,2000)
    # plt.ylim(0,1)
    # plt.show()

    