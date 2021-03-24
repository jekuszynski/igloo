import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

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

def plot_abs(dic,xdata='nm',title='[PH]'):
    fig,ax=plt.subplots(figsize=(8,4))
    for name, df in dic.items():
        sns.lineplot(data=df,x=xdata,y=' A', linewidth=1, label=name, legend='brief')
    if xdata is 'nm':
        plt.xlabel('Wavelength (nm)')
    if xdata is 'eV':
        plt.xlabel('Energy (eV)')
    plt.title(title +' Absorbance')
    plt.ylabel('Absorbance (a.u.)')
    # plt.xlim(1.55,.75)
    # plt.ylim()
    plt.style.use('seaborn-paper')
    plt.savefig(title + '_abs',dpi=200,bbox_inches='tight')
    plt.show()

directory = '/mnt/c/users/roflc/Downloads/3-23-21_scan/'
# sample_list_path = directory + 'Sample Table.csv'

sample_files, ref_files = parse_folder(directory)
print(sample_files)

plot_abs(sample_files,title='CFS Tests')

# xlabel='Wavelength (nm)'
# ylabel='Abs'

# x=df[xlabel]
# y=df[ylabel]

# plt.plot(x,y,color='black')
# plt.xlabel(xlabel)
# plt.ylabel(ylabel)
# plt.xlim(275,1100)
# plt.ylim(0,1)
# plt.show()

# save_directory='/home/jkusz/github/igloo/abs/'
# plt.savefig(save_directory+name+'.png',dpi=300)