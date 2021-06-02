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
from scipy import signal,optimize,constants as const
import math
import matplotlib as mpl
from pylab import cm

plt.rcParams.update({'figure.max_open_warning': 0}) #Remove figure creation RuntimeWarning.

def parse_mcd(path):
    d={}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "end" in str.lower(name):
                field = re.search('(.*)(?=T).',name).group(0) #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                f = field + ', ' + str(num) #differentiate repeated scans at same field
                # print("Adding", t + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['temperature']=float(temperature) #add column for field
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                d[f] = df #send dataframe to dictionary
                pass
            if "test" not in str.lower(name): #remove any test files performed during data acqusition
                temperature = re.search('(.*)(?=K)',name).group(0).replace('-', '.') #search for beginning of file name "#T_" and set as name for df. Man this was difficult to figure out syntactically. 
                t = temperature + ', ' + str(num) #differentiate repeated scans at same field
                # print("Adding", t + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['temperature']=float(temperature) #add column for field
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                d[t] = df #send dataframe to dictionary


                # df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                # df['temperature']=float(75) #add column for field
                # df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                # df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                # d['0T' + str(num)] = df #send dataframe to dictionary
    return d

# uncomment below to test dict format
    # print(d['-4T_0']['field'])
    # sys.exit()

def plot_mcd(dic,op='avg',x_axis='Energy (eV)',title='[PH]',xdata='energy',ydata='mdeg'):
    plt.clf()
    fig,ax=plt.subplots(figsize=(4,2))
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(0,76,76),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='coolwarm_r',norm=norm) 
    fig.colorbar(sm,ticks=range(0,76,10),label='T (K)') #make color bar based on H (T) for plot
    for df in dic.values():
        #Dr. Seaborn or: How I Learned to Stop Worrying and Love sns.lineplot. Such efficiency. Much wow.
        sns.lineplot(data=df,x=xdata,y=ydata, linewidth=0.6,
                    hue='temperature',hue_norm=(0,76),
                    palette=sns.color_palette('coolwarm_r',as_cmap=True),
                    legend=None)
    if x_axis=='Energy (eV)':
        ax.set_xlabel(x_axis)
    if op=='raw':
        plt.title("Raw MCD " + title)
    if op=='avg':
        plt.title("Averaged MCD " + title)
    ax.set_ylabel('MCD (mdeg)')

    # ax.xaxis.set_tick_params(which='major', size=5, width=1, direction='in', top='on')
    # ax.xaxis.set_tick_params(which='minor', size=2, width=1, direction='in', top='on')
    # ax.yaxis.set_tick_params(which='major', size=5, width=1, direction='in', right='on')
    # ax.yaxis.set_tick_params(which='minor', size=2, width=1, direction='in', right='on')

    # ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    # ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
    # ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
    # ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10))

    ax.set_xlim(1.75,.55)
    ax.set_ylim(-20,5)
    
    # ax2=ax.twiny()
    # x_1, x_2 = ax.get_xlim()
    # ax2.set_xlim(eV_to_nm(x_1),eV_to_nm(x_2))

    # ax2.xaxis.set_tick_params(which='major', size=5, width=1, direction='in')
    # ax2.xaxis.set_tick_params(which='minor', size=2, width=1, direction='in')

    # ax2.xaxis.set_major_locator(mpl.ticker.FixedLocator(nm_to_eV(np.linspace(500, 1500, 3))))
    # ax2.xaxis.set_minor_locator(mpl.ticker.FixedLocator(nm_to_eV(np.linspace(300, 1700, 15))))

    # ax2.set_xticklabels(['500', '1000', '1500'])

    # ax2.set_xlabel(r'$\mathregular{\lambda}$ (nm)')
    # ax2.set_xlim(387.5, 2255)

    plt.style.use('seaborn-paper')
    plt.savefig(op + '_mcd_' + title,dpi=300,transparent=False,bbox_inches='tight')
    plt.show()

def calc_raw_avg_mcd(dic): #need to define this before finding the mcd difference
    df_avgs={}
    for name, df in dic.items():
        temperature = re.search('(.*)(?=,)',name).group(0) #set variable 'field' equal to int(field) from original dictionary
        if temperature not in df_avgs:
            df_avgs[temperature] = pd.DataFrame() #if field is not in new dictionary, create an empty dataframe
        if temperature in df_avgs:
            df_concat=pd.concat([df_avgs[temperature], df]).groupby(['wavelength'], as_index=False) #concatenate field entry with new df entry, so long as field is matching
            df_avgs[temperature] = df_concat #update dictionary entry with newly concatenated one
        df_avgs[temperature]=df_avgs[temperature].mean() #take the average of the concatenated df
    return df_avgs

def calc_diff_mcd(dic,op='sub'):
    df_diff={}
    for name, df in dic.items():
        df_diff[name] = pd.DataFrame()
        if name == '0': 
            # if op=='add':
            #     df_diff[name] = df + dic[list(name)[list(name).index(name)+3]]
            # elif op=='sub':
            #     df_diff[name] = df - dic[list(name)[list(name).index(name)+3]]
            del df_diff[name] #placeholder for now. In future would like to plot 0T field difference, but previous function looks like it deletes a version of 0 so no pair to subtract.
            pass
        elif '-' not in name: #loop only positive half of dictionary
            if op=='add':
                df_diff[name] = df + dic['-' + name] #add positive and negative dictionary entries
                df_diff[name]['energy'] = df['energy'] #fix energy back to original values
                df_diff[name]['field'] = df['field'] #fix field back to original values
            elif op=='sub':
                df_diff[name] = df - dic['-' + name] #subtract positive and negative dictionary entries
                df_diff[name]['energy'] = df['energy'] #fix field back to original values
                df_diff[name]['field'] = df['field']
        else:
            del df_diff[name]
            continue
    return df_diff

def plot_diff_mcd(dic,op='avg',x_axis='Energy (eV)'):
    fig,ax=plt.subplots(figsize=(4,2))
    # norm=plt.Normalize(-10,10) #optional to remove discrete H bar divisions
    norm=colors.BoundaryNorm(np.linspace(0,10,6),ncolors=256)
    sm=plt.cm.ScalarMappable(cmap='Greys',norm=norm) 
    fig.colorbar(sm,ticks=range(0,11,2),label='H (T)') #make color bar based on H (T) for plot
    for df in dic.values():
        #Dr. Seaborn or: How I Learned to Stop Worrying and Love sns.lineplot. Such efficiency. Much wow.
        sns.lineplot(data=df,x='energy',y='mdeg', linewidth=0.6,
                    hue='field',hue_norm=(0,10),
                    palette=sns.color_palette('Greys',as_cmap=True),
                    legend=None)
    if x_axis=='Energy (eV)':
        plt.xlabel(x_axis)
    if op=='raw':
        plt.title("Raw MCD")
    if op=='avg':
        plt.title("Difference MCD")
    plt.ylabel('MCD (mdeg)')
    plt.xlim(3.2,.55)
    baseline = lines.Line2D(range(6),np.zeros(1),c='black',ls='--',lw=0.6) #draw baseline at 0T
    ax.add_line(baseline) #add baseline to plot
    plt.style.use('seaborn-paper')
    plt.savefig('diff_mcd',dpi=100,transparent=True,bbox_inches='tight')
    plt.show()

def nm_to_eV(nm):
    eV = 1240 / nm
    return ["%.3f" % z for z in eV]

def eV_to_nm(eV):
    nm = 1240/eV 
    return "%.0f" % nm



'''-------------------------------FUNCTIONAL CODE BELOW-------------------------------'''


'''parse all data files'''
#Change these pathways if using from GitHub.
raw_mcd_dic = parse_mcd("/mnt/c/Users/roflc/Desktop/MCD DATA/5-1 CFS/VT-MCD 05-21-21 NIR 5-1/") #raw mcd data in dictionary
print(raw_mcd_dic)
'''fit raw and avg mcd straight from datafile - no workup'''
plot_mcd(raw_mcd_dic,'raw',title='sample') #plot raw experimental mcd data
df_avgs = calc_raw_avg_mcd(raw_mcd_dic) #
plot_mcd(df_avgs,'avg',title='sample')

'''mcd difference (no blank)'''
for name, df in df_avgs.items():
    df_avgs[name]['avg-0T'] = df_avgs[name]['mdeg'] - df_avgs['75K-0T']['mdeg']
plot_mcd(df_avgs,'avg',title='Diff_no_blank_0T_subbed',ydata='avg-0T')

sys.exit()

fit_diff=plot_CP_diff(df_abs['energy'],df_abs['absorbance'])
average_ev, std_dev_ev, average_m, std_dev_m, zdf = calc_effective_mass_and_plot(fit_diff,df_avgs)
print(zdf)
zdf.to_csv('zeeman_data.csv')

# m, b = np.polyfit(B_list,ev_list,1)
plt.clf() #Clear all previous plots
fig=plt.figure(figsize=(4,2))
ax1=fig.add_subplot(111) #need this to add separate series to same graph
ax1.scatter([x for x in zdf['B'] if x > 0], list(zdf.loc[zdf['B'] > 0,'E_Z']), label=r"$B(+)$",color="b")
ax1.scatter(np.absolute([x for x in zdf['B'] if x < 0]), list(zdf.loc[zdf['B'] < 0,'E_Z']), label=r"$B(-)$",color="r")
plt.legend(loc=0)
# plt.plot(B_list, m*B_list + b)
plt.xlabel('B (T)')
plt.xticks(np.arange(0,11,2))
ax1.set_xlim(1,11)
ax1.set_ylim(0,0.24)
plt.ylabel(r'$E_Z$ (meV)')
plt.savefig('mev_test_plot.png',dpi=200,bbox_inches='tight')
plt.show()

# '''write HTML file report'''
# # writeHTMLfile('mcd.html','11-11-2020')
writeHTMLfile_difference('mcd_difference.html','03-30-2021, Both Max Signals')

xldf = pd.DataFrame()
for field, d in diff_df.items():
    df = pd.DataFrame(d)
    df.dropna(axis=1,how='all',inplace=True)
    df.dropna(how='any',inplace=True)
    df.drop(['chpx','chpy','field','pemx','pemy','energy'],axis=1,inplace=True)
    rename_list = {'deltaA':'{}_deltaA'.format(field),'mdeg':'{}_mdeg'.format(field),'zero_subtracted':'{}_zero_subtracted'.format(field)}
    df.rename(columns=rename_list,inplace=True)
    # print(df)
    try:
        xldf = xldf.merge(df, how='inner', on='wavelength')
    except KeyError:
        xldf = df

xldf = xldf.reindex(sorted(list(xldf), key=lambda x: x.split('_')[-1]), axis=1)
xldf.insert(0,'energy', [1240/x for x in xldf['wavelength']])
xldf.set_index('wavelength',inplace=True)
xldf.to_csv('04-27-2021'+'_worked_up_diff_mcd.csv')

print("...\nDone!")
