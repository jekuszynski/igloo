import os
import pandas as pd
import numpy as np
from scipy import signal

def parse_abs(path,sg_smoothing=43):
    d_abs = {}
    for root, dirs, files in os.walk(path): #walk along all files in directory given a path
        for num, name in enumerate(files): #for each file in list of files...
            if "blank" in str.lower(name): #select for blank and ems file
                dic_name="Blank"
                # print("Adding", dic_name + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                d_abs[dic_name] = df #send dataframe to dictionary
            if "ems" or "sample" in str.lower(name):
                dic_name="Ems"
                # print("Adding", dic_name + "...") #uncomment to check files are being added/named correctly
                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df['energy']=1240/df['wavelength'] #calculate energy from wavelength
                d_abs[dic_name] = df #send dataframe to dictionary
    df_abs=pd.DataFrame(data=(d_abs['Ems']['wavelength'], d_abs['Ems']['energy'], d_abs['Ems']['chpx'], d_abs['Blank']['chpx'])).transpose() #make dataframe from ems/blank dictionary
    df_abs.columns=['wavelength','energy','ems','blank'] #setup columns
    df_abs=df_abs[(df_abs != 0).all(1)] #remove rows with 0's.
    df_abs['absorbance']=(2-np.log10(100 * df_abs['ems'] / df_abs['blank'])) #calculate absorbance from emission and blank data
    df_abs['smoothed_absorbance']=signal.savgol_filter(df_abs['absorbance'],sg_smoothing,2) #smooth absorbance plot using Savitzky-Golay
    df_abs = df_abs[df_abs.wavelength < 1700] #remove collection greater than 1700 nm (used for InGaAs errors mainly)
    return df_abs

def parse_lambda_950_abs(path):
    d_abs = {}
    df_abs=pd.read_csv(path, header=0, names=['wavelength','absorbance'])
    df_abs['energy']=1240/df_abs['wavelength'] #calculate energy from wavelength
    df_abs['smoothed_absorbance']=signal.savgol_filter(df_abs['absorbance'],59,3) #smooth absorbance plot using Savitzky-Golay
    return df_abs

# xldf = pd.DataFrame()
# for field, d in df_avgs.items():
#     df = pd.DataFrame(d)
#     df.dropna(axis=1,how='all',inplace=True)
#     df.dropna(how='any',inplace=True)
#     df.drop(['chpx','chpy','field','pemx','pemy','energy'],axis=1,inplace=True)
#     df['avg-0T-deltaA']=df['avg-0T'] / 32982 * 1000 #give back deltaA x10^-3
#     rename_list = {'deltaA':'{}_deltaA'.format(field),'mdeg':'{}_mdeg'.format(field),'avg-0T':'{}_avg-0T'.format(field),'avg-0T-deltaA':'{}_avg-0T-deltaA'.format(field)}
#     df.rename(columns=rename_list,inplace=True)
#     # print(df)
#     try:
#         xldf = xldf.merge(df, how='inner', on='wavelength')
#     except KeyError:
#         xldf = df

#Make this a function
# xldf = xldf.reindex(sorted(list(xldf), key=lambda x: x.split('_')[-1]), axis=1)
# xldf.insert(0,'energy', [1240/x for x in xldf['wavelength']])
# xldf.set_index('wavelength',inplace=True)
# xldf.to_csv('CFS'+'_worked_up_diff_mcd.csv')

df_5_CFS_NIR = parse_abs("/mnt/c/Users/roflc/Desktop/MCD DATA/5-1 CFS/ABS 05-17-21 5-1/NIR/Use/")
df_5_CFS_VIS = parse_abs("/mnt/c/Users/roflc/Desktop/MCD DATA/5-1 CFS/ABS 05-17-21 5-1/VIS/")
df_3_CFS_NIR = parse_abs("/mnt/c/Users/roflc/Desktop/MCD DATA/3-1 CFS/NIR/ABS 04-08-21 NIR/")
df_3_CFS_VIS = parse_abs("/mnt/c/Users/roflc/Desktop/MCD DATA/3-1 CFS/VIS/ABS 04-07-21 VIS/")
df_7_CFS_VIS = parse_abs("/mnt/c/Users/roflc/Desktop/MCD DATA/7-1 CFS/VIS/ABS 03-29-21/")

abs_list=[df_5_CFS_NIR,df_5_CFS_VIS,df_3_CFS_NIR,df_3_CFS_VIS,df_7_CFS_VIS]
name_list=['df_5_CFS_NIR','df_5_CFS_VIS','df_3_CFS_NIR','df_3_CFS_VIS','df_7_CFS_VIS']

for num, x in enumerate(abs_list):
    df_abs = pd.DataFrame(x)
    df_abs.set_index('wavelength', inplace=True)
    df_abs.dropna(axis=1,how='all',inplace=True)
    df_abs.dropna(how='any',inplace=True)
    df_abs.drop(['ems','blank'],axis=1,inplace=True)
    df_abs.to_csv(str(name_list[num])+'_abs.csv')
    print(df_abs)

print("...\nDone!")