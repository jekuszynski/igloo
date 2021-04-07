'''
SOURCE CODE CREDITED TO UNDERGRADUATE JASON KING
'''

#setup 
import numpy as np
import pandas as pd 
# from tkinter import filedialog
import matplotlib.pyplot as plt
import os 
import glob
import csv
import itertools

#Dialog box 
initialdir='/home/jasonk0401/Desktop/UROP/' #assign a directory 
#filepaths= filedialog.askopenfilename(multiple=True, initialdir=initialdir,title='Select a File', filetypes=(("ASC Files",".ASC"),))
#^opens dialogbox to select multiple files that end in .ASC
#print(filepaths)

def list_files(filepath, filetype):
    paths=[]
    for root, dirs, files in os.walk(filepath):
        for name in files:
            if name.lower().endswith(filetype.lower()):
                paths.append(name)
    return(paths)

filenames= list_files(initialdir, '.asc')
print(filenames)

#remove .ASC from filenames for legend
#create new list
molecule_names=[] 
#for everyname in filenames split name and .asc 
#append name to new list 
for name in filenames:
    m=os.path.splitext(name)[0]
    molecule_names.append(m)
print(molecule_names)

#Compliling Files Into List of Lists
df=pd.DataFrame()
for element in filenames:
    df2=pd.read_csv(element, usecols=['Wavelength','Intensity'])
    df=pd.concat([df,df2], axis=1)
print(df)

#create individual names for each col
cols = []
count = 1
for column in df.columns:
    if column=='Wavelength':
        cols.append(f'Wavelength_{count}')
        count+=1
        continue
    cols.append(column)
df.columns = cols
cols = []
count = 1
for column in df.columns:
    if column=='Intensity':
        cols.append(f'Intensity_{count}')
        count+=1
        continue
    cols.append(column)
df.columns = cols
print(df)

#Plotting Raw Data
test='Intensity'
IntCols = [idx for idx in df if idx.lower().startswith(test.lower())]
#print(IntCols)
df.plot(x='Wavelength_1', y=IntCols)
plt.title('Raw Data')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Non-Normalized Intensity(a.u.)')
plt.legend(molecule_names)

#Skip first 300 rows to get rid of rayleigh line
df_no_rayleigh=df.drop(df.index[0:50])
#leaves rows as 250-575
#reset index becasuse it causes problems later when concating 
df_no_rayleigh.reset_index(drop=True, inplace=True)
print(df_no_rayleigh)

#Plotting Non-Normalized w/out Rayleigh line n Wavelength (nm)
test='Intensity'
IntCols = [idx for idx in df_no_rayleigh if idx.lower().startswith(test.lower())]
#print(IntCols)
df_no_rayleigh.plot(x='Wavelength_1', y=IntCols)
plt.title('Non-Normalized w/out Rayleigh line')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Non-Normalized Intensity(a.u.)')
plt.legend(molecule_names)

#Wavelength (nm) and Normalized Intensity 
#create list of columns that start with Intensity
test='Intensity'
test2='Wavelength'
IntCols = [idx for idx in df_no_rayleigh if idx.lower().startswith(test.lower())]
WavCols= [idx for idx in df_no_rayleigh if idx.lower().startswith(test2.lower())]
df_normalized=pd.DataFrame()
#for columns in df if they start with Intensity 
#apply calculation to columns that start with Intensity 
for col in df_no_rayleigh:
    if col.startswith(test):
        df_normalized[IntCols] = (df_no_rayleigh[IntCols].apply(lambda x: x/x.max())) #divide each number in a column by the max of that column
    if col.startswith(test2):
        df_normalized[WavCols] = df_no_rayleigh[WavCols]
print(df_normalized)

#Plotting Normalized w/out Rayleigh line in Wavelength (nm)
#probably dont need this everytime since IntCols is defined above
test='Intensity'
IntCols = [idx for idx in df_normalized if idx.lower().startswith(test.lower())]
#print(IntCols)
df_normalized.plot(x='Wavelength_1', y=IntCols)
plt.title('Normalized w/out Rayleigh line')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Normalized Intensity(a.u.)')
plt.legend(molecule_names)

#Raman Sh cm^-1 and Normalized Intensity 
#setup numbers for conversion 
test='Intensity'
test2='Wavelength'
IntCols = [idx for idx in df_no_rayleigh if idx.lower().startswith(test.lower())]
WavCols= [idx for idx in df_no_rayleigh if idx.lower().startswith(test2.lower())]
Excitation_Wavelength=532
a=10**7/Excitation_Wavelength
#for columns in df if they start with Intensity 
#apply calculation to columns that start with Intensity 
df_raman_shift=pd.DataFrame()
for col in df:
    if col.startswith(test):
        df_raman_shift[IntCols] = (df_normalized[IntCols])
    if col.startswith(test2):
        df_raman_shift[WavCols]= (df_normalized[WavCols].apply(lambda x:-10**7/x + a))
print(df_raman_shift)

#Plotting Normalized w/out Rayleigh line in Wavenumber (cm^-1)
test='Intensity'
IntCols = [idx for idx in df_raman_shift if idx.lower().startswith(test.lower())]
#print(IntCols)2
df_raman_shift.plot(x='Wavelength_1', y=IntCols)
plt.title('Final')
plt.xlabel('Raman Shift (cm^-1)')
plt.ylabel('Normalized Intensity(a.u.)')
plt.legend(molecule_names)

IntCols1 = [idx for idx in df if idx.lower().startswith(test.lower())]
y1=IntCols1
IntCols2 = [idx for idx in df_no_rayleigh if idx.lower().startswith(test.lower())]
y2=IntCols2
IntCols3 = [idx for idx in df_normalized if idx.lower().startswith(test.lower())]
y3=IntCols3
IntCols4 = [idx for idx in df_raman_shift if idx.lower().startswith(test.lower())]
y4=IntCols4
#print

#figure, axes = plt.subplots(1, 4)
#df.plot(x='Wavelength_1', y=y1, ax=axes[0])
#plt.legend(molecule_names)
#df_no_rayleigh.plot(x='Wavelength_1', y=y2, ax=axes[1])
#plt.legend(molecule_names)
#df_normalized.plot(x='Wavelength_1', y=y3, ax=axes[2])
#plt.legend(molecule_names)
#df_raman_shift.plot(x='Wavelength_1', y=y4, ax=axes[3])
#plt.legend(molecule_names)
#plt.show()


#setup 4 subplots 
fig, ((ax1, ax2),(ax3, ax4))= plt.subplots(nrows=2, ncols=2,figsize=(10,5))
#raw data plot
ax1=plt.subplot(2,2,1)
df.plot(x='Wavelength_1', y=y1, ax=ax1)
ax1.set_title("Raw Raman Data")
plt.legend(molecule_names)
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (counts)')
#witout rayleigh line
ax2=plt.subplot(2,2,2)
df_no_rayleigh.plot(x='Wavelength_1', y=y2, ax=ax2)
ax2.set_title("Raman Data w/out Rayleigh Line")
plt.legend(molecule_names)
ax2.set_xlabel('Wavelength (nm)')
ax2.set_ylabel('Intensity (counts)')
#normalized
ax3=plt.subplot(2,2,3)
df_normalized.plot(x='Wavelength_1', y=y3, ax=ax3)
ax3.set_title("Normalized Intensity w/ Wavelength")
plt.legend(molecule_names)
ax3.set_xlabel('Wavelength (nm)')
ax3.set_ylabel('Normalized Intensity (counts)')
#raman shift
ax4=plt.subplot(2,2,4)
df_raman_shift.plot(x='Wavelength_1', y=y4, ax=ax4)
ax4.set_title("Normalized Raman Data w/ Raman Shift")
plt.legend(molecule_names)
ax4.set_xlabel('Raman Shift (cm^-1)')
ax4.set_ylabel('Normalized Intensity (counts)')
#add y x lables and titles
#fix spacing of graphs
#Intensity (counts)
#Normalized Intensity (counts)
plt.tight_layout(pad=0.4, w_pad=3.0, h_pad=3.5)
plt.savefig("/home/jasonk0401/Desktop/UROP/RamanGraphs.png", dpi=300)
plt.show()

