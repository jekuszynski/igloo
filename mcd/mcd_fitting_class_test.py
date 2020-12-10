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

class MCD():
    def __init__(self,folder):
        path="/mnt/c/Users/roflc/Desktop/" + folder + '/"'
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

            name=[]
            wavelength=[]
            energy=[]
            chopperx=[]
            choppery=[]
            pemx=[]
            pemy=[]
            da=[]
            mdeg=[]

            name.append(os.listdir(file))

            for line in lines:
                col=line.split()
                wavelength.append(int(col[0]))
                chopperx.append(float(col[1]))
                choppery.append(float(col[2]))
                pemx.append(float(col[3]))
                pemy.append(float(col[4]))
                da.append(float(col[5]))
            
            self.wavelength=np.array(wavelength)
            self.chopperx=np.array(chopperx)
            self.choppery=np.array(choppery)
            self.pemx=np.array(pemx)
            self.pemy=np.array(pemy)
            self.da=np.array(da)
            self.name=np.array(name)


    def plot(self):
        plt.plot(self.wavelength,self.da)
        plt.show()

if __name__ == "__main__":
    mcd=MCD("MCD 11-11-20")
    print(MCD.da)
    # mcd.plot()