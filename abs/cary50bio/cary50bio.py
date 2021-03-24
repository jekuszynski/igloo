from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

directory = '/mnt/c/users/roflc/Desktop/Python Test Data/'
name='Nickel(II)Sulfate'
filename = name + '.csv'
filepath = directory + filename

df=pd.read_csv(filepath,engine='python',skip_blank_lines=True,header=1,usecols=[0,1],skipfooter=48)
# print(df)

xlabel='Wavelength (nm)'
ylabel='Abs'

x=df[xlabel]
y=df[ylabel]

plt.plot(x,y,color='black')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xlim(275,1100)
plt.ylim(0,1)
plt.show()

save_directory='/home/jkusz/github/igloo/abs/'
plt.savefig(save_directory+name+'.png',dpi=300)