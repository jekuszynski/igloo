import os
import re
import pandas as pd
import sys

def parsePDOS(path,data_type,verbose=False):
    d={}
    for _, _, files in os.walk(path): #walk along all files in directory given a path
        for name in files: #for each file in list of files...
            if "pdos_atm" in str.lower(name): #remove any test files performed during data acqusition
                atom, orbital = re.findall(r'\((.*?)\)',name) #parse filename to get atom and orbital
                raw_data = pd.read_table(path+name, sep='\s\s')
                column_data = raw_data.iloc[:,0:1].copy()
                if verbose == True:
                    print("Adding", atom + ' ' + orbital + " orbital...")
                print(column_data)
                sys.exit()

            

                df=pd.read_table(path+name, sep='\t',names=['wavelength','pemx','pemy','chpx','chpy','deltaA']) #create dataframe for each file
                df.insert(loc=1, column='energy', value=1240/df['wavelength']) #calculate energy from wavelength
                df['mdeg']=df['deltaA']*32982 # calculate mdeg from deltaA
                d[field + " " + scan_number + " " + data_type] = df #send dataframe to dictionary
    return d

if __name__ == '__main__':
    workingPath = '/home/jkusz/github/igloo/quantumESPRESSO/workingDir/'
    os.chdir(workingPath)
    
    material = 'Cu5FeS4'

    parsePDOS(workingPath, material, verbose=True)
    
    sys.exit()
    pdos = open

    pdos.to_csv(material + '_pdos.csv')

    print("Finished")