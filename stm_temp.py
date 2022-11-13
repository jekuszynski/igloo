import re
import os
import numpy as np
import linecache as lc
import pandas as pd
import scipy.stats as st
import sys

def parseData(path, polarization): 
    with open(path, 'rt') as lines: 
        #for each line in file...
        for lineNum, line in enumerate(lines, 1): 
            '''Begin line capture and x-axis (energy) assignments'''
            if 'EnergyScan ori' in line:
                centerScanEnergy = float(re.split('-|\s', line)[-2])
                # print(centerScanEnergy)
            if 'Energy Scanwidth (keV) ori' in line:
                energyScanWidth = float(re.split('-|\s', line)[-2])
                # print(energyScanWidth)
                startScanEnergy = centerScanEnergy - energyScanWidth/2 
                endScanEnergy = centerScanEnergy + energyScanWidth/2 
                # print(startScanEnergy, endScanEnergy)
            if 'Energy Pixels' in line:
                energyPixels = float(re.split('-|\s', line)[-2])
                # print(energyPixels)
                energyScanStepSize = energyScanWidth/energyPixels
                xAxis = np.arange(startScanEnergy, endScanEnergy, \
                                      energyScanStepSize)
                xAxis = np.around(xAxis, decimals=5)
                # print(xAxis)
            '''Setup y-axis data assignments'''
            #pull all data lines after header
            if polarization=='1' or polarization=='0':
                if '1D array data pol{}'.format(polarization) in line:
                    dataLines = range(lineNum + 1, lineNum + 34) 
                    headersList = []
                    yDataList = []
                    for i, dataLine in enumerate(dataLines):
                        newLine = lc.getline(str(path), dataLine)
                        if i%3 == 0:
                            header = re.split('\n', newLine)[0]
                            headersList.append(header)
                        if i%3 == 1:
                            data = re.split('\t|\n', newLine)[:-1]
                            yDataList.append(data)
                    yData = dict(zip(headersList, yDataList))
                    # print(yData)
                    yDataframe = pd.DataFrame.from_dict(yData, orient='columns', \
                                                            dtype=np.float64)
                    yDataframe.insert(0, 'Energy', xAxis)
            if polarization=='X':
                # print("Test1")
                if '1D array data pol0' in line:
                    # print("test2")
                    dataLines = range(lineNum + 1, lineNum + 34) 
                    headersList = []
                    yDataList = []
                    for i, dataLine in enumerate(dataLines):
                        newLine = lc.getline(str(path), dataLine)
                        if i%3 == 0:
                            header = re.split('\n', newLine)[0]
                            headersList.append(header)
                        if i%3 == 1:
                            data = re.split('\t|\n', newLine)[:-1]
                            yDataList.append(data)
                    yData_pol0 = dict(zip(headersList, yDataList))
                    
                if '1D array data pol1' in line:
                    # print("test3")
                    dataLines = range(lineNum + 1, lineNum + 34) 
                    headersList = []
                    yDataList = []
                    for i, dataLine in enumerate(dataLines):
                        newLine = lc.getline(str(path), dataLine)
                        if i%3 == 0:
                            header = re.split('\n', newLine)[0]
                            headersList.append(header)
                        if i%3 == 1:
                            data = re.split('\t|\n', newLine)[:-1]
                            yDataList.append(data)
                    yData_pol1 = dict(zip(headersList, yDataList))

                # print(pol0Dataframe, pol1Dataframe)
    if polarization=='1' or polarization=='0':
        try:
            fileDataframe = yDataframe.set_index('Energy')
            pol0Dataframe = 0
            pol1Dataframe = 0
        except UnboundLocalError:
            # messagebox.showinfo(message='Unexpected files found in folder')
            print('Unexpected files found in folder')
    elif polarization=='X':
        pol0Dataframe = pd.DataFrame.from_dict(yData_pol0, orient='columns', \
                                                   dtype=np.float64)
        pol0Dataframe.insert(0, 'Energy', xAxis)
        pol0Dataframe.set_index('Energy', inplace=True)
        pol1Dataframe = pd.DataFrame.from_dict(yData_pol1, orient='columns', \
                                               dtype=np.float64)
        pol1Dataframe.insert(0, 'Energy', xAxis)
        pol1Dataframe.set_index('Energy', inplace=True)
        fileDataframe = pol1Dataframe - pol0Dataframe
    else:
        fileDataframe, pol0Dataframe, pol1Dataframe = 0, 0, 0
        # messagebox.showinfo(message='Error with reading file')
        print('Error with reading file')
    return fileDataframe, pol0Dataframe, pol1Dataframe

def dataStats(dataset, confidence):
    # print(dataset)
    mean = dataset.groupby(['Energy'], as_index=True).mean()
    count = dataset.groupby(['Energy'], as_index=True).count()
    std = dataset.groupby(['Energy'], as_index=True).std()
    stderr = std / np.sqrt(count)
    zScore = st.norm.ppf(1-(1-confidence)/2)
    marginOfError = stderr * zScore
    upperCI = mean + marginOfError
    lowerCI = mean - marginOfError
    upperCI = upperCI.add_suffix('_UpperCI')
    lowerCI = lowerCI.add_suffix('_LowerCI')
    return mean, upperCI, lowerCI


#--------------------------------------------------#

folder = 'FF XMCD Yscan 2'
material = 'Cu7FeS4-dropcast'
pol = "X"
confidence = 0.95

filePath = 'C:/Users/roflc/Dropbox/Research/FSU/Strouse/Projects/CuFeS/SXSTM/Strouse_092022/' + material + '/' + folder + '/'
savePath = 'C:/Users/roflc/Dropbox/Research/FSU/Strouse/Projects/CuFeS/SXSTM/Strouse_092022/'

pol0Data = []
pol1Data = []
polXData = []
pol0count = 0
pol1count = 0
polXcount = 0


if __name__ == '__main__':
    for _, _, fileList in os.walk(filePath): 
        #open and read files in file list sequentially
        for fileName in fileList:
            xmcdData, rcpData, lcpData = parseData(filePath + fileName, pol)
            # xmcdData = (xmcdData - xmcdData.min()) / (xmcdData.max() - xmcdData.min())
            # rcpData = (rcpData - rcpData.min()) / (rcpData.max() - rcpData.min())
            # lcpData = (lcpData - lcpData.min()) / (lcpData.max() - lcpData.min())
            # xmcdData = xmcdData / xmcdData.abs().max()
            # rcpData = rcpData / rcpData.abs().max()
            # lcpData = lcpData / lcpData.abs().max()
                #bkgd subtraction
            xmcdData = xmcdData - xmcdData.iloc[-2]
            rcpData = rcpData - rcpData.iloc[-2]
            lcpData = lcpData - lcpData.iloc[-2]
            # print(data)
            if pol == '0':
                pol0Data.append(rcpData)
                pol0count += 1
            elif pol == '1':
                pol1Data.append(lcpData)
                pol1count += 1
            elif pol =='X':
                polXData.append(xmcdData)
                pol0Data.append(rcpData)
                pol1Data.append(lcpData)
                polXcount += 1
    xmcddataSet = pd.concat(polXData, axis=0, ignore_index=False).reset_index()
    lcpdataSet = pd.concat(pol1Data, axis=0, ignore_index=False).reset_index()
    rcpdataSet = pd.concat(pol0Data, axis=0, ignore_index=False).reset_index()

    finalxmcdData = pd.concat(dataStats(xmcddataSet, confidence), axis=1).sort_index(axis=1)
    finallcpData = pd.concat(dataStats(lcpdataSet, confidence), axis=1).sort_index(axis=1)
    finalrcpData = pd.concat(dataStats(rcpdataSet, confidence), axis=1).sort_index(axis=1)

    try:
        finalxmcdData.to_csv(savePath + folder + '_' + material + '_xmcd.csv')
        finallcpData.to_csv(savePath + folder + '_' + material + '_lcp.csv')
        finalrcpData.to_csv(savePath + folder + '_' + material + '_rcp.csv')
        print("Done!")
    except PermissionError: 
        print("You've left a goddamn file open'!")
    