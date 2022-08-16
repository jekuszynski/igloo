# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 17:10:46 2022

@author: Jason E. Kuszynski
"""
version = 'Ver: 1.0.0'

from tkinter import *
from tkinter import ttk, filedialog, messagebox
from pathlib import Path, PurePath
import re
import os
import numpy as np
import linecache as lc
import pandas as pd
import scipy.stats as st
import matplotlib

matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)


def findDir():
    dirname = Path(filedialog.askdirectory(title="Select Directory",initialdir=os.getcwd()))
    path.set(dirname)
    folderPath.delete(0, END) #deletes the current value
    folderPath.insert(0, path) #inserts new value assigned by 2nd parameter
    
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
            messagebox.showinfo(message='Unexpected files found in folder')
    elif polarization=='X':
        pol0Dataframe = pd.DataFrame.from_dict(yData_pol0, orient='columns', \
                                                   dtype=np.float64)
        pol0Dataframe.insert(0, 'Energy', xAxis)
        pol0Dataframe.set_index('Energy', inplace=True)
        pol1Dataframe = pd.DataFrame.from_dict(yData_pol1, orient='columns', \
                                               dtype=np.float64)
        pol1Dataframe.insert(0, 'Energy', xAxis)
        pol1Dataframe.set_index('Energy', inplace=True)
        fileDataframe = pol0Dataframe - pol1Dataframe  
    else:
        fileDataframe, pol0Dataframe, pol1Dataframe = 0, 0, 0
        messagebox.showinfo(message='Error with reading file')
    return fileDataframe, pol0Dataframe, pol1Dataframe
    
def dataStats(dataset, filecount, confidence):
    # print(dataset)
    DataAvg = sum(dataset)/filecount
    DataStdev = pd.DataFrame(np.dstack(dataset).std(axis=2),\
                                     index = DataAvg.index, \
                                     columns = DataAvg.columns)  
        
    zScore = st.norm.ppf(1-(1-confidence)/2)
    
    DataSterr = DataStdev / np.sqrt(filecount)
    
    DataMoE = DataSterr * zScore
    
    DataUpper = DataAvg + DataMoE
    DataLower = DataAvg - DataMoE
    return DataAvg, DataUpper, DataLower
    
def updateGraph():
    global avg0, upper0, lower0, avg1, upper1, lower1, avgX, upperX, lowerX
    
    pol0Data = []
    pol1Data = []
    polXData = []
    pol0count = 0
    pol1count = 0
    polXcount = 0
    ops = []
    
    activeOptions = datas.curselection()
    if not activeOptions:
        messagebox.showinfo(message='Please select at least one Y data option')
    for i in activeOptions:
        ops.append(datas.get(i))
        # print(ops)
    activePolarization = re.split(',', polarization.get().upper())
    # print(activePolarization)
    activeConfidence = float(confidence.get())
    # print(activeConfidence)
    for dirPath ,_ , fileList in os.walk(path.get()): 
        #open and read files in file list sequentially
        for fileName, pol in zip(fileList, activePolarization):
            # print(fileList, activePolarization)
            newfilePath = Path(dirPath) / Path(fileName)
            # print(newfilePath)
            data = parseData(newfilePath, pol)
            # print(data)
            if pol == '0':
                pol0Data.append(data[0])
                pol0count += 1
            elif pol == '1':
                pol1Data.append(data[0])
                pol1count += 1
            elif pol =='X':
                polXData.append(data)
                polXcount += 1
    alldata = (pol0Data, pol1Data, polXData)
    allcounts = (pol0count, pol1count, polXcount)
    # print(alldata, allcounts)
    for num, (data, count) in enumerate(zip(alldata, allcounts)):
        if count > 1:
            if num==0:
                avg0, upper0, lower0 = dataStats(data, count, activeConfidence)
            elif num==1:
                avg1, upper1, lower1 = dataStats(data, count, activeConfidence)
            elif num==2:
                avgX, upperX, lowerX = dataStats(data, count, activeConfidence)
    figure = Figure(figsize=(4, 4), dpi=100)
    axes = figure.add_subplot(1,1,1)
    canvas = FigureCanvasTkAgg(figure, mainframe)
    canvas.get_tk_widget().grid(row=4, column=2, rowspan=3)
    # print(graphdatatype.get())
    try: 
        for op in ops:
            if graphdatatype.get() == '0':
                axes.plot(avg0.index,avg0[op])
                axes.fill_between(upper0.index,upper0[op],lower0[op], alpha=0.5)
            if graphdatatype.get() == '1':
                axes.plot(avg1.index,avg1[op])
                axes.fill_between(upper1.index,upper1[op],lower1[op], alpha=0.5)
            if graphdatatype.get() == 'X':
                try:
                    axes.plot(avgX.index,avgX[op])
                    axes.fill_between(upperX.index,upperX[op],lowerX[op], alpha=0.5)
                except:
                    try:
                        avgX = polXData[0][0]
                        axes.plot(avg0.index,avg0[op])
                        upperX = 0
                        lowerX = 0
                    except NameError:
                        messagebox.showinfo(message='No XMCD found in folder')
    except NameError:
        messagebox.showinfo(message='No Folder Selected?')
    axes.set_xlabel('Energy (eV)')
    axes.set_ylabel('Voltage (V)')
            
def exportData():
    savepath = filedialog.asksaveasfilename(title="Save File As...",initialdir=os.getcwd(),
                         filetypes=(("Excel files", "*.xlsx"),("All files", "*.*")))
    try:
        avg0.to_csv(savepath + 'Avg0.csv')
        upper0.to_csv(savepath + 'Upper0.csv')
        lower0.to_csv(savepath + 'Lower0.csv')
    except NameError:
        messagebox.showinfo(message='No Pol 0 found in folder')
    try:
        avg1.to_csv(savepath + 'Avg1.csv')
        upper1.to_csv(savepath + 'Upper1.csv')
        lower1.to_csv(savepath + 'Lower1.csv')
    except NameError:
        messagebox.showinfo(message='No Pol 1 found in folder')
    try:
        avgX.to_csv(savepath + 'AvgX.csv')
        upperX.to_csv(savepath + 'UpperX.csv')
        lowerX.to_csv(savepath + 'LowerX.csv')
    except NameError:
        messagebox.showinfo(message='No XMCD found in folder')
    
root = Tk()
root.title("SXSTM File Processor")

mainframe = ttk.Frame(root, padding="3 6 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

path = StringVar()
buttonPath = ttk.Button(mainframe, text="Select Folder...", command=findDir).grid(column=0, row=0, sticky=W)
folderPath = ttk.Entry(mainframe, width=7, textvariable=path)
folderPath.grid(column=1, row=0, columnspan=1, sticky=(W,E))
folderPath.state(['disabled'])

ttk.Label(mainframe, text='Polarization (e.g. "0,1,X")').grid(column=0, row=1, sticky=W)
polarization = StringVar()
polarization.set("0,1,0,1")
polarization_entry = ttk.Entry(mainframe, width=7, textvariable=polarization)
polarization_entry.grid(column=1, row=1, sticky=(W,E))

ttk.Label(mainframe, text='Confidence (e.g. "0.95")').grid(column=0, row=3, sticky=W)
confidence = StringVar()
confidence.set("0.95")
confidence_entry = ttk.Entry(mainframe, width=7, textvariable=confidence)
confidence_entry.grid(column=1, row=3, sticky=(W,E))

ttk.Label(mainframe, text='Data Select').grid(column=0, row=4, sticky=W)
dataChoices = ['Sample_SR570 (V)', 'LIA Tip Ch1 (V)', 'LIA Sample Ch1 (V)']
dataChoicesVar = StringVar(value=dataChoices)
datas = Listbox(mainframe, listvariable=dataChoicesVar, selectmode='extended', height=5)
datas.select_set(0)
datas.grid(column=0, row=4, sticky=(W,E))

graphdatatype = StringVar()
pol0 = ttk.Radiobutton(mainframe, text='Pol 0', variable=graphdatatype, value='0') 
pol1 = ttk.Radiobutton(mainframe, text='Pol 1', variable=graphdatatype, value='1') 
polX = ttk.Radiobutton(mainframe, text='XMCD', variable=graphdatatype, value='X') 
graphdatatype.set("0")

ttk.Button(mainframe, text="Graph", command=updateGraph).grid(column=0, row=6, sticky=(W,E))

ttk.Button(mainframe, text="Export All Data in Folder", command=exportData).grid(column=1, row=6, sticky=E)

ttk.Label(mainframe, text=version).grid(column=2, row=0, sticky=(N,E))

for child in mainframe.winfo_children():
    child.grid_configure(padx=5,pady=5)

root.mainloop()
