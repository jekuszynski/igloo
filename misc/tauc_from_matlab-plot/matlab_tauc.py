import sys
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rcParamsDefault, gridspec, pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import re
from scipy.signal import savgol_filter
from scipy.stats import linregress

### modified from https://stackoverflow.com/questions/52458870/compute-and-plot-tangent-lines-along-a-curve-produced-by-polynomial-regression-u
def ModelAndScatterPlot(xData, yData, firstDer, axes):

    # first the raw data as a scatter plot
    # axes.plot(xData, yData,  'D')

    # create data for the fitted equation plot
    xModel = np.linspace(min(xData), max(xData))
    yModel = np.polyval(fittedParameters, xModel)

    # now the model as a line plot
    # axes.plot(xModel, yModel)

    axes.set_xlabel('X Data') # X axis data label
    axes.set_ylabel('Y Data') # Y axis data label

    # polynomial derivative from numpy
    deriv = np.polyder(fittedParameters)

    # for plotting
    minX = min(xData)
    maxX = max(xData)

    # value of derivative (slope) at a specific X value, so
    # that a straight line tangent can be plotted at the point
    # you might place this code in a loop to animate
    pointVal = xData[np.where(firstDer == np.amax(firstDer))[0][0]] # find index of x-axis from first derivative max
    y_value_at_point = np.polyval(fittedParameters, pointVal)
    slope_at_point = np.polyval(deriv, pointVal)

    ylow = (minX - pointVal) * slope_at_point + y_value_at_point
    yhigh = (maxX - pointVal) * slope_at_point + y_value_at_point

    # now the tangent as a line plot
    x, y = [minX, maxX], [ylow, yhigh]
    # axes.plot(x, y)
    return x, y

    # plt.show()
    # plt.close('all') # clean up after using pyplot

if __name__ == "__main__":

    plt.rcParams["figure.dpi"]=150
    plt.rcParams["figure.facecolor"]="white"
    # plt.rcParams["figure.figsize"]=(6.4, 4.8)
    plt.rcParams["font.size"] = 18
    plt.rcParams.update({'font.family':'sans-serif'})
    plt.rcParams.update({'font.sans-serif':'Arial'})
    plt.rcParams['axes.linewidth'] = 1.5

    working_path = 'C:/Users/roflc/Dropbox/Research/FSU/Strouse/Code/igloo/misc/tauc_from_matlab-plot'
    os.chdir(working_path)

    max_ev = 1.63
    min_ev = 1.5

    # newcolumnNames = []
    rawData = pd.read_excel('2022_11_16_2xMAFA Uvvis for Jason.xlsx')
    rawData = rawData.T.drop_duplicates().T
    rawData.columns = rawData.columns.str.replace("\s\s+", "")
    rawData.rename(columns = {'Unnamed: 0' : 'wavelength'}, inplace = True)
    rawData.insert(loc=1, column='energy', value=1240/rawData['wavelength'])
    # print(rawData)

    graphWidth = 640
    graphHeight = 480
    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=150)
    axes = f.add_subplot(111)

    xData=rawData.loc[rawData['energy'].between(min_ev, max_ev, inclusive='both'),'energy'].values
    # print(xData)
    print("# of data points is: " + str(len(xData)))

    # fig,axes = plt.subplots()
    # fig,axes = plt.subplots(3, 1, constrained_layout=True, figsize = (10,6))
    # gs = gridspec.GridSpec(3, 1, figure=fig)

    bandGapsCold = []
    bandGapsHot = []

    for colnum, col in enumerate(rawData.columns):
        if not colnum in (0, 1):
            yData=rawData.loc[rawData['energy'].between(min_ev, max_ev, inclusive='both'),col].values
            yData = (yData*xData)**2 #tauc data
            # plt.plot(xData, yData)

            firstDer = savgol_filter(yData, 13, 11, deriv=1)
            fittedParameters = np.polyfit(xData, yData, 7)
            tanx, tany = ModelAndScatterPlot(xData, yData, firstDer, axes)

            tanLine = linregress(tanx, tany)
            m, b = tanLine[0], tanLine[1]
            Eg = np.round(-b/m, 3)
            col = re.split('\s+', col)
            if col[-1] == 'h':
                if col[0] == 'RT':
                    bandGapsHot.append([300, Eg])
                else: 
                    bandGapsHot.append([int(col[0]), Eg])
            else:
                if col[0] == 'RT':
                    bandGapsCold.append([300, Eg])
                else: 
                    bandGapsCold.append([int(col[0]), Eg])
            # plt.plot(xData, firstDer)
        else: 
            pass
    # rawData.plot(x='wavelength', y=rawData.columns[1:], kind='line', ax=axes[0,0])

    # plt.show()

    bandGaps = pd.DataFrame(bandGapsHot, columns=['Temp', 'Eg'])
    bandGaps.to_csv('raw_tauc_warm_data.csv', index=False)

    plt.scatter(bandGaps['Temp'], bandGaps['Eg'], c=bandGaps['Temp'], cmap=cm.jet)
    plt.xlabel(r'Temperature, T (K)')
    plt.xlim([0,305])
    axes.xaxis.set_major_locator(MultipleLocator(50))
    axes.xaxis.set_minor_locator(AutoMinorLocator(5))

    plt.ylabel(r'Energy, $h\nu$ (eV)')
    axes.yaxis.set_minor_locator(AutoMinorLocator())

    plt.title("2022/11/16 2xMAFA")
    # plt.show()
    plt.savefig('WarmRun.png', bbox_inches='tight', edgecolor='w')
    # print(bandGaps)
    print("Done!")
