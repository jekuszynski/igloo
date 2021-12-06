from functionsForMCD import *
from icecream import ic
import sys

if __name__ == '__main__':
    workingPath = '/home/jkusz/github/igloo/mcd/fitting/'
    os.chdir(workingPath)
    
    material = 'Quartz-FringeField'

    '''parse data'''
    signalData = parseMCD('/mnt/c/Users/roflc/Dropbox/Research/FSU/Strouse/Projects/MCD/Fringe Field Testing/11192021_Quartz-400-850nm_outsidemagnet/','signal')
    blankData = parseMCD('/mnt/c/Users/roflc/Dropbox/Research/FSU/Strouse/Projects/MCD/Fringe Field Testing/11192021_Quartz_Blank-400-850nm_outsidemagnet/','blank')

    averagedSignalData = calcAverageMCD(signalData)
    averagedBlankData = calcAverageMCD(blankData)

    # print(averagedSignalData, averagedBlankData)

    differenceMCDSignal = calcDifference(averagedSignalData, averagedBlankData, spectra='MCD')
    absorptionData = calcDifference(averagedSignalData, averagedBlankData, spectra='ABS')

    MCDCSVConversion = zeroSubtract(convertToCSV(differenceMCDSignal, spectra='MCD'))
    ABSCSVConversion = zeroSubtract(convertToCSV(absorptionData, spectra='ABS'))

    MCDCSVConversion.to_csv(material + '_worked_up_mcd.csv')
    ABSCSVConversion.to_csv(material + '_worked_up_absorption.csv')

