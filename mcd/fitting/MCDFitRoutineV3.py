from functionsForMCD import *
import sys

if __name__ == '__main__':
    workingPath = '/home/jkusz/github/igloo/mcd/fitting/'
    os.chdir(workingPath)
    
    material = 'PEG_Au_Dropcast_Fringe_Field'

    '''parse data'''
    signalData = parseMCD('/mnt/c/Users/roflc/Dropbox/Research/FSU/Strouse/Projects/MCD/Fringe Field Testing/11152021_Au-PEG-Dropcast-400-850nm_outsidemagnet/','signal')
    blankData = parseMCD('/mnt/c/Users/roflc/Dropbox/Research/FSU/Strouse/Projects/MCD/Fringe Field Testing/11122021_Blank_400-850nm_filtersworking_PMT/','blank')

    averagedSignalData = calcAverageMCD(signalData)
    averagedBlankData = calcAverageMCD(blankData)

    # print(averagedSignalData, averagedBlankData)

    differenceMCDSignal = calcDifference(averagedSignalData, averagedBlankData, spectra='MCD')
    absorptionData = calcDifference(averagedSignalData, averagedBlankData, spectra='ABS')

    MCDCSVConversion = zeroSubtract(convertToCSV(differenceMCDSignal, spectra='MCD'))
    ABSCSVConversion = zeroSubtract(convertToCSV(absorptionData, spectra='ABS'))

    MCDCSVConversion.to_csv(material + '_worked_up_mcd.csv')
    ABSCSVConversion.to_csv(material + '_worked_up_absorption.csv')

    simulateMCDSpectra(ABSCSVConversion)
    average_ev, std_dev_ev, average_m, std_dev_m, zdf = calc_effective_mass_and_plot(ABSCSVConversion, MCDCSVConversion, 3, 1.5, 1)

    zdf.to_csv(material + '_zeeman_data.csv')

    peak1 = [2.5, 2]
    peak2 = [3.1, 2.6]
    peak3 = [1.9, 1.6]

    MCDPeaks = [peak1,peak2,peak3]
    ABSPeaks = []

    MCDAmplitudes = findAllPeaks(MCDCSVConversion, MCDPeaks, spectra='MCD')
    MCDAmplitudes.to_csv(material + '_MCD_Amp_Peaks.csv', index=False)

    # ABSAmplitudes = findAllPeaks(ABSCSVConversion, ABSPeaks, spectra='ABS')
    # ABSAmplitudes.to_csv(material + '_ABS_Amp_Peaks.csv', index=False)