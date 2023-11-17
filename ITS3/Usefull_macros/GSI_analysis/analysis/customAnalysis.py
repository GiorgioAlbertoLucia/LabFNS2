#!/usr/bin/env python3

'''
    Script to perform all steps of a analysis. Steps can be performed separately and everything is 
    controlled from a yml configuration file.

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''

import numpy as np
import yaml
import uproot
from ROOT import TTree, TFile, gDirectory, TCanvas

import sys
sys.path.append('../..')
from utils.dataPreprocessor import NpyFileManager
from utils.spectrumAnalyzer import SpectrumAnalyzer
from utils.energyCal import EnergyCalibration
from utils.peakMonitoring import PeakMonitor

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def drawDistribution(inputTree:TTree, branchName:str, histSpec, canvasOutPath:str):
    '''
    Create a histogram for the charge distribution from given TTree and saves it on a canvas

    Parameters
    ----------
        - inputTree (TTree): TTree to extract the histogram from
        - branchName (str): name of the branch to draw the histogram from
        - histSpec ([int, float, float]): [nbins, xmin, xmax] - specifics of the histogram to draw
        - canvasOutPath (str): file to draw the canvas to
    '''

    print('\n'+ color.BOLD + color.CYAN + 'Data Display...' + color.END)
    inputTree.Draw(f'{branchName}>>hist({histSpec[0]},{histSpec[1]},{histSpec[2]})')
    hist = gDirectory.Get('hist')
    cv = TCanvas()
    hist.Draw('hist')
    cv.SaveAs(canvasOutPath)
    del hist
    del cv

def dataPreprocessing(config):
    '''
    Full preprocessing of data. Create a TTree from a .npy file (with or without calibration) and save it in a .root file.
    Configurations are set from a yaml configuration file.

    Parameters
    ----------
        - config (yaml): loaded yaml configuration file
    '''

    print('\n'+ color.BOLD + color.CYAN + 'Data Preprocessing...' + color.END)
    reader = NpyFileManager(config['npyFilePath'], config['threshold'])
    reader.loadData()
    if config['do_display']:    reader.display()
    if config['do_calib']:      reader.calibrate(config['npzFilePath'], config['savePreproc'])
    reader.plantTree(config['treeName'])
    reader.fillTree(option=config['option'], v_reset=config['v_reset'])
    reader.saveTree(config['rootFile'])

def fitAnalysis(config):
    '''
    Perform analysis of the signal spectrum by identifying its main peaks and fitting those structures with
    custom functions (multiple structures can be fitted at once). More information on the README or the spectrumAnalyzer
    implementation file.

    Parameters
    ----------
        - config (yaml): loaded yaml configuration file
    '''
    print('\n'+ color.BOLD + color.CYAN + 'Fit Analysis...' + color.END)

    spectrumAnalyzer = SpectrumAnalyzer(config['configAnalysis'][0])
    for cfgFile in config['configAnalysis']:    spectrumAnalyzer.run_analysis(cfgFile)

def energyCalibration(config):
    '''
    Draws the energy vs signal amplitude curve for many (or a single) set of points. Everything is plotted on the same canvas.
    Single TGraphs are also saved. Configurations are set from a yaml configuration file.

    Parameters
    ----------
        - config (yaml): loaded yaml configuration file
    '''
    
    print('\n'+ color.BOLD + color.CYAN + 'Energy Calibration...' + color.END)
    energy_calibration = EnergyCalibration(config['configEnCalib'])
    for run in range(config['ECruns']):   energy_calibration.calibration(run)
    energy_calibration.buildMultigraph()

def peakPositionMonitoring(config):
    '''
    Monitors the position of a structure of the spectrum over time to spot systematic effects that would determine a shift.
    To do this, the input dataset is divided in subdatasets and those are inspected separetely. Configurations are set 
    from a yaml configuration file.

    Parameters
    ----------
        - config (yaml): loaded yaml configuration file
    '''

    print('\n'+ color.BOLD + color.CYAN + 'Peak position monitoring...' + color.END)
    peak_monitor = PeakMonitor(config['configPeakMon'])
    peak_monitor.monitorPeakOverTime()
    peak_monitor.plotPeakPositionOverTime()
    peak_monitor.plotPeaksOverTime()
    

if __name__ == '__main__':

    config_file = '/home/curved/apts_new/apts-dpts-ce65-daq-software/configs/config_customAnalysis.yml'
    with open(config_file) as f:   config = yaml.load(f, Loader=yaml.FullLoader)

    if config['do_preproc']:    dataPreprocessing(config)
    
    inFile = TFile(config['rootFile'])
    inTree = inFile.Get(config['treeName'])

    if config['do_display']:
        for branchName, histSpec, canvasOutPath  in zip(config['branchNames'], config['histSpecs'], config['canvasOutPaths']):    
            drawDistribution(inTree, branchName, histSpec, canvasOutPath)

    if config['do_analysis']:   fitAnalysis(config)
    if config['do_en_calib']:   energyCalibration(config)
    if config['do_peak_mon']:   peakPositionMonitoring(config)




