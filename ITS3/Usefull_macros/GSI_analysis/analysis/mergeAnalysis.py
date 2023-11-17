#!/usr/bin/env python3
'''
    Script to analyse data from different data acquisition by merging them. 
    The data to be merged are intended to be taken at different vbb values and should undergo an energy calibration
    using results from the calibration curves obtained earlier.

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''

import os
import yaml
from array import array

import uproot
from ROOT import TFile, TTree, gDirectory, TF1, TH1D, TCanvas, gStyle, TLatex, TLegend

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

TTReeDict = {'D': 'd',
             'I': 'i'}

class MergeAnalysis:
    '''
    Class to analyse data from different data acquisition by merging them. 
    The data to be merged are intended to be taken at different vbb values and should undergo an energy calibration
    using results from the calibration curves obtained earlier.
    
    Parameters
    ----------
        - _outputPath (str): location for the output file
        - _outputTreeName (str): name of the TTree storing output
        - _outputFile (TFile): output file
        - _outputTree (TTree): output TTree
    '''

    def __init__(self, outputPath, outputTreeName):

        self._outputPath = outputPath
        self._outputTreeName = outputTreeName
        self._outputFile = None
        self._outputTree = None

        self._copied_branches = {}
        self._new_variable = None

    def preprocessing(self, inputPaths, inputTrees, branchesToCopy, newVar, processFunction, inputsList):
        '''
            Function to copy data from multiple TTrees in a given one. 
            For all inputs the same branches will be copied and the same function to define the variable in a new
            branch will be executed with given input parameters.

            LIMITS: the only variables that can be copied are D and I (reference: https://root.cern.ch/doc/master/classTTree.html#a76927d17a94d091d4ef1ed0ba7ef3383)
        
            Parameters
            ----------
                - inputPaths ([str]): paths to the input files
                - inputTrees ([str]): names of the input trees
                - branchesToCopy ({str: str}): dictionary (name: name/format) of the names of the branches to copy in the new file
                - newVar ([str]): [name, type] new variable name and type
                - processFunction (func): function to execute to define the new variable
                - inputsList ([[...]]): list of lists of arguments for the processing function
        '''

        print('\n' + color.BOLD + 'Planting tree...' + color.END)
        self._outputFile = TFile(self._outputPath, 'recreate')
        print(f'File created at {color.DARKCYAN + self._outputPath + color.END}')
        self._outputTree = TTree(self._outputTreeName, 'Output Tree')

        for branchName, branchType in branchesToCopy.items():
            if TTReeDict[branchType] == 'd':    self._copied_branches[branchName] = array(TTReeDict[branchType], [0.])
            elif TTReeDict[branchType] == 'i':  self._copied_branches[branchName] = array(TTReeDict[branchType], [0])
            self._outputTree.Branch(branchName, self._copied_branches[branchName], f'{branchName}/{branchType}')

        self._new_variable = array(TTReeDict[newVar[1]], [0])
        self._outputTree.Branch(newVar[0], self._new_variable, f'{newVar[0]}/{newVar[1]}')

        for inputPath, inputTree, inputs in zip(inputPaths, inputTrees, inputsList):
            self.copyAndProcessData(inputPath, inputTree, branchesToCopy, processFunction, *inputs)
        
        self._outputFile.cd()
        self._outputTree.Write()
        self._outputFile.Close()
    

    def copyAndProcessData(self, inputPath, inputTree, branchesToCopy, processFunction, *args):
        '''
            Function to copy data from multiple TTrees in a given one, as well as defining new variables from
            old data and storing them in a new file.
        
            Parameters
            ----------
                - inputPath (str): path to the input file
                - inputTree (str): name of the input tree
                - branchesToCopy ({str: str}): dictionary (name: format) of the names of the branches to copy in the new file
                - processFunction (func): function to execute to define the new variable
                - *args: arguments for the processing function
        '''

        inputTree = uproot.open(inputPath)[inputTree]
        inputData = inputTree.arrays(library='pd')

        for entry in range(len(inputData)):

            for branchName, branchType in branchesToCopy.items():   
                self._copied_branches[branchName][0] = inputData[branchName][entry]
            
            #self._new_variable[0] = processFunction(inputData['seed_signal_mV'][entry], *args)
            self._new_variable[0] = processFunction(inputData['seed_signal_mV'][entry], *args)
            self._outputTree.Fill()

    def mergedAnalysis(self, configPath):
        '''
            Produce a histogram from a branch of a TTree and then fit it.
            Configurations are uploaded from a yaml comfiguration file.

            Parameters
            ----------
                - configPath (str): path to the configuration file
        '''
        
        print('\n' + color.BOLD + 'Running merge analysis...' + color.END)
        assert self._outputTree is not None, 'Fill your tree first.'
        with open(configPath) as f:   config = yaml.load(f, Loader=yaml.FullLoader)

        tree = uproot.open(config['dataFile'])[config['treeName']]
        data = tree.arrays(library='pd')
        for selection in config['selections']:  data.query(selection, inplace=True)
        
        gStyle.SetOptFit(0)
        spectrum = TH1D('spectrum', config['histTitle'], config["histSpec"][0], config["histSpec"][1], config["histSpec"][2])
        for x in data[config['branchName']]:    spectrum.Fill(x)
        spectrum.SetFillColorAlpha(config['histFillAlpha'][0], config['histFillAlpha'][1])
        spectrum.SetStats(0)

        fitFuncs = []
        
        for fit in config['fitFuncs']:
            fitCfg = config[fit]
            fitFunc = TF1(fit, fitCfg['fitFunc'], fitCfg['fitRange'][0], fitCfg['fitRange'][1])
            fitFunc.SetTitle(fitCfg['funcName'])

            for idx, par in enumerate(fitCfg['initialParameters']):    fitFunc.SetParameter(idx, par)
            for idx, parLim in fitCfg['limitParameters'].items():      fitFunc.SetParLimits(idx, parLim[0], parLim[1])
            fitFunc.SetLineColor(fitCfg['color'])
            spectrum.Fit(fitFunc, 'rm+')
            fitFuncs.append(fitFunc)
            print('#chi^{2} / NDF ='+f' {fitFunc.GetChisquare():#.0f} / {fitFunc.GetNDF()}')

        canvas = TCanvas('canvas', '', 900, 900)

        spectrum.Draw('hist')
        for fitFunc, funcName in zip(fitFuncs, config['fitFuncs']):
            fitCfg = config[funcName]
            fitFunc.SetRange(fitCfg['drawRange'][0], fitCfg['drawRange'][1])
            fitFunc.Draw('same')
        
        
        titleText = None
        if config['titleText'] is not None:
            titleText = TLatex(config['titleText'][0], config['titleText'][1], config['titleText'][2])
            titleText.SetNDC()
            titleText.SetTextSize(gStyle.GetTextSize())
            titleText.SetTextFont(42)
            titleText.Draw()

        texts = []
        for text in config['texts']:   
            Ltext = TLatex(text[0], text[1], text[2])
            Ltext.SetNDC()
            Ltext.SetTextSize(gStyle.GetTextSize()*0.7)
            Ltext.SetTextFont(42)
            texts.append(Ltext)
        for text in texts:  text.Draw()

        leg = TLegend(config['legendPos'][0], config['legendPos'][1], config['legendPos'][2], config['legendPos'][3])
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetTextSize(gStyle.GetTextSize()*0.7)
        leg.SetFillStyle(0)

        leg.AddEntry(spectrum, 'Sensor data', 'lf')
        for fitFunc, funcName in zip(fitFuncs, config['fitFuncs']): leg.AddEntry(fitFunc, config[funcName]['funcName'], 'lf')
        leg.Draw('same')

        assert os.path.exists(self._outputPath), 'The file has not beed created yet.'
        self._outputFile = TFile(self._outputPath, 'update')
        canvas.Write()
        print(f'Analysis plot saved at {color.DARKCYAN + self._outputPath + color.END}')
        self._outputFile.Close()

def energyCalibration(x, *args):
    '''
        Convert to energy

        Parameters
        ----------
            - x (float): value to convert
            - args[0] (float): offset
            - args[1] (float): slope
    
        NOTE: slope and offset are supposed to be obtained from a fit (variable) vs energy, not the inverse
    '''
    assert args[1] != 0, 'The slope is 0, conversion not possible'
    return (x - args[0])/args[1]


if __name__ == "__main__":

    outputPath = '/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/mergeAnalysis.root'
    analysis = MergeAnalysis(outputPath, 'dataFeEnergy')
    
    inputPaths = [#'/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/0.0/source/apts_DAQ-000901010542292A_20230810_182606.root',
                  '/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/1.2/source/apts_DAQ-000901010542292A_20230810_182704.root',
                  '/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/2.4/source/apts_DAQ-000901010542292A_20230810_182752.root',
                  '/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/3.6/source/apts_DAQ-000901010542292A_20230810_182842.root',
                  '/home/curved/apts_new/apts-dpts-ce65-daq-software/Data/dataAcqFe/SecondAcq_10082023/4.8/source/apts_DAQ-000901010542292A_20230810_182931.root'
                  ]
    inputTrees = ['dataFe',
                  'dataFe',
                  'dataFe',
                  'dataFe',
                  'dataFe']
    branchesToCopy = {'seed_signal_mV': 'D',
                      'signal_over_th_mV': 'D',
                      'seed': 'I',
                      'cluster_size': 'I',
                      'cluster_shape': 'I'}
    inputsList = [#[-0.003, 7.7], 
                  [-0.001, 13.4],
                  [-0.0008, 17.3],
                  [-0.001, 19.1], 
                  [-0.002, 20.1]]
    analysis.preprocessing(inputPaths, inputTrees, branchesToCopy, ['energy', 'D'], energyCalibration, inputsList)

    analysisCfgPath = '/home/curved/apts_new/apts-dpts-ce65-daq-software/configs/config_MergeAnalysis.yml'
    analysis.mergedAnalysis(analysisCfgPath)
