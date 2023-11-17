'''
    Class to analyse a signal or energy spectrum, fitting peaks with custom functions and performing 
    this operation with different possible selections in the dataset, even recursively.

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''



import numpy as np
import pandas as pd
import yaml

import uproot
from ROOT import TH1D, TF1, TFile, TCanvas, TLegend, TLatex, gStyle, TMath
from ROOT import kGreen, kAzure, kOrange

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


class SpectrumAnalyzer:
    '''
    Class to analyse a signal or energy spectrum, fitting peaks with custom functions and performing 
    this operation with different possible selections in the dataset, even recursively.

    Parameters
    ----------
        - configFilePath (str): path to a yml configuration file
        - config (yaml object): configuration object that stores all informations from the configuration file
        - dataset (pd.DataFrame): uploaded data for the analysis
        - outputFile (TFile): output root file where all produced output is written
    
    '''
    
    def __init__(self, configFilePath):

        self.config = self.loadConfig(configFilePath)
        tree = uproot.open(self.config['dataFile'])[self.config['treeName']]
        self.dataset = tree.arrays(library='pd')
        self.outputFile = TFile(self.config['outputFile'], 'recreate')
    
    def loadConfig(self, configFilePath: str):
        '''
        Load the configurations from a .yml configuration file

        Parameters
        ----------
            - configFilePath (str): path to the configuration file
        '''
        with open(configFilePath) as f:   return yaml.load(f, Loader=yaml.FullLoader)
    
    def applyCuts(self, condition: str):
        '''
        Apply cuts to the data with a query string.

        Parameters
        ----------
            - condition (str): condition based on which elements of the dataset will be selected

        Returns
        -------
            - pd.DataFrame: Data with applied cuts.
        '''
        
        return self.dataset.query(condition, inplace=False)
    
    def analyzeDataset(self, data: pd.Series, suffix: str = None):
        '''
        Analyse dataset and perform fits.

        Args:
            - data (pd.Series): data to analyse
            - suffix (str): suffix to name of the produced objects
            - configFilePath (str): path to configuration file 
        '''

        gStyle.SetOptFit(0)

        if suffix is None:  
            cfg = self.config
            suffix = ''
        else:               cfg = self.config[suffix]

        histSpecs = [int(data.max()-data.min()), data.min(), data.max()]
        if self.config['histSpec'] is not None:     histSpecs = [self.config['histSpec'][0], self.config['histSpec'][1], self.config['histSpec'][2]]
        spectrum = TH1D('spectrum'+suffix, self.config['histTitle'], histSpecs[0], histSpecs[1], histSpecs[2])
        for x in data:   spectrum.Fill(x)
        spectrum.Rebin(self.config['rebin'])
        spectrum.SetStats(0)

        fitFuncs = []
        for fitFunc, fitName, fitTitle, fitRange, initPars, limPars, fitCol, drawRange in zip(cfg['fitFunction'], cfg['fitName'], cfg['fitTitle'], cfg['fitRange'], cfg['initialParameters'], cfg['limitParameters'], cfg['fitColor'], cfg['drawRange']):
            fitFunc = TF1(fitName+suffix, fitFunc, fitRange[0], fitRange[1])
            for idx, par in enumerate(initPars):    fitFunc.SetParameter(idx, par)
            for idx, parLim in limPars.items():      fitFunc.SetParLimits(idx, parLim[0], parLim[1])
            fitFunc.SetTitle(fitTitle)
            fitFunc.SetLineColor(fitCol)
            spectrum.Fit(fitFunc, 'rm+')
            print('#chi^{2} / NDF ='+f' {fitFunc.GetChisquare():#.0f} / {fitFunc.GetNDF()}', flush=True)
            fitFunc.SetRange(drawRange[0], drawRange[1])
            fitFuncs.append(fitFunc)

        canvas = TCanvas(self.config['prefix']+suffix, '', 900, 900)
        funcs = []
        start_par = 0
        for sing_idx, singleFuncs in enumerate(cfg['singleFuncs']):
            for idx, (expr, pars) in enumerate(singleFuncs.items()):   
                foo = TF1(f'{cfg["fitName"][sing_idx]}_{idx}', expr, cfg['drawRange'][sing_idx][0], cfg['drawRange'][sing_idx][1])
                foo.SetTitle(cfg['singleNames'][sing_idx][idx])
                for par in range(pars): foo.SetParameter(start_par+par, fitFuncs[sing_idx].GetParameter(start_par+par))
                foo.SetLineColor(cfg['singleColors'][sing_idx][idx])
                funcs.append(foo)
                start_par += pars

        spectrum.SetFillColorAlpha(self.config['histFillAlpha'][0], self.config['histFillAlpha'][1])
        spectrum.Draw('hist')
        for fitFunc in fitFuncs:    fitFunc.Draw('same')
        for foo in funcs:           foo.Draw('same')
        
        titleText = None
        if cfg['titleText'] is not None:
            titleText = TLatex(cfg['titleText'][0], cfg['titleText'][1], cfg['titleText'][2])
            titleText.SetNDC()
            titleText.SetTextSize(gStyle.GetTextSize())
            titleText.SetTextFont(42)
            titleText.Draw()
    
        #text2 = TLatex(cfg['titleText'][0], cfg['titleText'][1]-0.06,f'Fit range: ({cfg["fitRange"][0]}, {cfg["fitRange"][1]})')
        #text2.SetNDC()
        #text2.SetTextSize(gStyle.GetTextSize()*0.7)
        #text2.SetTextFont(42)
        #text2.Draw()
        #text3 = TLatex(cfg['titleText'][0], cfg['titleText'][1]-0.1,'#chi^{2} / NDF ='+f' {fitFunc.GetChisquare():#.0f} / {fitFunc.GetNDF()}')
        #text3.SetNDC()
        #text3.SetTextSize(gStyle.GetTextSize()*0.7)
        #text3.SetTextFont(42)
        #text3.Draw()

        texts = []
        for text in cfg['texts']:   
            Ltext = TLatex(text[0], text[1], text[2])
            Ltext.SetNDC()
            Ltext.SetTextSize(gStyle.GetTextSize()*0.7)
            Ltext.SetTextFont(42)
            texts.append(Ltext)
        for text in texts:  text.Draw()

        leg = TLegend(cfg['legendPos'][0], cfg['legendPos'][1], cfg['legendPos'][2], cfg['legendPos'][3])
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetTextSize(gStyle.GetTextSize()*0.7)
        leg.SetFillStyle(0)

        leg.AddEntry(spectrum, 'Sensor data', 'lf')
        for fitFunc in fitFuncs:    leg.AddEntry(fitFunc, fitFunc.GetTitle(), 'lf')
        for foo in funcs:   leg.AddEntry(foo, foo.GetTitle(), 'lf')
        leg.Draw('same')

        self.outputFile.cd()
        canvas.Write()

        
    def compare_results(self, results1, results2):
        '''
        Compare results obtained with different cuts.

        Args:
            results1: Results obtained from the first dataset.
            results2: Results obtained from the second dataset.

        Returns:
            None
        '''
        # Implement your comparison logic here
        pass
    
    def run_analysis(self, configFilePath: str = None):
        '''
        Perform data analysis for each dataset and save or compare results.

        Args:
            - configFilePath (str): path to configuration file 

        Returns:
            None
        '''

        if configFilePath is not None:
            self.config = self.loadConfig(configFilePath)
            tree = uproot.open(self.config['dataFile'])[self.config['treeName']]
            self.dataset = tree.arrays(library='pd')
        print('\n' + color.BOLD + f'Analyzing {self.config["prefix"]}...' + color.END, flush=True)

        self.analyzeDataset(self.dataset[self.config['spectrumBranch']])
        
        for suffix, condition in self.config['cuts'].items():
            cut_data = self.applyCuts(condition)
            self.analyzeDataset(cut_data[self.config['spectrumBranch']], suffix)
            
            # Save or compare results as needed
            #self.compare_results(cut_data, fitted_peaks)

        #self.outputFile.Close()

if __name__ == "__main__":
    analyzer = SpectrumAnalyzer('config.yaml')
    analyzer.run_analysis()
