'''
    Class to perform energy calibration with settings loaded from a configuration file.
    Intended to be imported and used outside of this macro.

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''


import pandas as pd
import numpy as np
import yaml

from ROOT import TGraphErrors, TFile, TF1, TCanvas, TMultiGraph, TLatex, TLegend, gStyle

class EnergyCalibration:
    def __init__(self, configFilePath:str):
        
        self.config = self.loadConfig(configFilePath)
        self.outputFile = TFile(self.config['outputFile'], 'recreate')

        self.multigraph = TMultiGraph('mg', self.config['graphTitle'])
        self.fitFuncs = []
        self.graphs = []
        self.params = []
    
    def loadConfig(self, configFilePath: str):
        '''
        Load the configurations from a .yml configuration file

        Parameters
        ----------
            - configFilePath (str): path to the configuration file
        '''
        with open(configFilePath) as f:   return yaml.load(f, Loader=yaml.FullLoader)

    def calibration(self, idx: int):
        '''
        Energy calibration for one of the runs

        Parameters
        ----------
            - idx (int): run index
        '''

        fitCfg = self.config['fitCfg'][idx]

        data = pd.read_csv(self.config['inputFile'], comment='#')
        data.query(f'{self.config["runColumn"]} == {self.config["runs"][idx]}', inplace=True)
    
        graph = TGraphErrors(len(data[self.config['xCol']]), np.array(data[self.config['xCol']], dtype=float), np.array(data[self.config['yCol']], dtype=float), np.zeros(len(data[self.config['xCol']])), np.array(data[self.config['yErr']], dtype=float))
        graph.SetName(f'cal_graph_{idx}')
        graph.SetTitle(f'{self.config["runColumn"]} = {self.config["runs"][idx]} V, '+self.config['graphTitle'])
        graph.SetMarkerStyle(20)

        fitFunc = TF1(f'cal_curve_{idx}', self.config['fitFunc'])
        fitFunc.SetTitle(f'Fit, {self.config["runColumn"]} = {self.config["runs"][idx]} V')
        fitFunc.SetLineColor(fitCfg['lineColor'])

        graph.Fit(fitFunc)
        self.graphs.append(graph)
        self.fitFuncs.append(fitFunc)

        canvas = TCanvas(f'canvas_cal_{idx}', f'{self.config["runColumn"]} = {self.config["runs"][idx]} V, '+self.config['graphTitle'])
        graph.Draw('ap same')
        fitFunc.Draw('same')
        print('#chi^{2} / NDF ='+f' {fitFunc.GetChisquare():#.0f} / {fitFunc.GetNDF()}')

        titleText = None
        if fitCfg['titleText'] is not None:
            titleText = TLatex(fitCfg['titleText'][0], fitCfg['titleText'][1], fitCfg['titleText'][2])
            titleText.SetNDC()
            titleText.SetTextSize(gStyle.GetTextSize())
            titleText.SetTextFont(42)
            titleText.Draw()
    
        #text1 = TLatex(fitCfg['titleText'][0], fitCfg['titleText'][1]-0.1,'#chi^{2} / NDF ='+f' {fitFunc.GetChisquare():#.0f} / {fitFunc.GetNDF()}')
        #text1.SetNDC()
        #text1.SetTextSize(gStyle.GetTextSize()*0.7)
        #text1.SetTextFont(42)
        #text1.Draw()
        text2 = TLatex(fitCfg['titleText'][0], fitCfg['titleText'][1]-0.14,f'Fit curve: {self.config["fitFunc"]}')
        text2.SetNDC()
        text2.SetTextSize(gStyle.GetTextSize()*0.7)
        text2.SetTextFont(42)
        text2.Draw()
        textParams = []
        for idx in range(self.config['nParams']):
            text = TLatex(fitCfg['titleText'][0], fitCfg['titleText'][1]-0.18-0.04*idx,f'[{idx}] = {fitFunc.GetParameter(idx):#.2f} #pm {fitFunc.GetParError(idx):#.2f}')
            text.SetNDC()
            text.SetTextSize(gStyle.GetTextSize()*0.7)
            text.SetTextFont(42)
            textParams.append(text)
        for text in textParams: text.Draw()


        leg = TLegend(fitCfg['legendPos'][0], fitCfg['legendPos'][1], fitCfg['legendPos'][2], fitCfg['legendPos'][3])
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetTextSize(gStyle.GetTextSize()*0.7)
        leg.SetFillStyle(0)

        leg.AddEntry(graph, 'Peak position', 'lp')
        leg.AddEntry(fitFunc, 'Fit function', 'lf')
        leg.Draw('same')
    
        self.outputFile.cd()
        canvas.Write()
    
        self.params.append([fitFunc.GetParameter(i) for i in range(self.config['nParams'])])

        return canvas

    def buildMultigraph(self):
        '''
        Build and save a multigraph using the previous graphs from calibration
        '''
    
        canvas = TCanvas('canvas_cal', self.config['graphTitle'])
        canvas.DrawFrame(-0.4, 0, 7.2, 140)
        for graph in self.graphs:   self.multigraph.Add(graph)
        self.multigraph.Draw('ap')
        for fitFunc in self.fitFuncs:   fitFunc.Draw('same')

        titleText = TLatex(self.config['legendPos'][0], self.config['legendPos'][3]+0.1, 'Energy calibration:')
        titleText.SetNDC()
        titleText.SetTextSize(gStyle.GetTextSize()*0.8)
        titleText.SetTextFont(42)
        titleText.Draw()

        text = TLatex(self.config['legendPos'][0], self.config['legendPos'][3]+0.04, 'Amplitude = A + B * Energy')
        text.SetNDC()
        text.SetTextSize(gStyle.GetTextSize()*0.8)
        text.SetTextFont(42)
        text.Draw()

        leg = TLegend(self.config['legendPos'][0], self.config['legendPos'][1], self.config['legendPos'][2], self.config['legendPos'][3])
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetTextSize(gStyle.GetTextSize()*0.7)
        leg.SetFillStyle(0)
        for fitFunc in self.fitFuncs:   leg.AddEntry(fitFunc, fitFunc.GetTitle(), 'lf')
        leg.Draw('same')

        self.outputFile.cd()
        canvas.Write()

# Example usage
if __name__ == '__main__':
    
    calibration = EnergyCalibration(degree=2)
    calibration.calibrate('energy_calibration_data.csv')
