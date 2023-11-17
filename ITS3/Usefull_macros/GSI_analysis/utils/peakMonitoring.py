'''
    Class to monitor the position of the peak of a distribution in time. 
    Used to study systematic shifts due to e.g. temperature changes.

    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''

import yaml

import uproot
from ROOT import TH1D, TF1, TFile, TGraphErrors, TCanvas, TLegend, gStyle

class PeakMonitor:
    def __init__(self, configFilePath: str):

        self.config = self.loadConfig(configFilePath)
        tree = uproot.open(self.config['dataFile'])[self.config['treeName']]
        self.dataset = tree.arrays(library='pd')
        self.outputFile = TFile(self.config['outputFile'], 'recreate')

        self.intervals = self.createIntervals(len(self.dataset), self.config['nIntervals'])
        self.peak_positions = [] # as (start, stop, value, err)

    def loadConfig(self, configFilePath: str):
        '''
        Load the configurations from a .yml configuration file

        Parameters
        ----------
            - configFilePath (str): path to the configuration file
        '''
        with open(configFilePath) as f:   return yaml.load(f, Loader=yaml.FullLoader)

    def fitPeak(self, idx: int, data):
        hist = TH1D(f"hist_{idx}",  f'; {self.config["histXaxis"]}; {self.config["histYaxis"]}', self.config['histSpec'][0], self.config['histSpec'][1], self.config['histSpec'][2])
        for value in data[self.config['treeBranch']]:      hist.Fill(value)
        
        fit_func = TF1("fit_func", "gaus", self.config['fitRange'][0], self.config['fitRange'][1])
        for i, par in self.config['initFitPars'].items():   fit_func.SetParameter(i, par)
        
        hist.Fit(fit_func, "rm+")
        print('#chi^{2} / NDF ='+f' {fit_func.GetChisquare():#.0f} / {fit_func.GetNDF()}')
        self.outputFile.cd()
        hist.Write()
        return fit_func.GetParameter(1), fit_func.GetParameter(2)  # Return the mean of the Gaussian fit
    
    def createIntervals(self, totalEvents: int, numIntervals: int):
        '''
            Intended to be used inside the class. Creates equally spaced intervals.

            Parameters
            ----------
                - totalEvents (int): total number of events to divide in qeually spaced intervals
                - numIntervals (int): number of desired intervals
            
            Returns
            -------
                - [(int, int)]: list of intervals (start, stop)
        '''

        baseIntervalLength = totalEvents // numIntervals
        remainder = totalEvents % numIntervals

        intervals = []
        start = 0

        for i in range(numIntervals):
            
            end = start + baseIntervalLength

            if remainder > 0:
                end += 1
                remainder -= 1

            intervals.append((start, end))
            start = end

        return intervals

    def monitorPeakOverTime(self):

        for idx, (start, end) in enumerate(self.intervals):
            data_subset = self.dataset[start:end]
            peak_position, peak_err = self.fitPeak(idx, data_subset)
            self.peak_positions.append((start, end, peak_position, peak_err))

    def plotPeakPositionOverTime(self):

        graph = TGraphErrors(len(self.peak_positions))
        for i, (start, end, position, error) in enumerate(self.peak_positions):
            graph.SetPoint(i, (start + end) / 2, position)
            graph.SetPointError(i, 0, error)
        
        canvas = TCanvas("canvas", "Peak Position over Time", 1200, 900)
        graph.SetTitle("Peak Position over Time; Event index; Peak position (mV)")
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(797)
        graph.Draw("AP")
        
        self.outputFile.cd()
        canvas.Write()

    def plotPeaksOverTime(self):
        '''
            Plot histograms for different moment during the data acquisition
        '''

        canvas = TCanvas("hist_canvas", "Peak Position over Time", 1200, 800)

        leg = TLegend(self.config['legendPos'][0], self.config['legendPos'][1], self.config['legendPos'][2], self.config['legendPos'][3])
        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetTextSize(gStyle.GetTextSize()*0.7)
        leg.SetFillStyle(0)
    
        for idx in range(self.config['nIntervals']):    
            obj = self.outputFile.Get(f'hist_{idx}')
            obj.SetLineColor(self.config['colors'][idx])
            obj.SetFillColorAlpha(self.config['colors'][idx], 0.3)
            obj.SetTitle(f'Peak position over time; {self.config["histXaxis"]}; {self.config["histYaxis"]}')
            obj.SetStats(0)
            canvas.cd()
            obj.Draw('hist same')
            leg.AddEntry(obj, f'Events ({self.intervals[idx][0]}, {self.intervals[idx][1]})', 'lf')

        leg.Draw('same')
        
        self.outputFile.cd()
        canvas.Write()


if __name__ == '__main__':
    # Example usage

    configInputPath = ''

    peak_monitor = PeakMonitor(configFilePath=configInputPath)
    peak_monitor.monitor_peak_over_time()
    peak_monitor.plot_peak_position_over_time()

