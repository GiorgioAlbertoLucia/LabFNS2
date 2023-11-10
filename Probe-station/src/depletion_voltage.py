#!python3

'''
    Script to measure the depletion voltage from fits of a 1/C^2 vs V plot
    To run from the LabFNS2 folder:
        python3 -m Probe-station.src.depletion_voltage --input <input_file> --output <output_file>
        ./Probe-station/src/depletion_voltage.py --input INPUT_FILE --output OUTPUT_FILE

'''

import os
import argparse
import yaml

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import use as mpl_use
mpl_use('Agg')                              # Force matplotlib to not use any Xwindows backend (solves issues with ROOT)

#import uncertainties as unc

from ROOT import TFile, TGraphErrors, TF1, TCanvas, TLegend, gStyle, TLatex, TGraph2DErrors
from ROOT import kGreen, kAzure, kOrange
from ROOT import TColor


class CustomROOTPalette:
    def __init__(self):
        nRGBs = 5
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red = [1.00, 0.00, 0.00, 0.00, 1.00]
        green = [1.00, 0.00, 1.00, 0.00, 0.00]
        blue = [1.00, 1.00, 0.00, 0.00, 0.00]
        TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, 100)
        gStyle.SetNumberContours(100)


        
    


def yaml_load(file_path):
    """
    Load a YAML file and convert scientific notation values to float.

    Args:
        file_path (str): The path to the YAML file.

    Returns:
        dict: The YAML file contents as a dictionary.
    """
    with open(file_path, 'r') as f:
        yaml_data = yaml.safe_load(f)

    def convert_scientific_notation(value):
        try:
            return float(value)
        except ValueError:
            return value

    def recursive_convert_scientific_notation(data):
        if isinstance(data, dict):
            return {k: recursive_convert_scientific_notation(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [recursive_convert_scientific_notation(v) for v in data]
        else:
            return convert_scientific_notation(data)

    return recursive_convert_scientific_notation(yaml_data)

class DepletionAnalysis:

    def __init__(self, df: pd.DataFrame, args: argparse.Namespace):
        
        self.df = df
        self.sensor = args.sensor
        self.config = yaml_load(args.config)['detector'][args.sensor]
        self.plot_config = yaml_load(args.config_plot)['detector'][args.sensor]

        self.inversion = -1.
        if self.sensor == 'Strip' or self.sensor == 'strip':   self.inversion = 1.

        self.outFile = TFile(args.output, 'recreate')

        output_path = os.path.splitext(args.output)[0] + '.pdf'
        self.iC2vV_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'1C2_{args.sensor}_vs_V.pdf')
        self.zoom_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'1C2_{args.sensor}_vs_V_zoom.pdf')
        self.der_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'derivative_{args.sensor}.pdf')
        self.dcon_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'doping_concentration_{args.sensor}.pdf')
        self.dpro_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'doping_profile_{args.sensor}.pdf')
        self.dprocol_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'doping_profile_col_{args.sensor}.pdf')
        self.depth_outputPath = os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), f'depletion_depth_{args.sensor}.pdf')

        self.graph = None

        # standard texts for the plots
        self.color = int(self.config['color'])
        gStyle.SetOptStat(0)
        gStyle.SetOptFit(0)

        self.text1 = 'Area = 1.0 #times 1.0 mm^{2}'
        self.text2 = f'{args.sensor} sensor'

    def close(self):
        '''
            Close the output file 
        '''

        self.outFile.Close()

    #####################
    # DATA PREPROCESSING

    def preprocess_data(self):
        ''''
            Preprocess the data to calculate the doping concentration and the depletion width
        '''
        print("Preprocessing data...")

        eSi = 11.7 * 8.854e-12                                          # F/m - dielectric constant of silicon
        A = 1e-6                                                        # m^2 - detector area
        q = 1.602e-19                                                   # C - electron charge
        Vbi = -0.6                                                      # V - built-in voltage  

        self.df['V_abs'] = np.abs(self.df['V'] + self.inversion*Vbi)    # V - absolute value of the bias voltage
        self.df['1_C2'] = 1.0 / (self.df['C'] * self.df['C'])           # pF^-2
        self.df['1_C2_err'] = 2.0 * self.df['C_err'] / (self.df['C'] * self.df['C'] * self.df['C'])

        self.df['1_C2_F'] = self.df['1_C2'] * 1e24                      # F^-2
        self.df['1_C2_err_F'] = self.df['1_C2_err'] * 1e24              # F^-2

        self.df['derivative'] = np.gradient(self.df['1_C2_F'], self.inversion*self.df['V']) # F^-2 V^-1
        
        # error on the derivative
        self.df['derivative2_C'] = np.gradient(self.df['derivative'], self.df['1_C2_F'])
        self.df['derivative2_V'] = np.gradient(self.df['derivative'], self.inversion*self.df['V'])
        self.df['derivative_err'] = np.sqrt((self.df['derivative2_C']*self.df['1_C2_err_F'])**2 + (self.df['derivative2_V']*self.df['V_err'])**2)

        self.evaluate_W_NB()

    #####################
    # DATA ANALYSIS

    def print_zoom(self, xmin, xmax, ymin, ymax):
        '''
            Print the zoomed region of the plot

            Parameters
            ----------
            graph (TGraphErrors):   Graph with the data
            xmin (float):           Minimum x value
            xmax (float):           Maximum x value
            ymin (float):           Minimum y value
            ymax (float):           Maximum y value
            output_path (str):      Output path of the plot

            Returns
            -------
            None
        '''

        if self.graph is None:
            print('ERROR: graph is None')
            return

        self.graph.GetXaxis().SetRangeUser(xmin, xmax)
        self.graph.GetYaxis().SetRangeUser(ymin, ymax)

        canvas = TCanvas('1C2_zoom', 'canvas_zoom', 800, 600)
        self.graph.Draw('AP')    
        canvas.SetLeftMargin(0.15)

        legend = TLegend(0.45, 0.15, 0.65, 0.3)
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.AddEntry(self.graph, 'LGAD Data', 'pe')
        legend.Draw()

        latex = TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.04)
        latex.SetNDC()
        #latex.DrawLatex(0.45, 0.35, '#splitline{Zoom of non-linearity in the}{  depletion of the gain layer}')

        canvas.SaveAs(self.zoom_outputPath)
        self.outFile.cd()
        canvas.Write()

    def inverseC2_vs_V(self):
        '''
            Create a plot with 1/C^2 vs bias voltage

            Returns
            -------
            None
        '''

        plt_cfg = self.plot_config['ic2vV']
        
        graph = TGraphErrors(len(self.df['V']), 
                             np.asarray(self.inversion*self.df['V'], dtype=float), np.asarray(self.df['1_C2'], dtype=float), 
                             np.asarray(self.df['V_err'], dtype=float), np.asarray(self.df['1_C2_err'], dtype=float))
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(kAzure+1)
        graph.SetTitle('1/C^{2} vs V '+f'{args.sensor}'+'; Reverse bias (V); 1/C^{2} (pF^{-2})')

        self.graph = graph.Clone()
        
        # Fit the data
        
        # LGAD
        fit1 = TF1()
        if self.sensor == 'LGAD':   fit1 = self.fit(graph, self.config['fits'][0]['fit_range'][0], self.config['fits'][0]['fit_range'][1], 
                                                    init_fit_pars=[self.config['fits'][0]['parameters'][0]['value'], 
                                                                   self.config['fits'][0]['parameters'][0]['value']], 
                                                    lim_fit_pars=[[self.config['fits'][0]['parameters'][0]['min'], 
                                                                   self.config['fits'][0]['parameters'][0]['max']],
                                                                   [self.config['fits'][0]['parameters'][1]['min'], 
                                                                   self.config['fits'][0]['parameters'][1]['max']]],
                                                    line_color=797)
        fit2 = self.fit(graph, self.config['fits'][1]['fit_range'][0], self.config['fits'][1]['fit_range'][1], 
                        init_fit_pars=[self.config['fits'][1]['parameters'][0]['value'], 
                        self.config['fits'][1]['parameters'][0]['value']], line_color=862)
        fit3 = self.fit(graph, self.config['fits'][2]['fit_range'][0], self.config['fits'][2]['fit_range'][1], 
                        init_fit_pars=[self.config['fits'][2]['parameters'][0]['value'], 
                        self.config['fits'][2]['parameters'][0]['value']], line_color=kGreen+2)   

        if self.sensor == 'LGAD':   fit1.SetRange(self.config['fits'][0]['draw_range'][0], self.config['fits'][0]['draw_range'][1])
        fit2.SetRange(self.config['fits'][1]['draw_range'][0], self.config['fits'][1]['draw_range'][1])
        fit3.SetRange(self.config['fits'][2]['draw_range'][0], self.config['fits'][2]['draw_range'][1]) 

        # Find the intersection point
        intersection1 = 0.
        if self.sensor == 'LGAD':   
            intersection1 = self.find_intersection(fit2, fit1)
            intersection1_err = np.sqrt((fit1.GetParError(0))**2 + (fit2.GetParError(0))**2 + (np.mean(self.df['1_C2_err'])/(fit2.GetParameter(1)-fit1.GetParameter(1)))**2)
        intersection2 = self.find_intersection(fit3, fit2)
        intersection2_err = np.sqrt((fit2.GetParError(0))**2 + (fit3.GetParError(0))**2 + (np.mean(self.df['1_C2_err'])/(fit3.GetParameter(1)-fit2.GetParameter(1)))**2)

        canvas = TCanvas('1C2', 'canvas', 800, 600)
        canvas.DrawFrame(self.config['canvas_limits']['xmin'], self.config['canvas_limits']['ymin'],
                         self.config['canvas_limits']['xmax'], self.config['canvas_limits']['ymax'], '1/C^{2} vs V '+f'{args.sensor}'+'; Reverse bias [V]; 1/C^{2} [pF^{-2}]')
        #canvas.SetGrid()

        graph.Draw('P')
        if self.sensor == 'LGAD':   fit1.Draw('same')
        fit2.Draw('same')
        fit3.Draw('same')

        legend = TLegend(plt_cfg['legend'][0], plt_cfg['legend'][1], plt_cfg['legend'][2], plt_cfg['legend'][3])
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.AddEntry(graph, 'Data', 'pe')
        if self.sensor == 'LGAD':   legend.AddEntry(fit1, 'Gain layer depleted', 'l')
        legend.AddEntry(fit2, 'Bulk depletion', 'l')
        legend.AddEntry(fit3, 'Sensor fully depleted', 'l')
        legend.Draw()

        latex = TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.04)
        latex.SetNDC()
        if self.sensor == 'LGAD':
            latex.DrawLatex(plt_cfg['latex'][0][0], plt_cfg['latex'][0][1], f'Sensor: ({intersection2:.2f}'+'#pm'+f'{intersection2_err:.2f}) V')
            latex.DrawLatex(plt_cfg['latex'][1][0], plt_cfg['latex'][1][1], f'Gain layer: ({intersection1:.2f}'+'#pm'+f'{intersection1_err:.2f}) V')
            latex.DrawLatex(plt_cfg['latex'][2][0], plt_cfg['latex'][2][1], 'Depletion voltage:')
        else:   latex.DrawLatex(plt_cfg['latex'][0][0], plt_cfg['latex'][0][1], f'Sensor depletion: ({intersection2:.2f}'+'#pm'+f'{intersection2_err:.2f}) V')
        latex.DrawLatex(plt_cfg['latex'][3][0], plt_cfg['latex'][3][1], self.text1)
        latex.DrawLatex(plt_cfg['latex'][4][0], plt_cfg['latex'][4][1], self.text2)

        self.outFile.cd()
        canvas.Write()
        
        canvas.SaveAs(self.iC2vV_outputPath)

    # DOPING CONCENTRATION

    def derivative_plot(self):
        '''
            Plot the derivative of 1/C^2 vs bias voltage
        '''

        plt_cfg = self.plot_config['derivative']

        graph = TGraphErrors(len(self.df['V']), 
                             np.asarray(self.inversion*self.df['V'], dtype=float), np.asarray(self.df['derivative']*1e-6, dtype=float), 
                             np.asarray(self.df['V_err'], dtype=float), np.asarray(self.df['derivative_err']*1e-6, dtype=float) )
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(self.color)
        graph.SetTitle('Derivative '+f'{args.sensor}'+'; Reverse bias (V); #frac{d(1/C^2)}{dV} (cm^{-3} V^{-1})')

        canvas = TCanvas('derivative', 'canvas', 800, 600)
        canvas.DrawFrame(self.config['canvas_limits']['xmin'], self.config['canvas_limits']['der_ymin'], 
                         self.config['canvas_limits']['xmax'], self.config['canvas_limits']['der_ymax'], 
                         'Derivative '+f'{args.sensor}'+'; Reverse bias (V); #frac{d(1/C^2)}{dV} (cm^{-3} V^{-1})')
        canvas.SetLogy()

        graph.Draw('P')


        legend = TLegend(plt_cfg['legend'][0], plt_cfg['legend'][1], plt_cfg['legend'][2], plt_cfg['legend'][3])
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.AddEntry(graph, 'Data', 'pe')
        legend.Draw()

        latex = TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.04)
        latex.SetNDC()
        latex.DrawLatex(plt_cfg['latex'][0][0], plt_cfg['latex'][0][1], self.text1)
        latex.DrawLatex(plt_cfg['latex'][1][0], plt_cfg['latex'][1][1], self.text2)

        self.outFile.cd()
        canvas.Write()
        
        canvas.SaveAs(self.der_outputPath)

    def doping_concentration(self):
        '''
            Plot the doping concentration vs bias voltage
        '''

        plt_cfg = self.plot_config['dop_conc']

        graph = TGraphErrors(len(self.df['V']), 
                             np.asarray(self.inversion*self.df['V'], dtype=float), np.asarray(self.df['NB']*1e-6, dtype=float), 
                             np.asarray(self.df['V_err'], dtype=float), np.asarray(self.df['NB_err']*1e-6, dtype=float) )
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(self.color)
        graph.SetTitle('Doping concentration '+f'{args.sensor}'+'; Reverse bias (V); N_{B} (cm^{-3})')

        canvas = TCanvas('doping_conc', 'canvas', 800, 600)
        canvas.DrawFrame(self.config['canvas_limits']['xmin'], self.config['canvas_limits']['conc_ymin'], 
                         self.config['canvas_limits']['xmax'], self.config['canvas_limits']['conc_ymax'], 
                         'Doping concentration '+f'{args.sensor}'+'; Reverse bias (V); N_{B} (cm^{-3})')
        canvas.SetLogy()

        graph.Draw('P')


        legend = TLegend(plt_cfg['legend'][0], plt_cfg['legend'][1], plt_cfg['legend'][2], plt_cfg['legend'][3])
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.AddEntry(graph, 'Data', 'pe')
        legend.Draw()

        latex = TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.04)
        latex.SetNDC()
        latex.DrawLatex(plt_cfg['latex'][0][0], plt_cfg['latex'][0][1], self.text1)
        latex.DrawLatex(plt_cfg['latex'][1][0], plt_cfg['latex'][1][1], self.text2)

        self.outFile.cd()
        canvas.Write()
        
        canvas.SaveAs(self.dcon_outputPath)

    def doping_profile(self):
        '''
            Plot the doping profile vs depth of the depleted region

        '''

        plt_cfg = self.plot_config['dop_prof']

        graph = TGraphErrors(len(self.df['W']), 
                             np.asarray(self.df['W']*1e6, dtype=float), np.asarray(self.df['NB']*1e-6, dtype=float), 
                             np.asarray(self.df['W_err']*1e6, dtype=float), np.asarray(self.df['NB_err']*1e-6, dtype=float) )
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(self.color)
        graph.SetTitle('Doping profile '+f'{args.sensor}'+'; Depth (#mum); N_{B} (cm^{-3})')
        graph.GetYaxis().SetRangeUser(-0.1, 0.24)

        canvas = TCanvas('doping_prof', 'canvas', 800, 600)
        canvas.DrawFrame(self.config['canvas_limits']['depth_min'], self.config['canvas_limits']['conc_ymin'], 
                         self.config['canvas_limits']['depth_max'], self.config['canvas_limits']['conc_ymax'], 
                         'Doping profile '+f'{args.sensor}'+'; Depth (#mum); N_{B} (cm^{-3})')
        canvas.SetLogy()
        #canvas.SetGrid()

        graph.Draw('P')

        legend = TLegend(plt_cfg['legend'][0], plt_cfg['legend'][1], plt_cfg['legend'][2], plt_cfg['legend'][3])
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.AddEntry(graph, 'Data', 'pe')
        legend.Draw()

        latex = TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.04)
        latex.SetNDC()
        latex.DrawLatex(plt_cfg['latex'][3][0], plt_cfg['latex'][3][1], self.text1)
        latex.DrawLatex(plt_cfg['latex'][4][0], plt_cfg['latex'][4][1], self.text2)

        self.outFile.cd()
        canvas.Write()
        
        canvas.SaveAs(self.dpro_outputPath)

    def depletion_depth(self):
        '''
            Plot the depletion depth vs bias voltage
        '''

        plt_cfg = self.plot_config['depl_depth']

        graph = TGraphErrors(len(self.df['V']), 
                             np.asarray(self.inversion*self.df['V'], dtype=float), np.asarray(self.df['W']*1e6, dtype=float), 
                             np.asarray(self.df['V_err'], dtype=float), np.asarray(self.df['W_err']*1e6, dtype=float) )
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(self.color)
        graph.SetTitle('Depletion depth '+f'{args.sensor}'+'; Reverse bias (V); W (#mum)')

        canvas = TCanvas('depletion_depth', 'canvas', 800, 600)
        canvas.DrawFrame(self.config['canvas_limits']['xmin'], self.config['canvas_limits']['depth_min'], 
                         self.config['canvas_limits']['xmax'], self.config['canvas_limits']['depth_max'], 
                         'Depletion depth '+f'{args.sensor}'+'; Reverse bias (V); W (#mum)')

        graph.Draw('P')


        legend = TLegend(plt_cfg['legend'][0], plt_cfg['legend'][1], plt_cfg['legend'][2], plt_cfg['legend'][3])
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.AddEntry(graph, 'Data', 'pe')
        legend.Draw()

        latex = TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.04)
        latex.SetNDC()
        latex.DrawLatex(plt_cfg['latex'][3][0], plt_cfg['latex'][3][1], self.text1)
        latex.DrawLatex(plt_cfg['latex'][4][0], plt_cfg['latex'][4][1], self.text2)

        self.outFile.cd()
        canvas.Write()
        
        canvas.SaveAs(self.depth_outputPath)

    def doping_profile_color(self):
        '''
            Plot the doping profile vs depth of the depleted region

        '''

        # Define the color map
        cmap = plt.colormaps['cool']
        plt.scatter(self.df['W']*1e6, self.df['NB']*1e-6, c=self.df['V'], cmap=cmap)
        #plt.errorbar(self.df['W'], self.df['NB'], xerr=self.df['W_err'], yerr=self.df['NB_err'], fmt='none', ecolor='black')

        cbar = plt.colorbar()
        cbar.set_label('Reverse bias (V)')

        # Add labels and title
        plt.xlabel('Depth ($\mu$m)')
        plt.ylabel('$N_{B}$ ($cm^{-3}$)')
        plt.yscale('log')
        plt.title('Doping profile - '+f'{args.sensor}')

        plt.savefig(self.dprocol_outputPath)
        print('Plot saved in', self.dprocol_outputPath)
        plt.close()
        
    #####################
    # PRIVATE-LIKE METHODS

    def fit(self, graph, fit_min, fit_max, init_fit_pars=None, lim_fit_pars=None, line_color=None):
        '''
            Fit a linear function to the data and return the fitted function 

            Parameters
            ----------
            graph : TGraphErrors
                Graph with the data
            fit_min : float
                Minimum x value for the fit
            fit_max : float
                Maximum x value for the fit
            line_color : int
                Color of the fitted function    

            Returns
            -------
            fit : TF1
                Fitted function 
        ''' 

        # Fit the graph
        fit = TF1('fit', '[0] + [1]*x', fit_min, fit_max)
        if init_fit_pars is not None:   fit.SetParameters(init_fit_pars[0], init_fit_pars[1])
        if lim_fit_pars is not None:    
            fit.SetParLimits(0, lim_fit_pars[0][0], lim_fit_pars[0][1])
            fit.SetParLimits(1, lim_fit_pars[1][0], lim_fit_pars[1][1])
        fit.SetParNames('a', 'b')
        if line_color is not None:      fit.SetLineColor(line_color)
        fit.SetLineStyle(1)
        fit.SetLineWidth(2)
        fit.SetNpx(1000)
        graph.Fit(fit, 'rm+')
        print(f'Chi2/NDF: {fit.GetChisquare():.2f} / {fit.GetNDF():.2f}')
        return fit

    def find_intersection(self, fit1, fit2):
        '''
            Find the intersection of two linear functions

            Parameters
            ----------
            fit1 : TF1
                First linear function
            fit2 : TF1
                Second linear function

            Returns
            -------
            intersection : float
                Intersection point of the two linear functions
        '''
        # Set up a new function that is the difference between the two input functions
        diff_func = TF1("diff_func", "([0]+[1]*x)", fit1.GetXmin(), fit2.GetXmax())
        diff_func.SetParameters((fit1.GetParameter(0) - fit2.GetParameter(0)), (fit1.GetParameter(1) - fit2.GetParameter(1)))

        # Find the intersection point of the difference function
        intersection = diff_func.GetX(0)

        return intersection

    def evaluate_W_NB(self, Vbi=-0.6):
        '''
            Evaluate the depletion width and the doping concentration at the depletion voltage

            Parameters
            ----------
            Vbi (float) :   Built-in voltage (V)

            Returns
            -------
        '''

        eSi = 11.7 * 8.854e-12                                          # F/m - dielectric constant of silicon
        A = 1e-6                                                        # m^2 - detector area
        q = 1.602e-19                                                   # C - electron charge
            
        self.df['V_abs'] = np.abs(self.df['V'] + self.inversion*Vbi)    # V - absolute value of the bias voltage

        # doping concentration
        self.df['NB'] = 2 / (eSi * A*A * q * self.df['derivative'])     # m^-3
        self.df['NB_err'] = 2 / (eSi * A*A * q * self.df['derivative']**2) * self.df['derivative_err'] # m^-3

        # depletion width
        self.df['W'] = np.sqrt(2 * self.df['V_abs'] * eSi / (q * self.df['NB'])) # m
        self.df['W_err'] = np.sqrt((self.df['V_err']*self.df['W']/(2*self.df['V_abs']))**2 + (self.df['NB_err']*self.df['W']/(2*self.df['NB']))**2)
        





if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='depletion_voltage.py', description='Script to measure the depletion voltage from fits of a 1/C^2 vs V plot')
    parser.add_argument('--input', type=str, help='Input file with the 1/C^2 vs V data', required=True)
    parser.add_argument('--output', type=str, help='Output file with the plot', required=True)
    parser.add_argument('--config', type=str, default='Probe-station/src/depletion_voltage_conf.yml', help='Configuration file with the sensor parameters')
    parser.add_argument('--config_plot', type=str, default='Probe-station/src/plot_depletion_voltage_conf.yml', help='Configuration file with the sensor parameters')
    parser.add_argument('--sensor', type=str, help='Sensor name', required=True)
    parser.add_argument('--verbose', action='store_true', help='Verbose mode')
    args = parser.parse_args()

    df = pd.read_csv(args.input, comment='#')
    depletion_analysis = DepletionAnalysis(df=df, args=args)

    depletion_analysis.preprocess_data()
    depletion_analysis.inverseC2_vs_V()
    depletion_analysis.print_zoom(0., 33., -0.2e-4, 1.4e-4)

    depletion_analysis.derivative_plot()
    depletion_analysis.doping_concentration()
    depletion_analysis.doping_profile()
    depletion_analysis.depletion_depth()
    depletion_analysis.doping_profile_color()

    depletion_analysis.close()

    if args.verbose:    
        with open(os.path.splitext(args.output)[0] + '_verbose.txt', 'w') as f:
            f.write(depletion_analysis.df.to_string())
        print('Dataframe saved in', os.path.splitext(args.output)[0] + '_verbose.txt')
    