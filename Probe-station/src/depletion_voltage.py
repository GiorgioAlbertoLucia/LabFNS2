#!python3

'''
    Script to measure the depletion voltage from fits of a 1/C^2 vs V plot
    To run from the LabFNS2 folder:
        python3 -m Probe-station.src.depletion_voltage --input <input_file> --output <output_file>
        ./Probe-station/src/depletion_voltage.py --input INPUT_FILE --output OUTPUT_FILE

'''

import os
import argparse

import numpy as np
import pandas as pd

#import uncertainties as unc

from ROOT import TFile, TGraphErrors, TF1, TCanvas, TLegend, gStyle, TText
from ROOT import kGreen, kAzure

def compute_1c2(df, y: str, y_err: str):
    '''
        Compute 1/C^2 from the data
    
        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe with the data
        y : str 
            Name of the column with the capacitance
        y_err : str 
            Name of the column with the capacitance error
    
        Returns
        -------
        df : pandas.DataFrame
            Dataframe with the 1/C^2 and its error
    '''

    df['1_C2'] = 1.0 / (df[y] * df[y])
    df['1_C2_err'] = 2.0 * df[y_err] / (df[y] * df[y] * df[y])

    return df

def fit(graph: TGraphErrors, fit_min, fit_max, init_fit_pars=None, lim_fit_pars=None, line_color=None):
    '''
        Fit a linear function to the data

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

def find_intersection(fit1: TF1, fit2: TF1):
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




def main():
    parser = argparse.ArgumentParser(prog='depletion_voltage.py',
                                     description='''Script to measure the depletion voltage from fits of a 1/C^2 vs V plot''')
    parser.add_argument('--input', type=str, help='Input file with the 1/C^2 vs V data', required=True)
    parser.add_argument('--output', type=str, help='Output file with the plot', required=True)
    parser.add_argument('--sensor', type=str, help='Sensor name', required=True)
    
    args = parser.parse_args()

    # Read the data
    df = pd.read_csv(args.input, comment='#')
    df = compute_1c2(df, 'C', 'C_err')
    df['V_abs'] = np.abs(df['V'])

    # Plot the data
    graph = TGraphErrors(len(df['V_abs']), np.asarray(df['V_abs'], dtype=float), np.asarray(df['1_C2'], dtype=float), np.asarray(df['V_err'], dtype=float), np.asarray(df['1_C2_err'], dtype=float))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1)
    graph.SetMarkerColor(kAzure+1)
    graph.SetTitle('1/C^{2} vs |V| '+f'{args.sensor}'+'; Absolute reverse bias [V]; 1/C^{2} [pF^{-2}]')

    # Fit the data
    fit1 = fit(graph, 26.9, 45., init_fit_pars=[0.16, 0.0], line_color=797)
    fit2 = fit(graph, 24, 27, init_fit_pars=[0.0, 0.01], line_color=862)
    fit3 = fit(graph, 0., 24.1, init_fit_pars=[0., 0.0], #lim_fit_pars=[[-0.1, 0.1], [-0.03, -0.01]], 
               line_color=kGreen+2)    

    # Find the intersection point
    intersection1 = find_intersection(fit2, fit1)
    intersection2 = find_intersection(fit3, fit2)

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    canvas.SetGrid()
    graph.Draw('AP')
    fit1.Draw('same')
    fit2.Draw('same')
    fit3.Draw('same')
    

    legend = TLegend(0.15, 0.15, 0.35, 0.35)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.AddEntry(graph, 'Data', 'pe')
    legend.AddEntry(fit1, 'Gain layer depletion', 'l')
    legend.AddEntry(fit2, 'Substrate depletion', 'l')
    legend.AddEntry(fit3, 'Sensor depletion', 'l')
    legend.Draw()

    text = TText(0.15, 0.4, 'Depletion voltage of the sensor: {:.0f} V'.format(intersection2))
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.04)
    text.Draw()

    outFile = TFile(args.output, 'recreate')
    canvas.Write()

    output_path = os.path.splitext(args.output)[0] + '.pdf'
    canvas.SaveAs(os.path.join(os.path.join(os.path.dirname(output_path), 'Figures'), os.path.basename(output_path)))

    outFile.Close()


    

    

if __name__ == '__main__':
    main()