#!python3

'''
    Script to measure the depletion voltage from fits of a 1/C^2 vs V plot

'''

import argparse


import numpy as np
import pandas as pd

import uncertainties as unc

from ROOT import TFile, TGraphErrors, TF1, TCanvas, TLegend, gStyle, TText

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

def fit(graph: TGraphErrors, fit_min, fit_max, line_color=2):
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
    fit.SetParameters(0.0, 0.0)
    fit.SetParNames('a', 'b')
    fit.SetLineColor(line_color)
    fit.SetLineStyle(2)
    fit.SetLineWidth(2)
    fit.SetNpx(1000)
    graph.Fit(fit, 'rm+')

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
    diff_func = TF1("diff_func", "([0]-[1])", fit1.GetXmin(), fit1.GetXmax())
    diff_func.SetParameters(fit1.GetParameter(0), fit2.GetParameter(0))

    # Find the intersection point of the difference function
    intersection = diff_func.GetX(0, fit1.GetXmax())

    return intersection




def main():
    parser = argparse.ArgumentParser(prog='depletion_voltage.py',
                                     description='''Script to measure the depletion voltage from fits of a 1/C^2 vs V plot''')
    parser.add_argument('--input', type=str, help='Input file with the 1/C^2 vs V data', required=True)
    parser.add_argument('--output', type=str, help='Output file with the plot', required=True)
    
    args = parser.parse_args()

    # Read the data
    df = pd.read_csv(args.input)
    df = compute_1c2(df, 'C', 'C_err')
    print(df.describe())

    # Plot the data
    graph = TGraphErrors(len(df['V']), np.asarray(df['V'], dtype=float), np.asarray(df['1_C2'], dtype=float), np.asarray(df['V_err'], dtype=float), np.asarray(df['1_C2_err'], dtype=float))
    graph.SetTitle('#frac{1}{C^{2}} vs V LGAD; V [V]; #frac{1}{C^{2}} [pF^{-2}]')

    # Fit the data
    fit1 = fit(graph, 0.0, 100.0, 797)
    fit2 = fit(graph, 100.0, 200.0, 862)

    # Find the intersection point
    intersection = find_intersection(fit1, fit2)

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    canvas.SetGrid()
    graph.Draw('AP')
    fit1.Draw('same')
    fit2.Draw('same')

    legend = TLegend(0.15, 0.15, 0.35, 0.35)
    legend.SetBorderSize(0)
    legend.AddEntry(graph, 'Data', 'p')
    legend.AddEntry(fit1, 'Gain layer depletion', 'l')
    legend.AddEntry(fit2, 'Substrate depletion', 'l')
    legend.Draw()

    text = TText(0.15, 0.4, 'Depletion voltage: {:.2f} V'.format(intersection))

    outFile = TFile(args.output, 'recreate')
    canvas.Write()


    

    

if __name__ == '__main__':
    main()