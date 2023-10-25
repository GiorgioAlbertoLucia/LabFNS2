import pandas as pd
import numpy as np

from ROOT import TGraphErrors, TCanvas, TF1, TLegend, gStyle, TText
from ROOT import kGreen, kAzure

def main(inputFile: str, outputPDF: str):
    '''
        Main function

        Parameters
        ----------
        inputFile : str
            Path to the input csv file
        outputPDF : str
            Path to the output pdf file

        Returns
        -------
        None
    '''

    # Read the csv file
    df = pd.read_csv(inputFile, comment='#')

    # Compute the charge
    impedance = 50 # Ohm
    df['Charge'] = df['SignA'] / impedance
    df['Charge_err'] = df['SignA_err'] / impedance
    print(df)

    graph = TGraphErrors(len(df), np.asarray(df['Ampl'], dtype=float), np.asarray(df['Charge'], dtype=float), np.asarray(df['Ampl_err'], dtype=float), np.asarray(df['Charge_err'], dtype=float))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1)
    graph.SetMarkerColor(kGreen)
    graph.SetTitle('Charge vs Signal Amplitude; Signal Amplitude [mV]; Charge [fC]')

    # Fit the graph
    fit = TF1('fit', '[0] + [1] * x', 0, 300)
    fit.SetLineColor(kAzure)
    fit.SetLineStyle(2)
    fit.SetLineWidth(2)
    fit.SetParNames('a', 'b')
    fit.SetParameters(1, 1)
    graph.Fit(fit, 'RM+')

    # Plot the graph
    gStyle.SetOptFit(0)
    gStyle.SetOptStat(0)
    canvas = TCanvas('canvas', 'canvas', 800, 600)
    graph.Draw('AP')
    fit.Draw('same')

    # Add the legend
    legend = TLegend(0.15, 0.75, 0.35, 0.85)
    legend.AddEntry(graph, 'Data', 'p')
    legend.AddEntry(fit, 'Fit', 'l')
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.Draw()

    # Add the text
    text = TText(0.15, 0.65, f'a = {fit.GetParameter(0):.2f} fC/mV')
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.04)
    text.Draw()
    text2 = TText(0.15, 0.6, f'b = {fit.GetParameter(1):.2f} fC')
    text2.SetNDC()
    text2.SetTextFont(42)
    text2.SetTextSize(0.04)
    text2.Draw()    
    text3 = TText(0.15, 0.55, '#chi^2 / NDF ='+f'{fit.GetChisquare():.0f} / {fit.GetNDF()}')
    text3.SetNDC()
    text3.SetTextFont(42)
    text3.SetTextSize(0.04)
    text3.Draw()    
    
    # Save the plot
    canvas.SaveAs(outputPDF)


if __name__ == '__main__':
    main(inputFile='TCT/data/input/Charge_vs_Amplitude.csv', outputPDF='TCT/data/output/Figures/charge_vs_amplitude.pdf')
