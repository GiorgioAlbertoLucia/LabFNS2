import pandas as pd
import numpy as np

from ROOT import TGraphErrors, TCanvas, TF1, TLegend, gStyle, TText, TFile, TLatex
from ROOT import kGreen, kAzure

def calibrate(inputFile: str, outputPDF: str, output: TFile, color: int = kGreen, linecolor: int = kAzure):
    '''
        Calibrate the charge vs amplitude

        Parameters
        ----------
        inputFile (str) :   Path to the input csv file
        outputPDF (str) :   Path to the output pdf file
        ouput (TFile) :     Output ROOT file
        color (int) :       Color of the graph
        linecolor (int) :   Color of the fit

        Returns
        -------
        None
    '''

    # Read the csv file
    df = pd.read_csv(inputFile, comment='#')

    # Compute the charge
    impedance = 50 # Ohm
    amplification = 100
    pWb_to_fWb = 1000
    df['Charge'] = (np.abs(df[f'Asign']) - np.abs(df[f'Anoise'])) * pWb_to_fWb / (amplification * impedance)   # fC
    df['Charge_err'] = np.sqrt((df[f'Asign_err']/np.sqrt(1000))**2 + (df[f'Anoise_err']/np.sqrt(1000))**2) * pWb_to_fWb / (amplification * impedance)

    graph = TGraphErrors(len(df), np.asarray(df['Ampl_diode'], dtype=float), np.asarray(df['Charge'], dtype=float), 
                         np.asarray(df['Ampl_diode_err'], dtype=float), np.asarray(df['Charge_err'], dtype=float))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1)
    graph.SetMarkerColor(color)
    graph.SetTitle('Laser system calibration; Signal Amplitude (reference diode) (mV); Charge (fC)')

    # Fit the graph
    fit = TF1('fit', '[0] + [1] * x', 0, 300)
    fit.SetLineColor(linecolor)
    fit.SetLineStyle(2)
    fit.SetLineWidth(2)
    fit.SetParNames('a', 'b')
    fit.SetParameters(1, 1)
    graph.Fit(fit, 'RMF+')

    # Plot the graph
    gStyle.SetOptFit(0)
    gStyle.SetOptStat(0)
    canvas = TCanvas(f'diode_calibration', 'canvas', 800, 600)
    graph.Draw('AP')
    fit.Draw('same')

    # Add the legend
    legend = TLegend(0.55, 0.35, 0.75, 0.5)
    legend.AddEntry(graph, 'Data', 'p')
    legend.AddEntry(fit, 'Fit', 'l')
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.Draw()

    # Add the text
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.55, 0.3, f'a = ({fit.GetParameter(0):.2f} #pm {fit.GetParError(0):.2f}) fC')
    latex.DrawLatex(0.55, 0.25, f'b = ({fit.GetParameter(1):.3f} #pm {fit.GetParError(1):.3f}) fC'+' mV^{-1}')
    latex.DrawLatex(0.55, 0.2, '#chi^{2} / NDF ='+f'{fit.GetChisquare():.0f} / {fit.GetNDF()}')
    
    # Save the plot
    canvas.SaveAs(outputPDF)
    output.cd()
    canvas.Write()

    df.to_csv(f'TCT/data/output/charge_vs_amplitude.csv', sep='\t', index=False)

    return graph, fit

def main():
    '''
        Main function

        Parameters
        ----------
        None

        Returns
        -------
        None
    '''

    color = 2
    linecolor = 2
    
    inputFile = 'TCT/data/input/diode_calibration.csv'
    outputPDF = 'TCT/data/output/Figures/charge_vs_amplitude.pdf'
    output = TFile('TCT/data/output/charge_vs_amplitude.root', 'RECREATE')

    calibrate(inputFile, outputPDF, output, color, linecolor)

    output.Close()


if __name__ == '__main__':
    main()
