import pandas as pd
import numpy as np

from ROOT import TGraphErrors, TCanvas, TF1, TLegend, gStyle, TText, TFile, TLatex
from ROOT import kGreen, kAzure

def calibrate(inputFile: str, outputPDF: str, output: TFile, name: str, color: int = kGreen, linecolor: int = kAzure):
    '''
        Calibrate the charge vs amplitude

        Parameters
        ----------
        inputFile (str) :   Path to the input csv file
        outputPDF (str) :   Path to the output pdf file
        ouput (TFile) :     Output ROOT file
        name (str) :        Name of the sensor
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
    df['Charge'] = 1000 * np.abs(df[f'Atot{name}'] - df[f'Abas{name}']) / impedance   # fC
    df['Charge_err'] = 1000 * np.sqrt((df[f'Atot_err{name}']/np.sqrt(1000))**2 + (df[f'Abas_err{name}']/np.sqrt(1000))**2) / impedance

    graph = TGraphErrors(len(df), np.asarray(df['Ampl_diode'], dtype=float), np.asarray(df['Charge'], dtype=float), 
                         np.asarray(df['Ampl_diode_err']/np.sqrt(1000), dtype=float), np.asarray(df['Charge_err'], dtype=float))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1)
    graph.SetMarkerColor(color)
    graph.SetTitle('Charge vs Signal Amplitude; Signal Amplitude (reference diode) (mV); Charge (fC)')

    # Fit the graph
    fit = TF1('fit', '[0] + [1] * x', 0, 300)
    fit.SetLineColor(linecolor)
    fit.SetLineStyle(2)
    fit.SetLineWidth(2)
    fit.SetParNames('a', 'b')
    fit.SetParameters(1, 1)
    if name == 'LGAD':  fit.SetParameters(-3e5, 14e3)
    graph.Fit(fit, 'RMF+')

    # Plot the graph
    gStyle.SetOptFit(0)
    gStyle.SetOptStat(0)
    canvas = TCanvas(f'{name}', 'canvas', 800, 600)
    graph.Draw('AP')
    fit.Draw('same')

    # Add the legend
    legend = TLegend(0.65, 0.35, 0.85, 0.5)
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
    latex.DrawLatex(0.65, 0.3, f'a = {fit.GetParameter(0):.2f} fC/mV')
    latex.DrawLatex(0.65, 0.25, f'b = {fit.GetParameter(1):.2f} fC')
    latex.DrawLatex(0.65, 0.2, '#chi^2 / NDF ='+f'{fit.GetChisquare():.0f} / {fit.GetNDF()}')
    
    # Save the plot
    canvas.SaveAs(outputPDF)
    output.cd()
    canvas.Write()

    df.to_csv(f'TCT/data/output/charge_vs_amplitude_{name}.csv', sep='\t', index=False)

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

    colors = [2, 3]
    linecolors = [4, 6]
    names = ['LGAD', 'PIN']
    inputFiles = ['TCT/data/input/Gain_vs_Bias_LGAD.csv', 'TCT/data/input/Gain_vs_Bias_Pin.csv']
    outputPDFs = ['TCT/data/output/Figures/charge_vs_amplitude_LGAD.pdf', 'TCT/data/output/Figures/charge_vs_amplitude_PIN.pdf']
    output = TFile('TCT/data/output/charge_vs_amplitude.root', 'RECREATE')

    graphs = []
    fits = []

    for name, inputFile, outputPDF, color, linecolor in zip(names, inputFiles, outputPDFs, colors, linecolors):
        graph, fit = calibrate(inputFile, outputPDF, output, name, color, linecolor)
        graphs.append(graph)
        fits.append(fit)

    canvas = TCanvas('charge_vs_amplitude_all', 'canvas', 800, 600)
    canvas.DrawFrame(25, 10, 30, 5e4, 'Charge vs Amplitude; Signal Amplitude (reference diode) (mV); Charge (fC)')
    canvas.SetLogy()

    legend = TLegend(0.65, 0.5, 0.9, 0.7)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.15, 0.9, 'Laser system calibration')

    for name, graph, fit in zip(names, graphs, fits):
        graph.Draw('P same')
        fit.Draw('same')

        legend.AddEntry(graph, f'{name} data', 'p')
        legend.AddEntry(fit, f'{name} fit', 'l')

    legend.Draw('same')    
    
    canvas.SaveAs('TCT/data/output/Figures/charge_vs_amplitude_all.pdf')
    output.cd()
    canvas.Write()  

    output.Close()

def main_5MIP():
    '''
        Main function, 5 MIPs

        Parameters
        ----------
        None

        Returns
        -------
        None
    '''

    colors = [2, 3]
    linecolors = [4, 6]
    names = ['LGAD', 'PIN']
    inputFiles = ['TCT/data/input/Gain_vs_Bias_LGAD_5MIP.csv', 'TCT/data/input/Gain_vs_Bias_Pin_5MIP.csv']
    outputPDFs = ['TCT/data/output/Figures/charge_vs_amplitude_LGAD_5MIP.pdf', 'TCT/data/output/Figures/charge_vs_amplitude_PIN_5MIP.pdf']
    output = TFile('TCT/data/output/charge_vs_amplitude_5MIP.root', 'RECREATE')

    graphs = []
    fits = []

    for name, inputFile, outputPDF, color, linecolor in zip(names, inputFiles, outputPDFs, colors, linecolors):
        graph, fit = calibrate(inputFile, outputPDF, output, name, color, linecolor)
        graphs.append(graph)
        fits.append(fit)

    canvas = TCanvas('charge_vs_amplitude_all', 'canvas', 800, 600)
    canvas.DrawFrame(25, 10, 30, 5e4, 'Charge vs Amplitude; Signal Amplitude (reference diode) (mV); Charge (fC)')
    canvas.SetLogy()

    legend = TLegend(0.65, 0.5, 0.9, 0.7)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.15, 0.9, 'Laser system calibration')

    for name, graph, fit in zip(names, graphs, fits):
        graph.Draw('P same')
        fit.Draw('same')

        legend.AddEntry(graph, f'{name} data', 'p')
        legend.AddEntry(fit, f'{name} fit', 'l')

    legend.Draw('same')    
    
    canvas.SaveAs('TCT/data/output/Figures/charge_vs_amplitude_all_5MIP.pdf')
    output.cd()
    canvas.Write()  

if __name__ == '__main__':
    main()
    main_5MIP()
