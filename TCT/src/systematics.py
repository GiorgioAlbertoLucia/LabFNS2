
import pandas as pd
import numpy as np
from ROOT import TFile, TCanvas, TLegend, TGraphErrors, TLatex


def plot_systematics(inputFile, outputFile, color, name):

    df = pd.read_csv(inputFile, comment='#')

    graph = TGraphErrors(len(df['Ampl_diode']), 
                 np.asarray(df['Ampl_diode'], dtype=np.float64),
                 np.asarray(df['Bias'], dtype=np.float64),
                 np.asarray(df['Ampl_diode_err'], dtype=np.float64),
                 np.asarray(df['Bias_err'], dtype=np.float64))
    
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(color)
    graph.SetTitle('Systematics - ' + name + '; Signal amplitude (reference diode) (mV); Bias voltage (' + name +') (V)')

    outputFile.cd()
    graph.Write('systematics_'+name)

    return graph


def main():

    inputFiles = ['TCT/data/input/Gain_vs_Bias_LGAD.csv', 'TCT/data/input/Gain_vs_Bias_Pin.csv']
    names = ['LGAD', 'PIN']
    colors = [797, 857]
    outputFile = TFile('TCT/data/output/systematics.root', 'RECREATE')

    graphs = []
    for i in range(len(inputFiles)):    graphs.append(plot_systematics(inputFiles[i], outputFile, colors[i], names[i]))

    canvas = TCanvas('canvas', 'canvas', 800, 600)
    canvas.Divide(2, 1)
    canvas.SetLeftMargin(0.15)

    for i in range(len(graphs)):
        
        canvas.cd(i+1)
        graphs[i].Draw('AP')
        graphs[i].GetYaxis().SetTitleOffset(1.5)

        
        
                #latex = TLatex()
                #latex.SetTextSize(0.04)
                #latex.SetTextAlign(13)
                #latex.DrawLatexNDC(0.15, 0.85, 'Systematics in gain vs bias')
                #latex.DrawLatexNDC(0.14, 0.81, 'acquisition for the '+names[i]+' sensor')
                #latex.DrawLatexNDC(0.15, 0.77, 'Laser attenuation: 1.0')

    canvas.SaveAs('TCT/data/output/Figures/systematics.pdf')
    outputFile.cd()
    canvas.Write('systematics')
    outputFile.Close()

if __name__ == '__main__':
    main()
