from ROOT import TFile, TCanvas, TPad

def main():
    
    inFilePIN = TFile('Probe-station/data/output/depletion_voltage_PiN.root')
    inFileStrip = TFile('Probe-station/data/output/depletion_voltage_strip.root')
    
    c1 = inFileStrip.Get('1C2')
    c2 = inFilePIN.Get('1C2')

    canvas = TCanvas('canvas', 'canvas', 600, 600)

    
    #pad1 = TPad("pad1", "pad1", 0.05, 0.05, 0.48, 0.95)
    #pad1.Draw()

    gr = c1.GetListOfPrimitives().FindObject('Graph')
    canvas.cd()
    gr.Draw('ap')

    primitives1 = c1.GetListOfPrimitives()
    for primitive in primitives1: 
        clone = primitive.Clone()      
        canvas.cd()
        if primitive.ClassName() == 'TGraph':   print()
        else:                                   clone.Draw('same')
    

    #canvas.cd(2)
    #primitives2 = c2.GetListOfPrimitives()
    #for primitive in primitives2: 
    #    clone = primitive.Clone()      
    #    if primitive.ClassName() == 'TGraph':   clone.Draw('same')    
    #    else:                                   clone.Draw('same')

    canvas.Modified()
    canvas.Update()

    #pad2 = TPad("pad2", "pad2", 0.55, 0.1, 0.95, 0.9)
    ##pad2.DrawFrame(0, 0, 1000, 1000)
    #
    #c2.Draw()
    #canvas.cd(2)
    #pad2.Draw()

    canvas.SaveAs('Probe-station/data/output/Figures/depletion_voltage_strip_pin.pdf')

if __name__ == '__main__':
    main()