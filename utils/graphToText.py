from ROOT import TGraph, TFile

def graphToCsv(graph: TGraph, outFilePath: str):

    x = graph.GetX()
    y = graph.GetY()

    with open(outFilePath, 'w') as outFile:
        outFile.write("x,y\n")
        for i in range(graph.GetN()):
            outFile.write(f"{x[i]},{y[i]}\n")

if __name__ == '__main__':

    event = 1
    pixel = 9
    sampling = 25

    inFile = TFile.Open('ITS3/Data/run175174828_230428174901.root')
    graph = inFile.Get(f'grEv{event}Px{pixel}samp{sampling}')

    graphToCsv(graph, f'ITS3/Data/run175174828_230428174901_ev{event}_px{pixel}_samp{sampling}.csv')